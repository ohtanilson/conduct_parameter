
using LinearAlgebra, Distributions
using Statistics, Random, MultivariateStats
using JuMP, Ipopt
using DelimitedFiles, JLD, CSV, DataFrames
using Plots, Combinatorics, Dates, StatsPlots
using Parameters: @unpack, @with_kw


#---------------------------------------------------------------------------------------------

function termination_status_code(status)

    status_code = [
    OPTIMAL,
    INFEASIBLE,  
    DUAL_INFEASIBLE,  
    LOCALLY_SOLVED,  
    LOCALLY_INFEASIBLE,
    INFEASIBLE_OR_UNBOUNDED,  
    ALMOST_OPTIMAL,
    ALMOST_INFEASIBLE,
    ALMOST_DUAL_INFEASIBLE,
    ALMOST_LOCALLY_SOLVED, 
    ITERATION_LIMIT,  
    TIME_LIMIT,  
    NODE_LIMIT,  
    SOLUTION_LIMIT,  
    MEMORY_LIMIT,
    OBJECTIVE_LIMIT,  
    NORM_LIMIT,
    OTHER_LIMIT,
    SLOW_PROGRESS,
    NUMERICAL_ERROR,
    INVALID_MODEL,
    INVALID_OPTION,
    INTERRUPTED,
    OTHER_ERROR]

    for i = 1:length(status_code)
        if (status == status_code[i]) == true
            return i
            break
        end
    end
end


market_parameters_log = @with_kw (
    α_0 = 10, # Demand parameter
    α_1 = 1,
    α_2 = 1,
    α_3 = 1,
    γ_0 = 1,  # Marginal cost parameter
    γ_1 = 1,
    γ_2 = 1,
    γ_3 = 1,
    θ = 0.5,  # Conduct paramter
    σ = 1,    # Standard deviation of the error term
    T = 50,   # Number of markets
    S = 1000, # Number of simulation
)

# The following functions are for implementing the GMM estimation to the simultanoues equations.

function GMM_estimation_simultaneous(T, P, Z, X, X_s, X_d, Ω)
    
    L = size(Z, 2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "tol", 1e-15)
    set_optimizer_attribute(model, "max_iter", 1000)
    set_optimizer_attribute(model, "acceptable_tol", 1e-12)
    set_silent(model)
    @variable(model, β[k = 1:K_d+K_s-1])
    @variable(model, 0 <= θ <= 1)

    r = Any[];
    g = Any[];
    for t =1:T
        push!(r, @NLexpression(model, P[t]- sum(β[k] * X[2*t-1,k] for k = 1:K_d) ))
        push!(r, @NLexpression(model, P[t]- sum(β[k] * X[2*t,k] for k = K_d+1:K_d+K_s-1) + log(1 - θ *(β[2] + β[3] * X[2*t, end]) ) ) )
    end

    for l = 1:L
        push!(g, @NLexpression(model, sum(Z[t,l] * r[t] for t = 1:2*T)))
    end

    @NLobjective(model, Min, sum( g[l] *Ω[l,k] * g[k] for l = 1:L, k = 1:L))

    optimize!(model)
    
    α_hat = value.(β)[1:K_d]
    γ_hat = value.(β)[K_d+1:end]
    θ_hat = value.(θ)

    #=
    if termination_status(model) == OPTIMAL

        α_hat = value.(β)[1:K_d]
        γ_hat = value.(β)[K_d+1:end]
        θ_hat = value.(θ)

        return α_hat, γ_hat, θ_hat
    else
        error("The model was not solved correctly.")
    end
    =#

    return α_hat, γ_hat, θ_hat, termination_status_code(termination_status(model))
end


function estimation_nonlinear_2SLS_simultaneous(parameter, data)
    @unpack T = parameter

    Q  = data.Q
    w  = data.w
    r  = data.r
    z  = data.z
    iv_w = data.iv_w
    iv_r = data.iv_r
    p  = data.P
    y  = data.y

    iv = hcat(iv_w, iv_r)

    X_d = []
    X_s = []
    X   = []
    Z   = []
    Z_d = []
    Z_s = []
    P   = []

    for t = 1:T

        Z_dt = vcat(1, z[t], iv[t,:], y[t])
        Z_st = vcat(1, z[t], w[t], r[t], y[t])
        Z_t = [Z_dt zeros(length(Z_dt));  zeros(length(Z_st)) Z_st]'

        X_dt = vcat(1, -Q[t], -z[t].*Q[t], y[t])
        X_st = vcat(1, Q[t], w[t], r[t], z[t])
        X_t  = [X_dt zeros(length(X_dt));  zeros(length(X_st)) X_st]'

        push!(P, p[t])
        push!(X_d, X_dt')
        push!(X_s, X_st')
        push!(X, X_t)
        push!(Z, Z_t)
        push!(Z_d, Z_dt')
        push!(Z_s, Z_st')
    end

    Z   = reduce(vcat,(Z))
    X   = reduce(vcat,(X))
    X_d = reduce(vcat,(X_d))
    X_s = reduce(vcat,(X_s))
    Z_d = reduce(vcat,(Z_d))
    Z_s = reduce(vcat,(Z_s))

    Ω = inv(Z'Z/T)

    α_hat, γ_hat, θ_hat, status = GMM_estimation_simultaneous(T, P, Z, X, X_s, X_d, Ω)

    return α_hat, γ_hat, θ_hat, status
end


function simulation_nonlinear_2SLS_simultaneous(parameter, data) 

    @unpack T, S, σ = parameter 

    α_est = Vector{Float64}[]
    γ_est = Vector{Float64}[]
    θ_est = Float64[]
    status = Int64[]

    for s = 1:S
        data_s = data[(s-1)*T+1:s*T,:]
        α_est_s, γ_est_s, θ_est_s, status_s = estimation_nonlinear_2SLS_simultaneous(parameter, data_s)
    
        push!(α_est, α_est_s)
        push!(γ_est, γ_est_s)
        push!(θ_est, θ_est_s)
        push!(status,status_s)
    end
    
    α_est = reduce(vcat, α_est')
    γ_est = reduce(vcat, γ_est')
    θ_est = reduce(vcat, θ_est')
    status = reduce(vcat, status)

    estimation_result = DataFrame(
    T = T,
    σ = σ,
    α_0 = α_est[:,1],
    α_1 = α_est[:,2],
    α_2 = α_est[:,3],
    α_3 = α_est[:,4],
    γ_0 = γ_est[:,1],
    γ_1 = γ_est[:,2],
    γ_2 = γ_est[:,3],
    γ_3 = γ_est[:,4],
    θ = θ_est,
    status = status)

    return estimation_result
end

#---------------------------------------------------------------------------------------------

function GMM_estimation_separate(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d)
    

    QZ_hat = Z_d * inv(Z_d' * Z_d) * Z_d' * (Z_d[:,2] .* Q)
    Q_hat = Z_d * inv(Z_d' * Z_d) * Z_d' *  Q
    X_dd = hcat(ones(T), -Q_hat, -QZ_hat, X_d[:,end])
    α_hat = inv(X_dd' * X_dd) * (X_dd' * P)


    L = size(Z, 2)
    L_d = size(Z_d,2)
    L_s = size(Z_s,2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    Ω = inv(Z_s' * Z_s)/T

    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "tol", 1e-15)
    set_optimizer_attribute(model, "max_iter", 1000)
    set_optimizer_attribute(model, "acceptable_tol", 1e-12)
    set_silent(model)
    @variable(model, γ[k = 1:K_s-1])
    @variable(model, 0 <= θ <= 1)

    r = Any[];
    g = Any[];
    for t =1:T
        push!(r, @NLexpression(model, P[t]- sum(γ[k] * X_s[t,k] for k = 1:K_s-1) + log(1 - θ *(α_hat[2] + α_hat[3] * X_s[t, end]) ) ) )
    end

    for l = 1:L_s
        push!(g, @NLexpression(model, sum(Z_s[t,l] * r[t] for t = 1:T)))
    end

    @NLobjective(model, Min, sum( g[l] *Ω[l,k] * g[k] for l = 1:L_s, k = 1:L_s))

    optimize!(model)
    
    γ_hat = value.(γ)
    θ_hat = value.(θ)

    return α_hat, γ_hat, θ_hat, termination_status_code(termination_status(model))
end


function estimation_nonlinear_2SLS_separate(parameter, data)
    @unpack T = parameter

    Q  = data.Q
    w  = data.w
    r  = data.r
    z  = data.z
    iv_w = data.iv_w
    iv_r = data.iv_r
    p  = data.P
    y  = data.y

    iv = hcat(iv_w, iv_r)

    X_d = []
    X_s = []
    X   = []
    Z   = []
    Z_d = []
    Z_s = []
    P   = []

    for t = 1:T

        Z_dt = vcat(1, z[t], iv[t,:], y[t])
        Z_st = vcat(1, z[t], w[t], r[t], y[t])
        Z_t = [Z_dt zeros(length(Z_dt));  zeros(length(Z_st)) Z_st]'

        X_dt = vcat(1, -Q[t], -z[t].*Q[t], y[t])
        X_st = vcat(1, Q[t], w[t], r[t], z[t])
        X_t  = [X_dt zeros(length(X_dt));  zeros(length(X_st)) X_st]'

        push!(P, p[t])
        push!(X_d, X_dt')
        push!(X_s, X_st')
        push!(X, X_t)
        push!(Z, Z_t)
        push!(Z_d, Z_dt')
        push!(Z_s, Z_st')
    end

    Z   = reduce(vcat,(Z))
    X   = reduce(vcat,(X))
    X_d = reduce(vcat,(X_d))
    X_s = reduce(vcat,(X_s))
    Z_d = reduce(vcat,(Z_d))
    Z_s = reduce(vcat,(Z_s))


    α_hat, γ_hat, θ_hat, status = GMM_estimation_separate(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d)

    return α_hat, γ_hat, θ_hat, status
end


function simulation_nonlinear_2SLS_separate(parameter, data)

    @unpack T, S, σ = parameter 

    α_est = Vector{Float64}[]
    γ_est = Vector{Float64}[]
    θ_est = Float64[]
    status = Int64[]

    for s = 1:S
        data_s = data[(s-1)*T+1:s*T,:]
        α_est_s, γ_est_s, θ_est_s, status_s = estimation_nonlinear_2SLS_separate(parameter, data_s)
    
        push!(α_est, α_est_s)
        push!(γ_est, γ_est_s)
        push!(θ_est, θ_est_s)
        push!(status,status_s)
    end
    
    α_est = reduce(vcat, α_est')
    γ_est = reduce(vcat, γ_est')
    θ_est = reduce(vcat, θ_est')
    status = reduce(vcat, status)

    estimation_result = DataFrame(
    T = T,
    σ = σ,
    α_0 = α_est[:,1],
    α_1 = α_est[:,2],
    α_2 = α_est[:,3],
    α_3 = α_est[:,4],
    γ_0 = γ_est[:,1],
    γ_1 = γ_est[:,2],
    γ_2 = γ_est[:,3],
    γ_3 = γ_est[:,4],
    θ = θ_est,
    status = status)

    return estimation_result

end




#---------------------------------------------------------------------------------------------

# For testing the functions

parameter = market_parameters_log();
data = simulation_data_log(parameter);

@time simultaneous_test = simulation_nonlinear_2SLS_simultaneous(parameter, data)

@time separate_test = simulation_nonlinear_2SLS_separate(parameter, data)

describe(simultaneous_test)
describe(separate_test)



#---------------------------------------------------------------------------------------------------------
# Estimate the parameters for each number of markets and the value of the standard deviation of the error terms

for t = [50, 100, 200, 1000], sigma =  [0.001, 0.5, 1, 2]

    #=
    # Load the simulation data from the rds files
    filename_begin = "../conduct_parameter/output/data_linear_linear_n_"
    filename_end   = ".rds"

    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end

    filename = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end

    data = load(filename)
    data = DataFrames.sort(data, [:group_id_k])

    =#

    # Load the simulation data from the csv files
    filename_begin = "../conduct_parameter/output/data_loglinear_loglinear_n_"
    filename_end   = "_with_demand_shifter_y.csv"
    file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end
    data = DataFrame(CSV.File(file_name))

    # Set parameter values
    parameter = market_parameters(T = t, σ = sigma)
    
    # Estimation based on 2SLS
    estimation_result = simulation_nonlinear_2SLS_separate(parameter, data)

    # Uncomment below for checking the estimation result by the system of 2SLS
    # estimation_result = simulation_nonlinear_2SLS(parameter, data)

    # Save the estimation result as csv file. The file is saved at "output" folder
    filename_begin = "../conduct_parameter/output/parameter_hat_table_loglinear_loglinear_n_"
    filename_end   = "_with_demand_shifter_y.csv"
    file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end

    CSV.write(file_name, estimation_result)
end






#-----------------------------------------------------------------------------------------------------------------------



function GMM_function_separate(parameter, data)
    

    @unpack α_0, α_1, α_2, α_3,γ_0 , γ_1 ,γ_2 ,γ_3, θ,σ ,T = parameter

    γ = [γ_0 , γ_1 ,γ_2 ,γ_3]
    α = [α_0, α_1, α_2, α_3]

    Q  = data.Q
    w  = data.w
    r  = data.r
    z  = data.z
    iv_w = data.iv_w
    iv_r = data.iv_r
    p  = data.P
    y  = data.y

    iv = hcat(iv_w, iv_r)

    X_d = []
    X_s = []
    X   = []
    Z   = []
    Z_d = []
    Z_s = []
    P   = []

    for t = 1:T

        Z_dt = vcat(1, z[t], iv[t,:], y[t])
        Z_st = vcat(1, z[t], w[t], r[t], y[t])
        Z_t = [Z_dt zeros(length(Z_dt));  zeros(length(Z_st)) Z_st]'

        X_dt = vcat(1, -Q[t], -z[t].*Q[t], y[t])
        X_st = vcat(1, Q[t], w[t], r[t], z[t])
        X_t  = [X_dt zeros(length(X_dt));  zeros(length(X_st)) X_st]'

        push!(P, p[t])
        push!(X_d, X_dt')
        push!(X_s, X_st')
        push!(X, X_t)
        push!(Z, Z_t)
        push!(Z_d, Z_dt')
        push!(Z_s, Z_st')
    end

    Z   = reduce(vcat,(Z))
    X   = reduce(vcat,(X))
    X_d = reduce(vcat,(X_d))
    X_s = reduce(vcat,(X_s))
    Z_d = reduce(vcat,(Z_d))
    Z_s = reduce(vcat,(Z_s))

    L = size(Z, 2)
    L_d = size(Z_d,2)
    L_s = size(Z_s,2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    Ω = inv(Z_s' * Z_s)/T

    GMM_value = []


    for θ_hat = [0:0.01:0.501;], γ_hat = [-10:0.01:10;]

        r = P .- γ_hat .-sum(γ[k] .* X_s[:,k] for k = 2:K_s-1) + log.(1 .- θ_hat .*(α_1 .+ α_2 .* X_s[:, end]) ) 

        g = sum(Z_s[t,:] .* r[t] for t = 1:T)

        gmm_value = sum( g[l] *Ω[l,k] * g[k] for l = 1:L_s, k = 1:L_s)

        push!(GMM_value, gmm_value)
    end


#=
    for θ_hat = [0:0.01:1;]

        r = P .-sum(γ[k] .* X_s[:,k] for k = 1:K_s-1) + log.(1 .- θ_hat .*(α[2] .+ α[3] .* X_s[:, end]) ) 

        g = sum(Z_s[t,:] .* r[t] for t = 1:T)

        gmm_value = sum( g[l] *Ω[l,k] * g[k] for l = 1:L_s, k = 1:L_s)

        push!(GMM_value, gmm_value)
    end

    for γ_hat = [0:0.01:2;]

        r = P .- γ_hat .-sum(γ[k] .* X_s[:,k] for k = 2:K_s-1) + log.(1 .- θ.*(α[2] .+ α[3] .* X_s[:, end]) ) 

        g = sum(Z_s[t,:] .* r[t] for t = 1:T)

        gmm_value = sum( g[l] *Ω[l,k] * g[k] for l = 1:L_s, k = 1:L_s)

        push!(GMM_value, gmm_value)
    end

=#

    return reshape(GMM_value, (length([-10:0.01:10;]), length([0:0.01:0.501;])))
    #return GMM_value
end



data_s = data[1:50,:]

test = estimation_nonlinear_2SLS_simultaneous(parameter, data_s)

test = GMM_function_separate(parameter, data);
test_plot = plot(contour([0:0.01:0.501;], [-10:0.01:10;],test,
    xlabel="θ", ylabel="γ_0",
    title="Value of GMM, N =200, σ = 1"))
    vline!([0.5], linestyle=:dash)
    hline!([1], linestyle=:dash)










