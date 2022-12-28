using LinearAlgebra, Distributions
using Statistics, Random, MultivariateStats
using JuMP, Ipopt
using DelimitedFiles, JLD, CSV, DataFrames
using Plots, Combinatorics, Dates, StatsPlots
using Parameters: @unpack, @with_kw


#---------------------------------------------------------------------------------------------

function termination_status_code(status)

    status_code = [
    OPTIMAL,                   # 1
    INFEASIBLE,                # 2
    DUAL_INFEASIBLE,           # 3
    LOCALLY_SOLVED,            # 4
    LOCALLY_INFEASIBLE,        # 5
    INFEASIBLE_OR_UNBOUNDED,   # 6
    ALMOST_OPTIMAL,            # 7
    ALMOST_INFEASIBLE,         # 8
    ALMOST_DUAL_INFEASIBLE,    # 9
    ALMOST_LOCALLY_SOLVED,     # 10
    ITERATION_LIMIT,           # 11
    TIME_LIMIT,                # 12
    NODE_LIMIT,                # 13
    SOLUTION_LIMIT,            # 14
    MEMORY_LIMIT,              # 15
    OBJECTIVE_LIMIT,           # 16
    NORM_LIMIT,                # 17
    OTHER_LIMIT,               # 18
    SLOW_PROGRESS,             # 19
    NUMERICAL_ERROR,           # 20
    INVALID_MODEL,             # 21
    INVALID_OPTION,            # 22
    INTERRUPTED,               # 23
    OTHER_ERROR                # 24
    ]

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
    α_2 = 0.1,
    α_3 = 1,
    γ_0 = 1,  # Marginal cost parameter
    γ_1 = 1,
    γ_2 = 1,
    γ_3 = 1,
    θ = 0.3,  # Conduct paramter
    σ = 1,    # Standard deviation of the error term
    T = 50,   # Number of markets
    S = 1000, # Number of simulation
)

# The following functions are for implementing the GMM estimation to the simultanoues equations.

#---------------------------------------------------------------------------------------------------------

function GMM_estimation_simultaneous(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, parameter, estimation_method)
    
    L = size(Z, 2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    Ω = inv(Z' * Z)/T

    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "tol", 1e-15)
    set_optimizer_attribute(model, "max_iter", 1000)
    set_optimizer_attribute(model, "acceptable_tol", 1e-12)
    set_silent(model)
    @variable(model, β[k = 1:K_d+K_s-1])
    @variable(model, θ)

    r = Any[];
    g = Any[];
    for t =1:T
        push!(r, @NLexpression(model, P[t]- sum(β[k] * X[2*t-1,k] for k = 1:K_d) ))
        push!(r, @NLexpression(model, P[t]- sum(β[k] * X[2*t,k] for k = K_d+1:K_d+K_s-1) + log(1 - θ *(β[2] + β[3] * X[2*t, end]) ) ) )
    end

    for l = 1:L
        push!(g, @NLexpression(model, sum(Z[t,l] * r[t] for t = 1:2*T)))
    end

    if estimation_method[2] == :constraint
        for t = 1:T
            @NLconstraint(model, 0 <= 1 - θ *(β[2] + β[3] * X[2*t, end]))
        end
    end

    @NLobjective(model, Min, sum( g[l] *Ω[l,k] * g[k] for l = 1:L, k = 1:L))

    optimize!(model)
    
    α_hat = value.(β)[1:K_d]
    γ_hat = value.(β)[K_d+1:end]
    θ_hat = value.(θ)

    if  sum(1 .- θ_hat .*(α_hat[2] .+ α_hat[3] .* X_s[:,end]) .<= 0) == 0

        return α_hat, γ_hat, θ_hat, termination_status_code(termination_status(model))
    else 
        error("The estimation result violates the model assumption ")

    end
    
end


#---------------------------------------------------------------------------------------------

function GMM_estimation_separate(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, parameter, estimation_method)
    
    @unpack θ = parameter

    QZ_hat = Z_d * inv(Z_d' * Z_d) * Z_d' * (Z_d[:,2] .* Q)
    Q_hat = Z_d * inv(Z_d' * Z_d) * Z_d' *  Q
    X_dd = hcat(ones(T), -Q_hat, -QZ_hat, X_d[:,end])
    α_hat = inv(X_dd' * X_dd) * (X_dd' * P)

    L = size(Z, 2)
    L_d = size(Z_d,2)
    L_s = size(Z_s,2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    sample_violatiton_index = Int64[]

    # Pick up the index of market under which the inside of the log has a negative value
    for t = 1:T
        if 1 .- θ .*(α_hat[2] .+ α_hat[3] .* X_s[t,end]) <= 0
            push!(sample_violatiton_index, t)
        end
    end

    # If all markets do not satisfy the assumption, stop the supply estimation and rerutns missing values
    if length(sample_violatiton_index) == T

        γ_hat = repeat([missing], K_s-1)
        θ_hat = missing

        return α_hat, γ_hat, θ_hat, missing

    else
        
        # Drop the samples that violate the assumption
        if  1 <= length(sample_violatiton_index)

            sample_index = setdiff([1:T;], sample_violatiton_index)

            Z_s = Z_s[sample_index, :]
            X_s = X_s[sample_index, :]
            T = length(sample_index)
            L_s = size(Z_s,2)
            K_s = size(X_s, 2)
        end

        #Check if the weight matrix can be obtained
        if rank(Z_s' * Z_s) == L_s
            
            Ω = inv(Z_s' * Z_s)/T

            model = Model(Ipopt.Optimizer)
            set_optimizer_attribute(model, "tol", 1e-15)
            set_optimizer_attribute(model, "max_iter", 1000)
            set_optimizer_attribute(model, "acceptable_tol", 1e-12)
            set_silent(model)
            @variable(model, γ[k = 1:K_s-1])
            @variable(model, θ)

            r = Any[];
            g = Any[];
            for t =1:T
                push!(r, @NLexpression(model, P[t]- sum(γ[k] * X_s[t,k] for k = 1:K_s-1) + log(1 - θ *(α_hat[2] + α_hat[3] * X_s[t, end]) ) ) )
            end

            for l = 1:L_s
                push!(g, @NLexpression(model, sum(Z_s[t,l] * r[t] for t = 1:T)))
            end

            if estimation_method[2] == :constraint

                for t = 1:T
                    @NLconstraint(model, 0 <= 1 - θ *(α_hat[2] + α_hat[3] * X_s[t, end]) )
                end
            end

            @NLobjective(model, Min, sum( g[l] *Ω[l,k] * g[k] for l = 1:L_s, k = 1:L_s))

            optimize!(model)
            
            γ_hat = value.(γ)
            θ_hat = value.(θ)
        else

            # If the weight matrix can not be obtained, return missing values
            γ_hat = repeat([missing], K_s-1)
            θ_hat = missing

            return α_hat, γ_hat, θ_hat, missing
        end

        # Check if the supply estimation result satisfies the assumption
        if  sum(1 .- θ_hat .*(α_hat[2] .+ α_hat[3] .* X_s[:,end]) .<= 0) == 0

            return α_hat, γ_hat, θ_hat, termination_status_code(termination_status(model))
        else 
            error("The estimation result violates the model assumption ")
        end
    end

end


function estimation_nonlinear_2SLS(parameter, data, estimation_method)
    @unpack T = parameter

    Q  = data.logQ
    w  = data.w
    r  = data.r
    z  = data.z
    iv_w = data.iv_w
    iv_r = data.iv_r
    p  = data.logP
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
        Z_st = vcat(1, z[t], log(w[t]), log(r[t]), y[t])
        Z_t = [Z_dt zeros(length(Z_dt));  zeros(length(Z_st)) Z_st]'

        X_dt = vcat(1, -Q[t], -z[t].*Q[t], y[t])
        X_st = vcat(1, Q[t], log(w[t]), log(r[t]), z[t])
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

    if estimation_method[1] == :separate
        α_hat, γ_hat, θ_hat, status = GMM_estimation_separate(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, parameter, estimation_method)
    else
        α_hat, γ_hat, θ_hat, status = GMM_estimation_simultaneous(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, parameter, estimation_method)
    end

    return α_hat, γ_hat, θ_hat, status
end


function simulation_nonlinear_2SLS(parameter, data, estimation_method::Tuple{Symbol, Symbol})

    @unpack T, S, σ = parameter 

    α_est  = Vector{Float64}[]
    γ_est  = Vector{Union{Missing, Float64}}[]
    θ_est  = Union{Missing, Float64}[]
    status = Union{Missing, Int64}[]

    for s = 1:S
        data_s = data[(s-1)*T+1:s*T,:]
        α_est_s, γ_est_s, θ_est_s, status_s = estimation_nonlinear_2SLS(parameter, data_s, estimation_method)
    
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


#=

parameter = market_parameters_log(T = 1000, σ = 1);
data = load("../conduct_parameter/output/data_loglinear_loglinear_n_1000_sigma_1.rds");

describe(data)

@time simultaneous_test = simulation_nonlinear_2SLS(parameter, data, (:simultaneous, :non_constraint));
describe(simultaneous_test)
@time separate_test = simulation_nonlinear_2SLS(parameter, data, (:separate, :non_constraint));
describe(separate_test)

@time simultaneous_constraint_test = simulation_nonlinear_2SLS(parameter, data,  (:simultaneous, :constraint));
describe(simultaneous_constraint_test)
@time separate_constraint_test = simulation_nonlinear_2SLS(parameter, data, (:separate, :constraint));
describe(separate_constraint_test)


=#


#--------------------------------------------------------------------------------------------------------------
#estimation_methods = [(:separate,:non_constraint), (:separate,:constraint), (:simultaneous,:non_constraint), (:simultaneous,:constraint)];
estimation_methods = [(:separate,:non_constraint), (:separate,:constraint)];
#estimation_methods = [(:simultaneous,:non_constraint), (:simultaneous,:constraint)];


# Estimate the parameters for each number of markets and the value of the standard deviation of the error terms
for estimation_method = estimation_methods
    for t = [50, 100, 200, 1000], sigma =  [0.001, 0.5, 1, 2]

        # Load the simulation data from the rds files
        filename_begin = "../conduct_parameter/output/data_loglinear_loglinear_n_"
        filename_end   = ".rds"

        if sigma == 1 || sigma == 2
            sigma = Int64(sigma)
        end

        filename = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end

        data = load(filename)
        data = DataFrames.sort(data, [:group_id_k])

        #Uncomment the following lines to load the simulation data from the csv files
            #filename_begin = "../conduct_parameter/output/data_loglinear_loglinear_n_"
            #filename_end   = ".csv"
            #file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end
            #data = DataFrame(CSV.File(file_name))

        # Set parameter values
        parameter = market_parameters_log(T = t, σ = sigma)
        
        # Estimation based on 2SLS
        @time estimation_result = simulation_nonlinear_2SLS(parameter, data, estimation_method)

        # Save the estimation result as csv file. The file is saved at "output" folder
        filename_estimation = "_"*String(estimation_method[1])*"_"*String(estimation_method[2])

        filename_begin = "../conduct_parameter/output/parameter_hat_table_loglinear_loglinear_n_"
        filename_end   = ".csv"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end

        CSV.write(file_name, estimation_result, transform=(col, val) -> something(val, missing))

        println("")

    end
end

# Check the summary of the eatimation result
for estimation_method = estimation_methods
    
    for t = [50, 100, 200, 1000], sigma =  [0.001, 0.5, 1, 2]

        # Load the simulation data from the rds files
        filename_begin = "../conduct_parameter/output/data_loglinear_loglinear_n_"
        filename_end   = ".rds"

        if sigma == 1 || sigma == 2
            sigma = Int64(sigma)
        end

        filename = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end

        data = load(filename)
        data = DataFrames.sort(data, [:group_id_k])



        filename_estimation = "_"*String(estimation_method[1])*"_"*String(estimation_method[2])
        filename_begin = "../conduct_parameter/output/parameter_hat_table_loglinear_loglinear_n_"
        filename_end   = ".csv"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end
        estimation_result = DataFrame(CSV.File(file_name))

        parameter = market_parameters_log(T = t, σ = sigma)

        @unpack θ, S = parameter

        rate_satisfy_assumption = 1 - sum(sum(1 .- θ .*(estimation_result.α_1[s] .+ estimation_result.α_2[s] .* data.z[(s-1)*t+1:s*t,:]) .<= 0)  for s = 1:S)/(S*t)

        #@show rate_satisfy_assumption

        display(describe(estimation_result))
    end


    println("------------------------------------------------------------------------------------\n")
end


#-----------------------------------------------------------------------------------------------------------------------


#=
function GMM_function_separate(parameter, data)
    

    @unpack α_0, α_1, α_2, α_3,γ_0 , γ_1 ,γ_2 ,γ_3, θ,σ ,T = parameter

    γ = [γ_0 , γ_1 ,γ_2 ,γ_3]
    α = [α_0, α_1, α_2, α_3]

    Q  = data.logQ
    w  = data.w
    r  = data.r
    z  = data.z
    iv_w = data.iv_w
    iv_r = data.iv_r
    p  = data.logP
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
        Z_st = vcat(1, z[t], log(w[t]), log(r[t]), y[t])
        Z_t = [Z_dt zeros(length(Z_dt));  zeros(length(Z_st)) Z_st]'

        X_dt = vcat(1, -Q[t], -z[t].*Q[t], y[t])
        X_st = vcat(1, Q[t], log(w[t]), log(r[t]), z[t])
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


    for θ_hat = [0:0.01:0.6;], γ_hat = [-10:0.01:10;]

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

    return reshape(GMM_value, (length([-10:0.01:10;]), length([0:0.01:0.6;])))
    #return GMM_value
end


#---------------------------------------------------------------------------------------------
# For testing the functions (Need to load the functions in 01simulate_data file)

parameter = market_parameters_log(T = 1000, σ = 1);
data = simulation_data_log(parameter);

data = load("../conduct_parameter/output/data_loglinear_loglinear_n_1000_sigma_1.rds");

describe(data)


@time simultaneous_test = simulation_nonlinear_2SLS(parameter, data, (:simultaneous, :non_constraint));
describe(simultaneous_test)
@time simultaneous_constraint_test = simulation_nonlinear_2SLS(parameter, data,  (:simultaneous, :constraint));
describe(simultaneous_test)


@time separate_test = simulation_nonlinear_2SLS(parameter, data, (:separate, :non_constraint));
describe(separate_test)
@time separate_constraint_test = simulation_nonlinear_2SLS(parameter, data, (:separate, :constraint));
describe(separate_constraint_test)


# Draw the contour set of the GMM function



test = GMM_function_separate(parameter, data);
test_plot = plot(contour([0:0.01:0.6;], [-10:0.01:10;],test,
    xlabel="θ", ylabel="γ_0",
    title="Value of GMM, N =200, σ = 1"))
    vline!([0.5], linestyle=:dash)
    hline!([1], linestyle=:dash)




    test = rand(1000, 100)

    test_index = setdiff([1:1000;], sample([1:1000;], 90, replace = false))


    test= test[test_index, :]







    function GMM_estimation_separate(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, parameter)
    
        @unpack θ = parameter
    
        QZ_hat = Z_d * inv(Z_d' * Z_d) * Z_d' * (Z_d[:,2] .* Q)
        Q_hat = Z_d * inv(Z_d' * Z_d) * Z_d' *  Q
        X_dd = hcat(ones(T), -Q_hat, -QZ_hat, X_d[:,end])
        α_hat = inv(X_dd' * X_dd) * (X_dd' * P)
    
        L = size(Z, 2)
        L_d = size(Z_d,2)
        L_s = size(Z_s,2)
        K_s = size(X_s, 2)
        K_d = size(X_d, 2)
    
    
        if sum(1 .- θ .*(α_hat[2] .+ α_hat[3] .* X_s[:,end]) .<= 0) == 0
    
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
    
    
            if  sum(1 .- θ_hat .*(α_hat[2] .+ α_hat[3] .* X_s[:,end]) .<= 0) == 0
    
                return α_hat, γ_hat, θ_hat, termination_status_code(termination_status(model))
            else 
                error("The estimation result violates the model assumption ")
            end
        else
    
            γ_hat = repeat([missing], K_s-1)
            θ_hat = missing
    
            return α_hat, γ_hat, θ_hat, missing
        end
    
    
    end






#--------------------------------------------------------------------------------------------------------


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
    @variable(model, θ)

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

    if  sum(1 .- θ_hat .*(α_hat[2] .+ α_hat[3] .* X_s[:,end]) .<= 0) == 0

        return α_hat, γ_hat, θ_hat, termination_status_code(termination_status(model))
    else 
        error("The estimation result violates the model assumption ")

    end
    
end

function estimation_nonlinear_2SLS_simultaneous(parameter, data)
    @unpack T = parameter

    Q  = data.logQ
    w  = data.w
    r  = data.r
    z  = data.z
    iv_w = data.iv_w
    iv_r = data.iv_r
    p  = data.logP
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
        Z_st = vcat(1, z[t], log(w[t]), log(r[t]), y[t])
        Z_t = [Z_dt zeros(length(Z_dt));  zeros(length(Z_st)) Z_st]'

        X_dt = vcat(1, -Q[t], -z[t].*Q[t], y[t])
        X_st = vcat(1, Q[t], log(w[t]), log(r[t]), z[t])
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
    γ_est = Vector{Union{Missing, Float64}}[]
    θ_est  = Union{Missing, Float64}[]
    status = Union{Missing, Int64}[]

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

function GMM_estimation_separate(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, parameter)
    
    @unpack θ = parameter

    QZ_hat = Z_d * inv(Z_d' * Z_d) * Z_d' * (Z_d[:,2] .* Q)
    Q_hat = Z_d * inv(Z_d' * Z_d) * Z_d' *  Q
    X_dd = hcat(ones(T), -Q_hat, -QZ_hat, X_d[:,end])
    α_hat = inv(X_dd' * X_dd) * (X_dd' * P)

    L = size(Z, 2)
    L_d = size(Z_d,2)
    L_s = size(Z_s,2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    sample_violatiton_index = Int64[]

    # Pick up the index of market under which the inside of the log has a negative value
    for t = 1:T
        if 1 .- θ .*(α_hat[2] .+ α_hat[3] .* X_s[t,end]) <= 0
            push!(sample_violatiton_index, t)
        end
    end

    # If all markets do not satisfy the assumption, stop the supply estimation and rerutns missing values
    if length(sample_violatiton_index) == T

        γ_hat = repeat([missing], K_s-1)
        θ_hat = missing

        return α_hat, γ_hat, θ_hat, missing

    else
        
        # Drop the samples that violate the assumption
        if  1 <= length(sample_violatiton_index)

            sample_index = setdiff([1:T;], sample_violatiton_index)

            Z_s = Z_s[sample_index, :]
            X_s = X_s[sample_index, :]
            T = length(sample_index)
            L_s = size(Z_s,2)
            K_s = size(X_s, 2)
        end

        #Check if the weight matrix can be obtained
        if rank(Z_s' * Z_s) == L_s
            
            Ω = inv(Z_s' * Z_s)/T

            model = Model(Ipopt.Optimizer)
            set_optimizer_attribute(model, "tol", 1e-15)
            set_optimizer_attribute(model, "max_iter", 1000)
            set_optimizer_attribute(model, "acceptable_tol", 1e-12)
            set_silent(model)
            @variable(model, γ[k = 1:K_s-1])
            @variable(model, θ)

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
        else

            # If the weight matrix can not be obtained, return missing values
            γ_hat = repeat([missing], K_s-1)
            θ_hat = missing

            return α_hat, γ_hat, θ_hat, missing
        end

        # Check if the supply estimation result satisfies the assumption
        if  sum(1 .- θ_hat .*(α_hat[2] .+ α_hat[3] .* X_s[:,end]) .<= 0) == 0

            return α_hat, γ_hat, θ_hat, termination_status_code(termination_status(model))
        else 
            error("The estimation result violates the model assumption ")
        end
    end

end


function estimation_nonlinear_2SLS_separate(parameter, data)
    @unpack T = parameter

    Q  = data.logQ
    w  = data.w
    r  = data.r
    z  = data.z
    iv_w = data.iv_w
    iv_r = data.iv_r
    p  = data.logP
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
        Z_st = vcat(1, z[t], log(w[t]), log(r[t]), y[t])
        Z_t = [Z_dt zeros(length(Z_dt));  zeros(length(Z_st)) Z_st]'

        X_dt = vcat(1, -Q[t], -z[t].*Q[t], y[t])
        X_st = vcat(1, Q[t], log(w[t]), log(r[t]), z[t])
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


    α_hat, γ_hat, θ_hat, status = GMM_estimation_separate(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, parameter)

    return α_hat, γ_hat, θ_hat, status
end


function simulation_nonlinear_2SLS_separate(parameter, data)

    @unpack T, S, σ = parameter 

    α_est  = Vector{Float64}[]
    γ_est  = Vector{Union{Missing, Float64}}[]
    θ_est  = Union{Missing, Float64}[]
    status = Union{Missing, Int64}[]

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
