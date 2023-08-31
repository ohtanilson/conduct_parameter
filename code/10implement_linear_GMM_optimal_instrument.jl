"""

Monte Carlo Simulation with linear model and compute optimal instrument variable

"""


##

#----------------------------------------------------------------------------------
# We assume that the reader already downloaded the above packages
using LinearAlgebra
using Random
using Distributions
using CSV
using DataFrames
using Plots
using VegaLite
using Parameters: @with_kw, @unpack
using JuMP, Ipopt, Optim


##


# Set parameters for linear model
market_parameters = @with_kw (
    α_0 = 10, # Demand parameter
    α_1 = 1,
    α_2 = 1,
    α_3 = 1,
    γ_0 = 1,  # Marginal cost parameter
    γ_1 = 1,
    γ_2 = 1,
    γ_3 = 1,
    θ_0 = 0.5,  # Conduct paramter
    σ = 1,    # Standard deviation of the error term
    T = 50,   # Number of markets
    S = 1000,   # Number of simulation
    start_θ = 0.0,
    start_γ = [0.0, 0.0, 0.0, 0.0]
)

mutable struct SIMULATION_SETTING
    α_0::Float64
    α_1::Float64
    α_2::Float64 
    α_3::Float64
    γ_0::Float64
    γ_1::Float64
    γ_2::Float64
    γ_3::Float64
    θ_0::Float64
    start_θ::Float64
    start_γ::Vector{Float64}
    T::Int64
    σ::Float64
    data::DataFrame
    estimation_method::Tuple{Symbol, Symbol, Symbol,Symbol}
    starting_value::Symbol
    tol_level::Symbol
    simulation_index::Int64
end


##




function GMM_estimation_linear(T, P, Z, X, X_s, X_d, Ω, α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, starting_value, tol_level)
    
    """ 
    Estimate the demand and supply parameter given a market simultaneously
    """

    L = size(Z, 2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)


    if tol_level == :tight
        tol = 1e-15
        acceptable_tol = 1e-12
    elseif tol_level == :loose
        tol = 1e-6
        acceptable_tol = 1e-5
    end

    if starting_value == :true_value
        start_β = [α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3]
        start_θ = θ_0

    elseif starting_value == :random
        start_β = [α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3] .+ rand(Uniform(-10, 10), 8)
        start_θ = θ_0 + rand(Uniform(-10, 1))
    end

    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "tol", tol)
    set_optimizer_attribute(model, "max_iter", 1000)
    set_optimizer_attribute(model, "acceptable_tol", acceptable_tol)
    set_silent(model)
    @variable(model, β[k = 1:K_d+K_s-1])
    @variable(model, 0 <= θ <= 1)


    r = Any[];
    for t =1:T
        push!(r, @NLexpression(model, P[t] - sum(β[k] * X[2*t-1,k] for k = 1:K_d) ))
        push!(r, @NLexpression(model, P[t] - θ * (β[2] + β[3] * X[2*t, end]) * X[2*t,K_d + 2] - sum(β[k] * X[2*t,k] for k = K_d+1:K_d+K_s-1)))
    end

    g = Any[];
    for l = 1:L
        push!(g, @NLexpression(model, sum(Z[t,l] * r[t] for t = 1:2*T)))
    end

    @NLobjective(model, Min, sum( g[l] *Ω[l,k] * g[k] for l = 1:L, k = 1:L))
    optimize!(model)
    
    α_hat = value.(β)[1:K_d]
    γ_hat = value.(β)[K_d+1:end]
    θ_hat = value.(θ)

    return α_hat, γ_hat, θ_hat, termination_status_code(termination_status(model))
end






@everywhere function estimate_linear_optimal_GMM(simulation_setting::SIMULATION_SETTING)
    """
    Given data, reshape the data and pass it to the function that implement the GMM estimation

    """

    α_0 = simulation_setting.α_0
    α_1 = simulation_setting.α_1
    α_2 = simulation_setting.α_2
    α_3 = simulation_setting.α_3
    γ_0 = simulation_setting.γ_0
    γ_1 = simulation_setting.γ_1
    γ_2 = simulation_setting.γ_2
    γ_3 = simulation_setting.γ_3
    θ_0 = simulation_setting.θ_0
    start_θ = simulation_setting.start_θ
    start_γ = simulation_setting.start_γ
    T = simulation_setting.T
    data = simulation_setting.data
    estimation_method = simulation_setting.estimation_method
    starting_value = simulation_setting.starting_value
    tol_level = simulation_setting.tol_level
    simulation_index  = simulation_setting.simulation_index

    data_s = data[(simulation_index-1)*T+1:simulation_index*T,:]


    Q  = data_s.Q
    w  = data_s.w
    r  = data_s.r
    p  = data_s.P
    y  = data_s.y

    z  = data_s.z
    iv_w = data_s.iv_w
    iv_r = data_s.iv_r
    
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


    # The wight function for the GMM estimation
    Ω_initial = inv(Z' * Z)/T

    α_hat, γ_hat, θ_hat, status = GMM_estimation_linear(T, P, Z, X, X_s, X_d, Ω_initial, α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, starting_value, tol_level)


    Z_optimal = compute_optimal_instruments(T, α_hat, γ_hat, θ_hat, Z_d, Z_s, X_s, X_d, P )

    @show size(Z), size(Z_optimal)

    Ω_optimal = inv(Z_optimal' * Z_optimal)/T


    α_optimal, γ_optimal, θ_optimal, status = GMM_estimation_linear(T, P, Z_optimal, X, X_s, X_d, Ω_optimal, α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, starting_value, tol_level)

    @show α_optimal, γ_optimal, θ_optimal, status

    #return α_hat, γ_hat, θ_hat, status
    return  α_optimal, γ_optimal, θ_optimal, statusstatus
end



function compute_optimal_instruments(T, α_hat, γ_hat, θ_hat, Z_d, Z_s, X_s, X_d, P )

    K_d = size(X_d, 2)
    K_s = size(X_s, 2) - 1


    ε_d = P .- sum(α_hat[k] * X_d[:, k] for k = 1:K_d)                                                    # T × 1 vector
    ε_s = P .- sum(γ_hat[k] * X_s[:, k] for k = 1:K_s) .- θ_hat * (α_hat[2] .+ α_hat[3] .* Z_s[:,2])

    ε = hcat(ε_d, ε_s)  # T × 2 mateix 

    Ω_inverse = inv(sum(ε[t,:]' * ε[t,:] for t = 1:T)/T)


    Q_bar = (α_hat[1] .+ α_hat[4] .* X_d[:,4] .- γ_hat[1] .- γ_hat[3] .* X_s[:,3]  .- γ_hat[4] .* X_s[:,4] .+  ε_d .- ε_s)./((1 + θ_hat) * (α_hat[2] .+ α_hat[3] .* Z_s[:,2]) .+ γ_hat[2])
    Q_bar_Z = Z_s[:,2] .* Q_bar
    Q_bar_α_Z = (α_hat[2] .+ α_hat[3] .* Z_s[:,2]) .* Q_bar


    Q_bar = mean(Q_bar)
    Q_bar_Z = mean(Q_bar_Z)
    Q_bar_α_Z = mean(Q_bar_α_Z)

    D_z = [-1 Q_bar Q_bar_Z -mean(X_d[:,4]) 0 0 0 0 0; 0 (- θ_hat * Q_bar) (θ_hat*Q_bar_Z) 0 (-1)  (-Q_bar) (-mean(X_s[:,3])) (- mean(X_s[:,4])) -Q_bar_α_Z]

    Z_optimal = D_z * Ω_inverse

    return Z_optimal
end


@everywhere function iterate_esimation_linear_GMM(parameter, data, estimation_method::Tuple{Symbol, Symbol, Symbol, Symbol}, starting_value, tol_level)

    """
    Given the simulation data, run the estimation in each simulation index s = 1,..., 1000, and store the simulation results as a DataFrame file.
    """

    @unpack α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, start_θ, start_γ, T, S, σ = parameter

    α_est  = Vector{Float64}[]
    γ_est  = Vector{Union{Missing, Float64}}[]
    θ_est  = Union{Missing, Float64}[]
    status = Union{Missing, Int64}[]

    simulation_setting = [SIMULATION_SETTING(α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, start_θ, start_γ, T, σ, data, estimation_method, starting_value, tol_level,simulation_index) for simulation_index = 1:S]

    @elapsed @time simulation_mpec_result = pmap(estimate_linear_optimal_GMM, simulation_setting);

    for s = 1:S
        push!(α_est, simulation_mpec_result[s][1])
        push!(γ_est, simulation_mpec_result[s][2])
        push!(θ_est, simulation_mpec_result[s][3])
        push!(status,simulation_mpec_result[s][4])
    end 

    α_est = reduce(vcat, α_est')
    γ_est = reduce(vcat, γ_est')
    θ_est = reduce(vcat, θ_est')
    status = reduce(vcat, status)

    status_indicator = []

    for i = eachindex(status)
        if status[i] !== missing
            if status[i] <= 7
                push!(status_indicator, 1)
            else
                push!(status_indicator, 0)
            end
        else
            push!(status_indicator, missing)
        end
    end

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
    status = status,
    status_indicator = status_indicator)

    return estimation_result
end


##


estimation_methods = [
    #(:linear,:theta_constraint, :slope_constraint, :equilibrium_constraint), 
    (:linear_optimal,:no_constraint, :non_constraint, :non_constraint)
    ]
starting_value = :true_value
tol_level = :loose





# Estimate the parameters for each number of markets and the value of the standard deviation of the error terms
for estimation_method = estimation_methods
    for t = [100], sigma =  [1]
        # Load the simulation data from the rds files
        filename_begin = "../conduct_parameter/output/data_linear_linear_n_"
        filename_end   = ".rds"

        if sigma == 1 || sigma == 2
            sigma = Int64(sigma)
        end
        filename = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end

        data = load(filename)
        data = DataFrames.sort(data, [:group_id_k])
        # Set parameter values
        parameter = market_parameters_log(T = t, σ = sigma)
        
        # Estimation based on 2SLS
        @time estimation_result = iterate_esimation_linear_GMM(parameter, data, estimation_method, starting_value, tol_level)

        # Save the estimation result as csv file. The file is saved at "output" folder
        filename_estimation = "_"*String(estimation_method[1])*"_"*String(estimation_method[2])*"_"*String(estimation_method[3])

        filename_begin = "../conduct_parameter/output/parameter_hat_table_linear_linear_n_"
        filename_end   = ".csv"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end
        print("Simulate : $file_name \n")

        #CSV.write(file_name, estimation_result, transform=(col, val) -> something(val, missing))
    end
    println("\n")
    println("----------------------------------------------------------------------------------\n")
end