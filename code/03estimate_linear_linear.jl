
using LinearAlgebra, Distributions
using Statistics, Random, MultivariateStats
using JuMP, Ipopt
using DelimitedFiles, JLD, CSV, DataFrames
using Plots, Combinatorics, Dates, StatsPlots
using Parameters: @unpack, @with_kw

#-----------------------------------------------------------------------------------------------------

market_parameters = @with_kw (
    α_0 = 10,   # Demand parameter
    α_1 = 1,    
    α_2 = 1,    
    α_3 = 1,    
    γ_0 = 1,    # Marginal cost parameter
    γ_1 = 1,    
    γ_2 = 1,    
    γ_3 = 1,    
    θ = 0.5,    # Conduct parameter
    σ = 1,      # The standard deviation of the error terms
    T = 50,     # The number of market
    S = 1000,   # The number of simulation
)


# The following functions are for implementing 2SLS separately

function estimaton_linear_2SLS_separate(parameter, data)

    # Separately implement 2SLS to the demand equation and the supply equation

    Q  = data.Q
    w  = data.w
    r  = data.r
    z  = data.z
    iv_w = data.iv_w
    iv_r = data.iv_r
    P  = data.P
    y  = data.y

    iv = hcat(iv_w, iv_r)

    @unpack  T = parameter

    # Demand estimation by 2SLS

    Z_demand = hcat(ones(T), z, iv, y)
    QZ_hat = Z_demand * inv(Z_demand' * Z_demand) * Z_demand' * (z .* Q)
    Q_hat = Z_demand * inv(Z_demand' * Z_demand) * Z_demand' *  Q
    X_d = hcat(ones(T), -Q_hat, -QZ_hat, y)
    α_hat = inv(X_d' * X_d) * (X_d' * P)

    
    # Supply side estimation by 2SLS

    Z_supply = hcat(ones(T), w, r, z, y)
    Q_hat = Z_supply * inv(Z_supply' * Z_supply) * Z_supply' *  Q
    QZ_hat = Z_supply * inv(Z_supply' * Z_supply) * Z_supply' *  ((α_hat[2] .+ α_hat[3] .* z) .*Q)
    X_s = hcat(ones(T), Q_hat, w, r,  QZ_hat)
    γ_hat = inv(X_s' * X_s) * (X_s' * P)
    
    θ_hat = γ_hat[end]

    return α_hat, γ_hat[1:end-1], θ_hat
end


function simulation_2SLS_separate(parameter, data)

    α_est = Vector{Float64}[]
    γ_est = Vector{Float64}[]
    θ_est = Float64[]

    @unpack T, S , σ = parameter 

    for s = 1:S
        data_s = data[(s-1)*T+1:s*T,:]

        α_est_s, γ_est_s, θ_est_s = estimaton_linear_2SLS_separate(parameter, data_s)

        push!(α_est, α_est_s)
        push!(γ_est, γ_est_s)
        push!(θ_est, θ_est_s)
    end

    α_est = reduce(vcat, α_est')
    γ_est = reduce(vcat, γ_est')
    θ_est = reduce(vcat, θ_est')

    estimation_result = DataFrame(
    σ = σ,
    T = T,
    α_0 = α_est[:,1],
    α_1 = α_est[:,2],
    α_2 = α_est[:,3],
    α_3 = α_est[:,4],
    γ_0 = γ_est[:,1],
    γ_1 = γ_est[:,2],
    γ_2 = γ_est[:,3],
    γ_3 = γ_est[:,4],
    θ = θ_est)


    return estimation_result
end

#-----------------------------------------------------------------------------------------------------------------
# The following functions are for implementing the system of 2SLS

function estimaton_linear_2SLS_simultaneous(parameter, data)


    # Implement the system of 2SLS

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

    P = Vector{Float64}[]
    Z = Matrix{Float64}[]
    X = Matrix{Float64}[]

    for t = 1:T
        Z_td = vcat(1, z[t], y[t], iv[t,:])
        Z_ts = vcat(1, z[t], w[t], r[t], y[t])

        Z_t = [Z_td zeros(length(Z_td));  zeros(length(Z_ts)) Z_ts]'

        X_td = vcat(1, -Q[t],-z[t].*Q[t], y[t])
        X_ts = vcat(1, Q[t], w[t], r[t], z[t].*Q[t])

        X_t = [X_td zeros(length(X_td)); zeros(length(X_ts)) X_ts]'

        push!(P, vcat(p[t], p[t]))
        push!(X, X_t)
        push!(Z, Z_t)
    end

    Z = reduce(vcat,(Z))
    X = reduce(vcat,(X))
    P = reshape(reduce(vcat,transpose.(P))', (T * 2))

    # The system of 2SLS estimator
    β_hat = inv(X' * Z * inv(Z'Z) * Z' * X) * (X' * Z * inv(Z'Z) * Z' * P)

    # Save the estimation results
    α_0_hat, α_1_hat, α_2_hat, α_3_hat = β_hat[1], β_hat[2], β_hat[3], β_hat[4]
    γ_0_hat, γ_2_hat, γ_3_hat = β_hat[5], β_hat[7], β_hat[8]

    θ_hat = β_hat[9]/α_2_hat
    γ_1_hat = β_hat[6] - θ_hat * α_1_hat

    return vcat(α_0_hat, α_1_hat, α_2_hat, α_3_hat), vcat(γ_0_hat, γ_1_hat, γ_2_hat, γ_3_hat), θ_hat
end


function simulation_2SLS_simultaneous(parameter, data)

    @unpack T, S , σ = parameter 

    α_est = Vector{Float64}[]
    γ_est = Vector{Float64}[]
    θ_est = Float64[]

    for s = 1:S
        data_s = data[(s-1)*T+1:s*T,:]            
        
        α_est_s, γ_est_s, θ_est_s = estimaton_linear_2SLS_simultaneous(parameter, data_s)

        push!(α_est, α_est_s)
        push!(γ_est, γ_est_s)
        push!(θ_est, θ_est_s)
    end

    α_est = reduce(vcat, α_est')
    γ_est = reduce(vcat, γ_est')
    θ_est = reduce(vcat, θ_est')

    estimation_result = DataFrame(
    σ = σ,
    T = T,
    α_0 = α_est[:,1],
    α_1 = α_est[:,2],
    α_2 = α_est[:,3],
    α_3 = α_est[:,end],
    γ_0 = γ_est[:,1],
    γ_1 = γ_est[:,2],
    γ_2 = γ_est[:,3],
    γ_3 = γ_est[:,4],
    θ = θ_est)


    return estimation_result
end


#--------------------------------------------------------------------------------
# Estimate the parameters for each number of markets and the value of the standard deviation of the error terms

for t = [50, 100, 200, 1000], sigma =  [0.001, 0.5, 1, 2]


    # Load the simulation data from the rds files
    filename_begin = "../conduct_parameter/output/data_linear_linear_n_"
    filename_end   = ".rds"

    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end

    filename = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end

    data = load(filename)
    data = DataFrames.sort(data, [:group_id_k])

    
    #= Undo the following lines to load the simulation data generated by julia code

        filename_begin = "../conduct_parameter/output/data_linear_linear_n_"
        filename_end   = "_with_demand_shifter_y.csv"

        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end

        data = DataFrame(CSV.File(file_name))
    =#

    # Set paramter values
    parameter = market_parameters(T = t, σ = sigma)
    
    # Estimate the parameter based on 2SLS
    estimation_result = simulation_2SLS_separate(parameter, data)

    # Uncomment below for checking the estimation result by the system of 2SLS
    # estimation_result = simulation_2SLS_simultaneous(parameter, data)

    # Save the estimatin results as csv files. The file is found at "output" folder
    filename_begin = "../conduct_parameter/output/parameter_hat_table_linear_linear_n_"
    filename_end   = ".csv"

    file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end

    CSV.write(file_name, estimation_result)
end





#------------------------------------------------------------------------------------------------

# Estimation that uses the data without demand shifter

#------------------------------------------------------------------------------------------------


# The following functions are for implementing 2SLS separately

function estimaton_linear_2SLS_separate_wo_demand_shifter(parameter, data)

    # Separately implement 2SLS to the demand equation and the supply equation

    Q  = data.Q
    w  = data.w
    r  = data.r
    z  = data.z
    iv_w = data.iv_w
    iv_r = data.iv_r
    P  = data.P
    
    iv = hcat(iv_w, iv_r)

    @unpack  T = parameter

    # Demand estimation by 2SLS

    Z_demand = hcat(ones(T), z, iv)
    QZ_hat = Z_demand * inv(Z_demand' * Z_demand) * Z_demand' * (z .* Q)
    Q_hat = Z_demand * inv(Z_demand' * Z_demand) * Z_demand' *  Q
    X_d = hcat(ones(T), -Q_hat, -QZ_hat)
    α_hat = inv(X_d' * X_d) * (X_d' * P)

    
    # Supply side estimation by 2SLS

    Z_supply = hcat(ones(T), w, r, z)
    Q_hat  = Z_supply * inv(Z_supply' * Z_supply) * Z_supply' *  Q
    QZ_hat = Z_supply * inv(Z_supply' * Z_supply) * Z_supply' *  ((α_hat[2] .+ α_hat[3] .* z) .*Q)
    X_s = hcat(ones(T), Q_hat, w, r,  QZ_hat)

    display(X_s)
    γ_hat = inv(X_s' * X_s) * (X_s' * P)
    
    θ_hat = γ_hat[end]

    return α_hat, γ_hat[1:end-1], θ_hat
end


function simulation_2SLS_separate_wo_demand_shifter(parameter, data)

    α_est = Vector{Float64}[]
    γ_est = Vector{Float64}[]
    θ_est = Float64[]

    @unpack T, S , σ = parameter 

    for s = 1:S
        data_s = data[(s-1)*T+1:s*T,:]

        α_est_s, γ_est_s, θ_est_s = estimaton_linear_2SLS_separate_wo_demand_shifter(parameter, data_s)

        push!(α_est, α_est_s)
        push!(γ_est, γ_est_s)
        push!(θ_est, θ_est_s)
    end

    α_est = reduce(vcat, α_est')
    γ_est = reduce(vcat, γ_est')
    θ_est = reduce(vcat, θ_est')

    estimation_result = DataFrame(
    σ = σ,
    T = T,
    α_0 = α_est[:,1],
    α_1 = α_est[:,2],
    α_2 = α_est[:,3],
    γ_0 = γ_est[:,1],
    γ_1 = γ_est[:,2],
    γ_2 = γ_est[:,3],
    γ_3 = γ_est[:,4],
    θ = θ_est)


    return estimation_result
end

#-----------------------------------------------------------------------------------------------------------------
# The following functions are for implementing the system of 2SLS

function estimaton_linear_2SLS_simultaneous_wo_demand_shifter(parameter, data)


    # Implement the system of 2SLS

    @unpack T = parameter

    Q  = data.Q
    w  = data.w
    r  = data.r
    z  = data.z
    iv_w = data.iv_w
    iv_r = data.iv_r
    p  = data.P

    iv = hcat(iv_w, iv_r)

    P = Vector{Float64}[]
    Z = Matrix{Float64}[]
    X = Matrix{Float64}[]

    for t = 1:T
        Z_td = vcat(1, z[t],  iv[t,:])
        Z_ts = vcat(1, z[t], w[t], r[t])

        Z_t = [Z_td zeros(length(Z_td));  zeros(length(Z_ts)) Z_ts]'

        X_td = vcat(1, -Q[t],-z[t].*Q[t])
        X_ts = vcat(1, Q[t], w[t], r[t], z[t].*Q[t])

        X_t = [X_td zeros(length(X_td)); zeros(length(X_ts)) X_ts]'

        push!(P, vcat(p[t], p[t]))
        push!(X, X_t)
        push!(Z, Z_t)
    end

    Z = reduce(vcat,(Z))
    X = reduce(vcat,(X))
    P = reshape(reduce(vcat,transpose.(P))', (T * 2))

    # The system of 2SLS estimator
    β_hat = inv(X' * Z * inv(Z'Z) * Z' * X) * (X' * Z * inv(Z'Z) * Z' * P)

    # Save the estimation results
    α_0_hat, α_1_hat, α_2_hat = β_hat[1], β_hat[2], β_hat[3]
    γ_0_hat, γ_2_hat, γ_3_hat = β_hat[4], β_hat[6], β_hat[7]

    θ_hat = β_hat[8]/α_2_hat
    γ_1_hat = β_hat[5] - θ_hat * α_1_hat

    return vcat(α_0_hat, α_1_hat, α_2_hat), vcat(γ_0_hat, γ_1_hat, γ_2_hat, γ_3_hat), θ_hat
end


function simulation_2SLS_simultaneous_wo_demand_shifter(parameter, data)

    @unpack T, S , σ = parameter 

    α_est = Vector{Float64}[]
    γ_est = Vector{Float64}[]
    θ_est = Float64[]

    for s = 1:S
        data_s = data[(s-1)*T+1:s*T,:]            
        
        α_est_s, γ_est_s, θ_est_s = estimaton_linear_2SLS_simultaneous_wo_demand_shifter(parameter, data_s)

        push!(α_est, α_est_s)
        push!(γ_est, γ_est_s)
        push!(θ_est, θ_est_s)
    end

    α_est = reduce(vcat, α_est')
    γ_est = reduce(vcat, γ_est')
    θ_est = reduce(vcat, θ_est')

    estimation_result = DataFrame(
    σ = σ,
    T = T,
    α_0 = α_est[:,1],
    α_1 = α_est[:,2],
    α_2 = α_est[:,3],
    γ_0 = γ_est[:,1],
    γ_1 = γ_est[:,2],
    γ_2 = γ_est[:,3],
    γ_3 = γ_est[:,4],
    θ = θ_est)


    return estimation_result
end



#--------------------------------------------------------------------------------
# Estimate the parameters for each number of markets and the value of the standard deviation of the error terms


for t = [50, 100, 200, 1000], sigma =  [0.001, 0.5, 1, 2]


    # Load the simulation data from the rds files
    filename_begin = "../conduct_parameter/output/data_linear_linear_n_"
    filename_end   = "_without_demand_shifter_y.rds"

    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end

    filename = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end

    data = load(filename)
    data = DataFrames.sort(data, [:group_id_k])

    #Undo the following lines to load the simulation data generated by julia code

        #filename_begin = "../conduct_parameter/output/data_linear_linear_n_"
        #filename_end   = "_with_demand_shifter_y.csv"

        #file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end

        #data = DataFrame(CSV.File(file_name))


    # Set paramter values
    parameter = market_parameters(T = t, σ = sigma)
    
    # Estimate the parameter based on 2SLS
    estimation_result = simulation_2SLS_separate_wo_demand_shifter(parameter, data)

    # Uncomment below for checking the estimation result by the system of 2SLS
    # estimation_result = simulation_2SLS_simultaneous_wo_demand_shifter(parameter, data)

    # Save the estimatin results as csv files. The file is found at "output" folder
    filename_begin = "../conduct_parameter/output/parameter_hat_table_linear_linear_n_"
    filename_end   = "_with_demand_shifter_y.csv"

    file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end

    CSV.write(file_name, estimation_result)
end