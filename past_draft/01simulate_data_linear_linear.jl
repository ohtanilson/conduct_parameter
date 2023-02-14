using LinearAlgebra, Distributions
using Statistics, Random, MultivariateStats
using CSV, DataFrames
using Parameters: @unpack, @with_kw


#--------------------------------------------------------------------------------------------------------------
# Set parameters
market_parameters = @with_kw (
    α_0 = 10,   # Demand parameter
    α_1 = 1,    # Demand parameter
    α_2 = 1,    # Demand parameter
    α_3 = 1,    # Demand parameter
    γ_0 = 1,    # Marginal cost parameter
    γ_1 = 1,    # Marginal cost parameter
    γ_2 = 1,    # Marginal cost parameter
    γ_3 = 1,    # Marginal cost parameter
    θ = 0.5,    # Conduct parameter
    σ = 0.000001,      # The standard deviation of the error terms
    T = 50,     # The number of market
    S = 1000,   # The number of simulation
)

mutable struct market_data
    Q::Vector{Float64}
    W::Vector{Float64}
    R::Vector{Float64}
    Z::Vector{Float64}
    IV_W::Vector{Float64}
    IV_R::Vector{Float64}
    P::Vector{Float64}
    Y::Vector{Float64}
end

#------------------------------------------------------------------------------------------------------------

function simulation_data(parameter)
    """
    Generate 1000 simulation data with demand shifter Y_t
    """

    @unpack α_0, α_1, α_2, α_3,γ_0 , γ_1 ,γ_2 ,γ_3, θ,σ ,T, S = parameter
    
    Q = Float64[];
    P = Float64[];
    W = Float64[];
    R = Float64[];
    Z = Float64[];
    Y = Float64[];
    IV_W = Float64[];
    IV_R = Float64[];

    group_id = Int64[];

    for s = 1:S
        Random.seed!(s * 12345)

        for t = 1:T
            z_t = rand(Normal(10,1))
            w_t = rand(Normal(3,1))
            r_t = rand(Normal(0,1))
            iv_w_t = w_t + randn()
            iv_r_t = r_t + randn()
            Y_t = randn()

            ε_d = rand(Normal(0,σ))
            ε_c = rand(Normal(0,σ))

            Q_t = (α_0 + α_3 * Y_t - γ_0 - γ_2 * w_t - γ_3 * r_t + ε_d - ε_c)/((1+θ) * (α_1 + α_2 *z_t) + γ_1)

            p_t = α_0 - α_1 * Q_t  - α_2 * z_t * Q_t + α_3 * Y_t + ε_d 

            push!(group_id, s)
            push!(Q, Q_t)
            push!(P, p_t)
            push!(W, w_t)
            push!(R, r_t)
            push!(Z, z_t)
            push!(Y, Y_t)
            push!(IV_W, iv_w_t)
            push!(IV_R, iv_r_t)
        end
    end

    data = DataFrame(group_id_k = group_id, Q = Q, P = P, w = W, r = R, y = Y, z = Z, iv_w = IV_W, iv_r = IV_R)

    return data
end


function simulation_data_without_demand_shifter(parameter)
    """
    Generate 1000 simulation data without demand shifter Y_t
    """


    @unpack α_0, α_1, α_2, γ_0 , γ_1 ,γ_2 ,γ_3, θ,σ ,T, S = parameter
    
    Q = Float64[];
    P = Float64[];
    W = Float64[];
    R = Float64[];
    Z = Float64[];
    IV_W = Float64[];
    IV_R = Float64[];

    group_id = Int64[];

    for s = 1:S
        Random.seed!(s * 12345)

        for t = 1:T
            z_t = rand(Normal(10,1))
            w_t = rand(Normal(3,1))
            r_t = rand(Normal(0,1))
            iv_w_t = w_t + randn()
            iv_r_t = r_t + randn()

            ε_d = rand(Normal(0,σ))
            ε_c = rand(Normal(0,σ))

            Q_t = (α_0  - γ_0 - γ_2 * w_t - γ_3 * r_t + ε_d - ε_c)/((1+θ) * (α_1 + α_2 *z_t) + γ_1)

            p_t = α_0 - α_1 * Q_t  - α_2 * z_t * Q_t + ε_d 

            push!(group_id, s)
            push!(Q, Q_t)
            push!(P, p_t)
            push!(W, w_t)
            push!(R, r_t)
            push!(Z, z_t)
            push!(IV_W, iv_w_t)
            push!(IV_R, iv_r_t)
        end
    end

    data = DataFrame(group_id_k = group_id, Q = Q, P = P, w = W, r = R, z = Z, iv_w = IV_W, iv_r = IV_R)

    return data
end


#-----------------------------------------------------------------------------------------------------------

# Simulation with demand shifter
# Generate the simulation data and save the data as CSV file        
for t in [50, 100, 200, 1000],  sigma in [0.001, 0.5, 1, 2] 
        
    parameter = market_parameters(T = t, σ = sigma);
    data = simulation_data(parameter);


    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end


    @unpack α_1, α_2 = parameter;
    if sum(0 .<= -(α_1 .+ α_2 * data.z)) == 0

        file_name = "../conduct_parameter/output/data_linear_linear_n_"*string(t)*"_sigma_"*string(sigma)*".csv"

        CSV.write(file_name, data)
    else
        error("The demand function is not downward sloping")
    end
end


# Simulation without demand shifter y
# Generate the simulation data and save the data as CSV file
for t in [50, 100, 200, 1000],  sigma in [0.001, 0.5, 1, 2] 
        
    parameter = market_parameters(T = t, σ = sigma);
    data = simulation_data_without_demand_shifter(parameter);


    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end

    @unpack α_1, α_2 = parameter;
    if sum(0 .<= -(α_1 .+ α_2 * data.z)) == 0

        file_name = "../conduct_parameter/past_draft/data_linear_linear_n_"*string(t)*"_sigma_"*string(sigma)*"_without_demand_shifter_y"*".csv"

        CSV.write(file_name, data)
    else
        error("The demand function is not downward sloping")
    end
end