using LinearAlgebra, Distributions
using Statistics, Random, MultivariateStats
using JuMP, Ipopt
using DelimitedFiles, JLD, CSV, DataFrames
using Plots, Combinatorics, Dates, StatsPlots
using Parameters: @unpack, @with_kw


#---------------------------------------------------------------------------------------------------------

"""

Generate a simulation data under log demand and log marginal cost model

"""



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

mutable struct market_data_log
    Q::Vector{Float64}
    W::Vector{Float64}
    R::Vector{Float64}
    Z::Vector{Float64}
    P::Vector{Float64}
    Y::Vector{Float64}
    MC::Vector{Float64}
    IV::Matrix{Float64}
end


function simulation_data_log(parameter)


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
            r_t = rand(Uniform(0,1))
            z_t = rand(Uniform(0,1))
            w_t = rand(Uniform(1,3))
            iv_w_t = w_t + randn()
            iv_r_t = r_t + randn()

            y_t = randn()
            ε_c = rand(Normal(0,σ))
            ε_d = rand(Normal(0,σ))

            log_Q_t = (α_0 + log(1 - θ * (α_1 + α_2 * z_t))  + α_3 * y_t - ( γ_0 + γ_2 * log(w_t) + γ_3 * log(r_t)) + ε_d- ε_c)/ (γ_1 + α_1 + α_2 * z_t) # Equilibrium total quantity

            log_p_t = α_0 - (α_1 + α_2 * z_t )* log_Q_t  + α_3 * y_t + ε_d # The demand function

            push!(group_id, s)
            push!(Q, log_Q_t)
            push!(P,log_p_t)
            push!(W,log(w_t))
            push!(R,log(r_t))
            push!(Z, z_t)
            push!(Y, y_t)
            push!(IV_W, iv_w_t)
            push!(IV_R, iv_r_t)

        end
    end

    # Save the data as DataFrame
    data = DataFrame(group_id_k = group_id, Q = Q, P = P, w = W, r = R, y = Y, z = Z, iv_w = IV_W, iv_r = IV_R)
    return data
end


#---------------------------------------------------------------------------------------------------------------------


# parameter = market_parameters_log();
# data = simulation_data_log(parameter);

# Plot the scatter plot of P and Q
# plot(scatter(data.P, data.Q, ms=2, ma=0.2), ylabel = "P", xlabel = "Q")

# Generate the simulation data and save the data as CSV file        
for t in [50, 100, 200, 1000],  sigma in [0.001, 0.5, 1, 2] 
        
    parameter = market_parameters_log(T = t, σ = sigma);
    data = simulation_data_log(parameter);


    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end

    @unpack α_1, α_2 = parameter;

    # Save the simulation data if the demand function is downward sloping
    if sum(0 .<= -(α_1 .+ α_2 * data.z)) == 0

        file_name = "../conduct_parameter/output/data_loglinear_loglinear_n_"*string(t)*"_sigma_"*string(sigma)*"_with_demand_shifter_y"*".csv"

        CSV.write(file_name, data)
    else
        error("The demand function is not downward sloping")
    end
end

