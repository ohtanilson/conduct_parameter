"""

Monte Carlo Simulation with linear model and compute optimal instrument variable

"""


##

#----------------------------------------------------------------------------------
# We assume that the reader already downloaded the above packages
using LinearAlgebra, Distributed
using Random
using Distributions
using CSV
using DataFrames
using Plots
using VegaLite
using Parameters: @with_kw, @unpack
using JuMP, Ipopt, Optim
using Distributed

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
    estimation_method::Symbol
    starting_value::Symbol
    tol_level::Symbol
    simulation_index::Int64
end

##

estimation_methods = [
    #(:linear,:theta_constraint, :slope_constraint, :equilibrium_constraint), 
    :linear_optimal_separate,
    #:linear_optimal_simultaneous
    ]
starting_value = :true_value
tol_level = :loose

# Estimate the parameters for each number of markets and the value of the standard deviation of the error terms
for estimation_method = estimation_methods
    for t = [100], sigma =  [1], theta = [1], alpha = [1]
        # Load the simulation data from the rds files
        filename_begin = "../conduct_parameter/output/testing_project/data_linear_linear_n_"
        filename_end   = ".rds"

        if sigma == 1 || sigma == 2
            sigma = Int64(sigma)
        end
        filename = filename_begin*string(t)*"_theta_"*string(theta)*"_alpha2_"*string(alpha)*"_sigma_"*string(sigma)*filename_end

        data = load(filename)
        data = DataFrames.sort(data, [:group_id_k])
        # Set parameter values
        parameter = market_parameters(T = t, σ = sigma)
        
        # Estimation based on 2SLS
        @time estimation_result = simulation_GMM_optimal_instrument(parameter, data, estimation_method, starting_value, tol_level)

        # Save the estimation result as csv file. The file is saved at "output" folder

        filename_begin = "../conduct_parameter/output/testing_project/parameter_hat_table_linear_linear_n_"
        filename_end   = ".csv"
        file_name = filename_begin*string(t)*"_theta_"*string(theta)*"_alpha_"*string(alpha)*"_sigma_"*string(sigma)*"_"*String(estimation_method)*filename_end
        print("Simulate : $file_name \n")

        #CSV.write(file_name, estimation_result, transform=(col, val) -> something(val, missing))
    end
    println("\n")
    println("----------------------------------------------------------------------------------\n")
end