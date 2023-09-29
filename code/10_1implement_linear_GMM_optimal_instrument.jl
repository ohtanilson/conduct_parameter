##
#using Distributed
#number_workers = addprocs(Sys.CPU_THREADS-1)
##
using Distributed

##

Distributed.@everywhere include("../code/10_0implement_linear_GMM_optimal_instrument_setting.jl")
Distributed.@everywhere include("../code/10_0implement_linear_GMM_optimal_instrument_function.jl")

estimation_methods = [
    #(:linear,:theta_constraint, :slope_constraint, :equilibrium_constraint), 
    :linear_optimal_separate,
    #:linear_optimal_simultaneous
    ]
starting_value = :true_value
tol_level = :loose

## Estimate the parameters for each number of markets and the value of the standard deviation of the error terms
@time for estimation_method = estimation_methods
    for t = [50, 100, 200, 1000, 2000, 5000, 10000], sigma =  [1], theta = Union{Float64, Int64}[0.05, 0.1, 0.2, 0.33, 0.5, 1], alpha2 = Union{Float64, Int64}[0.1, 0.5, 1, 5, 20]
    #for t = [50], sigma =  [1], theta = Union{Float64, Int64}[0.05], alpha2 = Union{Float64, Int64}[0.1]
        # Load the simulation data from the rds files
        filename_begin = "../conduct_parameter/output/testing_project/data_linear_linear_n_"
        filename_end   = ".rds"

        if sigma == 1 || sigma == 2
            sigma = Int64(sigma)
        end
        filename = filename_begin*string(t)*"_theta_"*string(theta)*"_alpha2_"*string(alpha2)*"_sigma_"*string(sigma)*filename_end

        data = load(filename)
        data = DataFrames.sort(data, [:group_id_k])
        # Set parameter values
        parameter = market_parameters(T = t, σ = sigma, θ_0 = theta, α_2 = alpha2)
        
        # Estimation based on 2SLS
        estimation_result = simulation_GMM_optimal_instrument(parameter, data, estimation_method, starting_value, tol_level)

        # Save the estimation result as csv file. The file is saved at "output" folder

        filename_begin = "../conduct_parameter/output/testing_project/estimation_result/parameter_hat_table_linear_linear_n_"
        filename_end   = ".csv"
        
        file_name = filename_begin*string(t)*"_theta_"*string(theta)*"_alpha2_"*string(alpha2)*"_sigma_"*string(sigma)*"_"*String(estimation_method)*filename_end
        print("Simulate : $file_name \n")

        CSV.write(file_name, estimation_result, transform=(col, val) -> something(val, missing))
    end
    println("\n")
    println("----------------------------------------------------------------------------------\n")
end