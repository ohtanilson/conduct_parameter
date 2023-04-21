include("00setting_julia.jl")
include("00functions.jl")
parameter = market_parameters_log()
starting_value = :default
tol_level = :loose

# Check the performance of the MPEC method by using the linear specification
# The estimation outcome will be stored in "output"

for estimation_method = [(:mpec_linear, :non_constraint, :theta_constraint)]
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
        # Set parameter values
        parameter = market_parameters_log(T = t, σ = sigma)
        
        # Estimation based on 2SLS
        @time estimation_result = iterate_esimation_nonlinear_2SLS(parameter, data, estimation_method, start_value, tol_level)

        # Save the estimation result as csv file. The file is saved at "output" folder
        filename_estimation = "_"*String(estimation_method[1])*"_"*String(estimation_method[2])*"_"*String(estimation_method[3])

        filename_begin = "../conduct_parameter/output/parameter_hat_table_linear_linear_n_"
        filename_end   = ".csv"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end
        print("Simulate : $file_name \n")

        CSV.write(file_name, estimation_result, transform=(col, val) -> something(val, missing))
    end
    println("\n")
    println("----------------------------------------------------------------------------------\n")
end