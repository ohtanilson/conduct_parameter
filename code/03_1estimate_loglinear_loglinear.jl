using Distributed
Distributed.@everywhere include("../code/00setting_julia.jl")
Distributed.@everywhere include("../code/00functions.jl")
parameter = market_parameters_log()
estimation_methods = 
    [
    (:separate, :non_constraint, :non_constraint),
    (:separate, :non_constraint, :theta_constraint),
    #(:separate, :log_constraint, :theta_constraint),
    (:simultaneous, :non_constraint, :non_constraint),
    (:simultaneous, :non_constraint, :theta_constraint),
    #(:simultaneous, :log_constraint, :theta_constraint)
    #(:optim_nelder_mead, :non_constraint, :theta_constraint)
    ];
starting_value = :default
tol_level = :loose
#--------------------------------------------------------------------------------------------------------------
# Estimate the parameters for each number of markets and the value of the standard deviation of the error terms
#--------------------------------------------------------------------------------------------------------------

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
        # Set parameter values
        parameter = market_parameters_log(T = t, σ = sigma)
        
        # Estimation based on 2SLS
        @time estimation_result = iterate_esimation_nonlinear_2SLS(parameter, data, estimation_method, starting_value, tol_level)

        # Save the estimation result as csv file. The file is saved at "output" folder
        filename_estimation = "_"*String(estimation_method[1])*"_"*String(estimation_method[2])*"_"*String(estimation_method[3])

        filename_begin = "../conduct_parameter/output/parameter_hat_table_loglinear_loglinear_n_"
        filename_end   = ".csv"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end
        print("Simulate : $file_name \n")

        CSV.write(file_name, estimation_result, transform=(col, val) -> something(val, missing))
    end
    println("\n")
    println("----------------------------------------------------------------------------------\n")
end