include("00setting_julia.jl")
include("00functions.jl")
# change starting value
parameter = market_parameters_log()
estimation_methods = [(:separate,:non_constraint, :non_constraint)];
starting_value_list = [:true, :random];
tol_list = [:tight, :loose];

Random.seed!(1234)
#---------------------------------------------------------------------------------------------------------
# estimate with tight/loose tolerance and true/random starting value
#---------------------------------------------------------------------------------------------------------
for estimation_method = estimation_methods
    for t = [50, 100, 200, 1000], sigma = [0.001, 0.5, 1, 2]
        for starting_value_used = starting_value_list
            for tol_list_used = tol_list
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
                # Save the estimation result as csv file. The file is saved at "output" folder
                filename_estimation = "_"*String(estimation_method[1])*"_"*String(estimation_method[2])*"_"*String(estimation_method[3])
                @time estimation_result = iterate_esimation_nonlinear_2SLS(parameter, data, estimation_method, starting_value_used, tol_list_used)

                if starting_value_used == :true && tol_list_used == :tight
                    filename_begin = "../conduct_parameter/output/tight/parameter_hat_table_loglinear_loglinear_n_"
                    filename_end   = "_true_start.csv"
                elseif starting_value_used == :random && tol_list_used == :tight
                    filename_begin = "../conduct_parameter/output/tight/parameter_hat_table_loglinear_loglinear_n_"
                    filename_end   = "_random_start.csv"
                elseif starting_value_used == :true && tol_list_used == :loose
                    filename_begin = "../conduct_parameter/output/loose/parameter_hat_table_loglinear_loglinear_n_"
                    filename_end   = "_true_start.csv"
                elseif starting_value_used == :random && tol_list_used == :loose
                    filename_begin = "../conduct_parameter/output/loose/parameter_hat_table_loglinear_loglinear_n_"
                    filename_end   = "_random_start.csv"
                end
                file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end
                CSV.write(file_name, estimation_result, transform=(col, val) -> something(val, missing))
            end
        end
    end
    println("----------------------------------------------------------------------------------\n")
end

