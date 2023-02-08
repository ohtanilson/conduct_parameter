include("00setting_julia.jl")
include("00functions.jl")
# change starting value
parameter = market_parameters_log()
estimation_methods = [(:separate,:non_constraint, :non_constraint)];
starting_value_list = [:true, :random];
tol_list = [:tight, :loose];
#-----------------------------------------------------------------------------------------
# Draw histograms consisting of the estimation reuslt of θ
# Histograms for the simulation results that use different combinations of tolerances and starting values
#---------------------------------------------------------------------------------------------------------
for estimation_method = estimation_methods
    for starting_value_used = starting_value_list
        for tol_list_used = tol_list
            for t = [50, 100, 200, 1000], sigma = [0.001, 0.5, 1, 2]
                @unpack θ = parameter

                if sigma == 1 || sigma == 2
                    sigma = Int64(sigma)
                end
                
                # Load the estimation result
                if starting_value_used == :true && tol_list_used == :tight
                    filename_begin = "../conduct_parameter/output/tight/parameter_hat_table_loglinear_loglinear_n_"
                    filename_end   = "_true_start.csv"
                    histo_title = "true start, tight, T = $t, σ = $sigma"
                elseif starting_value_used == :random && tol_list_used == :tight
                    filename_begin = "../conduct_parameter/output/tight/parameter_hat_table_loglinear_loglinear_n_"
                    filename_end   = "_random_start.csv"
                    histo_title = "random start, tight, T = $t, σ = $sigma"
                elseif starting_value_used == :true && tol_list_used == :loose
                    filename_begin = "../conduct_parameter/output/loose/parameter_hat_table_loglinear_loglinear_n_"
                    filename_end   = "_true_start.csv"
                    histo_title = "true start, loose, T = $t, σ = $sigma"
                elseif starting_value_used == :random && tol_list_used == :loose
                    filename_begin = "../conduct_parameter/output/loose/parameter_hat_table_loglinear_loglinear_n_"
                    filename_end   = "_random_start.csv"
                    histo_title = "random start, loose, T = $t, σ = $sigma"
                end

                filename_estimation = "_"*String(estimation_method[1])*"_"*String(estimation_method[2])*"_"*String(estimation_method[3])
                file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end
                estimation_result = DataFrame(CSV.File(file_name))
        
                rate_convergence = count( x -> (x <= 7), dropmissing(estimation_result, :status).status)/length(dropmissing(estimation_result, :status).status)

                println("rate of convergence is ", round(rate_convergence * 100, digits = 1), " when t = $t, σ = $sigma, $starting_value_used, and $tol_list_used")

                # count the number of the estimation result out of [0, 1]
                estimation_result = dropmissing(estimation_result, :θ);
                estimation_result = filter(row -> (row.status_indicator == 1), estimation_result)
        
                number_non_missing = size(estimation_result.θ,1)
                number_out_range = number_non_missing - count(x -> (-10e-9 <= x <= 1 + 1e-8), estimation_result.θ)
                rate_out_range = round(number_out_range/number_non_missing * 100, digits = 3)
                rate_out_range = "$rate_out_range %"
                
                estimation_result = filter(x -> (-2 <= x <= 3), estimation_result.θ)
        
                label_title = "no constraint ($rate_out_range are out of [0, 1])"

                histo_result = Plots.histogram(xlims = [0, 1], title = histo_title, legend = :topright, size = (800, 600))
                Plots.vline!(histo_result, [θ], label = "true value : θ = $θ")
                Plots.histogram!(histo_result, estimation_result, xlims = [-2, 3], xlabel = "θ",ylabel = "counts", label = label_title, fill = true, fillalpha = 0.5, bins = -2:0.1:3)
        
                if starting_value_used == :true && tol_list_used == :tight
                    filename_begin = "../conduct_parameter/figuretable/tight/histogram_loglinear_loglinear_n_"
                    filename_end   = "_true_start.pdf"
                elseif starting_value_used == :random && tol_list_used == :tight
                    filename_begin = "../conduct_parameter/figuretable/tight/histogram_loglinear_loglinear_n_"
                    filename_end   = "_random_start.pdf"
                elseif starting_value_used == :true && tol_list_used == :loose
                    filename_begin = "../conduct_parameter/figuretable/loose/histogram_loglinear_loglinear_n_"
                    filename_end   = "_true_start.pdf"
                elseif starting_value_used == :random && tol_list_used == :loose
                    filename_begin = "../conduct_parameter/figuretable/loose/histogram_loglinear_loglinear_n_"
                    filename_end   = "_random_start.pdf"
                end
                file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*"_"*String(estimation_method[3])filename_end
        
                Plots.display(histo_result)
                Plots.savefig(histo_result, file_name)
            end
        end
    end
end
