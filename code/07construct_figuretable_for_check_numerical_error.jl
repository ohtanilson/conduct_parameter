include(joinpath(dirname(@FILE),"00setting_julia.jl"))
#-----------------------------------------------------------------------------------------

estimation_methods = [(:separate,:non_constraint, :non_constraint)];
# Draw histograms consisting of the estimation reuslt of θ

# Histograms for the simulation results where the starting values are true parameters
for t = [50, 100, 200, 1000], sigma =  [0.001, 0.5, 1, 2]

    @unpack θ = parameter

    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end
    
    histo_result = histogram(xlims = [0, 1], title = "true start, tight, n = $t, σ = $sigma", legend = :topright, size = (800, 600))
    vline!(histo_result, [θ], label = "true value : θ = $θ")

    for estimation_method = estimation_methods
        
        # Load the estimation result
        filename_estimation = "_"*String(estimation_method[1])*"_"*String(estimation_method[2])*"_"*String(estimation_method[3])
        filename_begin = "../conduct_parameter/output/tight/parameter_hat_table_loglinear_loglinear_n_"
        filename_end   = "_true_start.csv"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end
        estimation_result = DataFrame(CSV.File(file_name))

        histo_status = plot(histogram(estimation_result.status), xlims = [1,24], bins = 1:1:24, title = "true start, tight, n = $t, σ = $sigma")

        display(histo_status)
        rate_convergence = count( x -> (x <= 7), dropmissing(estimation_result, :status).status)/length(dropmissing(estimation_result, :status).status)


        # count the number of the estimation result out of [0, 1]
        estimation_result  = dropmissing(estimation_result, :θ);
        estimation_result = filter(row -> (row.status_indicator == 1), estimation_result)

        number_non_missing = size(estimation_result.θ,1)
        number_out_range   = number_non_missing  - count(x -> (-10e-9 <= x <= 1 + 1e-8), estimation_result.θ)
        rate_out_range     = round(number_out_range/number_non_missing * 100, digits = 3)
        rate_out_range     = "$rate_out_range %"
        
        estimation_result = filter(x -> (-2 <= x <= 3), estimation_result.θ)


        if estimation_method[2] == :log_constraint
            if estimation_method[3] == :theta_constraint
                label_title = "both constraints ($rate_out_range are out of [0, 1])"
            else
                label_title = "only log constraint ($rate_out_range are out of [0, 1])"
            end
        else
            if estimation_method[3] == :theta_constraint
                label_title = "with constraint ($rate_out_range are out of [0, 1])"
            else
                label_title = "no constraint ($rate_out_range are out of [0, 1])"
            end
        end

        histogram!(histo_result, estimation_result, xlims = [-2, 3], xlabel = "θ",ylabel = "counts", label = label_title, fill = true, fillalpha = 0.5, bins = -2:0.1:3)

        filename_begin = "../conduct_parameter/figuretable/tight/histogram_loglinear_loglinear_n_"
        filename_end   = "_true_start.pdf"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*"_"*String(estimation_method[3])filename_end

        display(histo_result)

        savefig(histo_result, file_name)

    end
end


# Histograms for the simulation results where the starting values are randomly drawn
for t = [50, 100, 200, 1000], sigma =  [0.001, 0.5, 1, 2]

    @unpack θ = parameter

    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end
    
    histo_result = histogram(xlims = [-10^5, 10^3], title = "random start, tight, n = $t, σ = $sigma", legend = :topright, size = (800, 600))
    vline!(histo_result, [θ], label = "true value : θ = $θ")

    for estimation_method = estimation_methods
        
        # Load the estimation result
        filename_estimation = "_"*String(estimation_method[1])*"_"*String(estimation_method[2])*"_"*String(estimation_method[3])
        filename_begin = "../conduct_parameter/output/tight/parameter_hat_table_loglinear_loglinear_n_"
        filename_end   = "_random_start.csv"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end
        estimation_result = DataFrame(CSV.File(file_name))

        histo_status = plot(histogram(estimation_result.status), xlims = [1,24], bins = 1:1:24, title = "random start, tight, n = $t, σ = $sigma")

        display(histo_status)
        rate_convergence = count( x -> (x <= 7), dropmissing(estimation_result, :status).status)/length(dropmissing(estimation_result, :status).status)


        # count the number of the estimation result out of [0, 1]
        estimation_result  = dropmissing(estimation_result, :θ);
        estimation_result = filter(row -> (row.status_indicator == 1), estimation_result)

        number_non_missing = size(estimation_result.θ,1)
        number_out_range   = number_non_missing  - count(x -> (-10e-9 <= x <= 1 + 1e-8), estimation_result.θ)
        rate_out_range     = round(number_out_range/number_non_missing * 100, digits = 3)
        rate_out_range     = "$rate_out_range %"

        estimation_result = filter(x -> (-10^5 <= x <= 10^3), estimation_result.θ)


        if estimation_method[2] == :log_constraint
            if estimation_method[3] == :theta_constraint
                label_title = "both constraints ($rate_out_range are out of [0, 1])"
            else
                label_title = "only log constraint ($rate_out_range are out of [0, 1])"
            end
        else
            if estimation_method[3] == :theta_constraint
                label_title = "with constraint ($rate_out_range are out of [0, 1])"
            else
                label_title = "no constraint ($rate_out_range are out of [0, 1])"
            end
        end

        histogram!(histo_result, estimation_result, xlims = [-10^5, 10^3], xlabel = "θ",ylabel = "counts", label = label_title, fill = true, fillalpha = 0.5, bins = -10^5:1000:10^3)

        filename_begin = "../conduct_parameter/figuretable/tight/histogram_loglinear_loglinear_n_"
        filename_end   = "_random_start.pdf"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*"_"*String(estimation_method[3])filename_end

        display(histo_result)

        savefig(histo_result, file_name)


    end
end





#--------------------------------------------------------------------------------------------------------

# Histograms for the simulation results where the starting values are the true parameters and the GMM minimization problem is solved with a loose tolerance
for t = [50, 100, 200, 1000], sigma =  [0.001, 0.5, 1, 2]

    @unpack θ = parameter

    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end
    
    histo_result = histogram(xlims = [-10^5, 10^3], title = "true start, loose, n = $t, σ = $sigma", legend = :topright, size = (800, 600))
    vline!(histo_result, [θ], label = "true value : θ = $θ")

    for estimation_method = estimation_methods
        
        # Load the estimation result
        filename_estimation = "_"*String(estimation_method[1])*"_"*String(estimation_method[2])*"_"*String(estimation_method[3])
        filename_begin = "../conduct_parameter/output/loose/parameter_hat_table_loglinear_loglinear_n_"
        filename_end   = "_true_start_loose.csv"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end
        estimation_result = DataFrame(CSV.File(file_name))

        histo_status = plot(histogram(estimation_result.status), xlims = [1,24], bins = 1:1:24, title = "true start, loose, n = $t, σ = $sigma")

        display(histo_status)
        rate_convergence = count( x -> (x <= 7), dropmissing(estimation_result, :status).status)/length(dropmissing(estimation_result, :status).status)


        # count the number of the estimation result out of [0, 1]
        estimation_result  = dropmissing(estimation_result, :θ);
        estimation_result = filter(row -> (row.status_indicator == 1), estimation_result)

        number_non_missing = size(estimation_result.θ,1)
        number_out_range   = number_non_missing  - count(x -> (-10e-9 <= x <= 1 + 1e-8), estimation_result.θ)
        rate_out_range     = round(number_out_range/number_non_missing * 100, digits = 3)
        rate_out_range     = "$rate_out_range %"
        
        estimation_result = filter(x -> (-10^5 <= x <= 10^3), estimation_result.θ)


        if estimation_method[2] == :log_constraint
            if estimation_method[3] == :theta_constraint
                label_title = "both constraints ($rate_out_range are out of [0, 1])"
            else
                label_title = "only log constraint ($rate_out_range are out of [0, 1])"
            end
        else
            if estimation_method[3] == :theta_constraint
                label_title = "with constraint ($rate_out_range are out of [0, 1])"
            else
                label_title = "no constraint ($rate_out_range are out of [0, 1])"
            end
        end

        histogram!(histo_result, estimation_result, xlims = [-10^5, 10^3], xlabel = "θ",ylabel = "counts", label = label_title, fill = true, fillalpha = 0.5, bins = -10^5:1000:10^3)

        filename_begin = "../conduct_parameter/figuretable/loose/histogram_loglinear_loglinear_n_"
        filename_end   = "_true_start.pdf"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*"_"*String(estimation_method[3])filename_end

        display(histo_result)

        savefig(histo_result, file_name)


    end
end

# Histograms for the simulation results where the starting values are the randomly drawn and the GMM minimization problem is solved with a loose tolerance
for t = [50, 100, 200, 1000], sigma =  [0.001, 0.5, 1, 2]

    @unpack θ = parameter

    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end
    
    histo_result = histogram(xlims = [-10^5, 10^3], title = "random start, loose, n = $t, σ = $sigma", legend = :topright, size = (800, 600))
    vline!(histo_result, [θ], label = "true value : θ = $θ")

    for estimation_method = estimation_methods
        
        # Load the estimation result
        filename_estimation = "_"*String(estimation_method[1])*"_"*String(estimation_method[2])*"_"*String(estimation_method[3])
        filename_begin = "../conduct_parameter/output/loose/parameter_hat_table_loglinear_loglinear_n_"
        filename_end   = "_random_start_loose.csv"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end
        estimation_result = DataFrame(CSV.File(file_name))

        histo_status = plot(histogram(estimation_result.status), xlims = [1,24], bins = 1:1:24, title = "random start, loose, n = $t, σ = $sigma")

        display(histo_status)
        rate_convergence = count( x -> (x <= 7), dropmissing(estimation_result, :status).status)/length(dropmissing(estimation_result, :status).status)


        # count the number of the estimation result out of [0, 1]
        estimation_result  = dropmissing(estimation_result, :θ);
        estimation_result = filter(row -> (row.status_indicator == 1), estimation_result)

        number_non_missing = size(estimation_result.θ,1)
        number_out_range   = number_non_missing  - count(x -> (-10e-9 <= x <= 1 + 1e-8), estimation_result.θ)
        rate_out_range     = round(number_out_range/number_non_missing * 100, digits = 3)
        rate_out_range     = "$rate_out_range %"
        
        estimation_result = filter(x -> (-10^5 <= x <= 10^3), estimation_result.θ)


        if estimation_method[2] == :log_constraint
            if estimation_method[3] == :theta_constraint
                label_title = "both constraints ($rate_out_range are out of [0, 1])"
            else
                label_title = "only log constraint ($rate_out_range are out of [0, 1])"
            end
        else
            if estimation_method[3] == :theta_constraint
                label_title = "with constraint ($rate_out_range are out of [0, 1])"
            else
                label_title = "no constraint ($rate_out_range are out of [0, 1])"
            end
        end

        histogram!(histo_result, estimation_result, xlims = [-10^5, 10^3], xlabel = "θ",ylabel = "counts", label = label_title, fill = true, fillalpha = 0.5, bins = -10^5:1000:10^3)

        filename_begin = "../conduct_parameter/figuretable/loose/histogram_loglinear_loglinear_n_"
        filename_end   = "_random_start_loose.pdf"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*"_"*String(estimation_method[3])filename_end

        display(histo_result)

        savefig(histo_result, file_name)
    end
end