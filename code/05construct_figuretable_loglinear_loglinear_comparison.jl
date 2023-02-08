include("00setting_julia.jl")
include("00functions.jl")
parameter = market_parameters_log()
estimation_methods = [(:separate,:non_constraint, :non_constraint), (:separate,:non_constraint, :theta_constraint)];

#-----------------------------------------------------------------------------------------
# Draw histograms of the estimation reuslt of θ 
#-----------------------------------------------------------------------------------------
for estimation_method = estimation_methods
    for t = [50, 100, 200, 1000], sigma = [0.001, 0.5, 1, 2]

        @unpack θ = parameter
        if sigma == 1 || sigma == 2
            sigma = Int64(sigma)
        end

        histo_result = Plots.histogram(xlims = [0, 1], title = " T = $t, σ = $sigma", legend = :topright, size = (800, 600))
        Plots.vline!(histo_result, [θ], label = "true value : θ = $θ")


            
            # Load the estimation result
            filename_estimation = "_"*String(estimation_method[1])*"_"*String(estimation_method[2])*"_"*String(estimation_method[3])
            filename_begin = "../conduct_parameter/output/parameter_hat_table_loglinear_loglinear_n_"
            filename_end = ".csv"
            file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end
            estimation_result = DataFrame(CSV.File(file_name))

            # count the number of the estimation result out of [0, 1]
            estimation_result = dropmissing(estimation_result, :θ);
            estimation_result = filter(row -> (row.status_indicator == 1), estimation_result)

            @show minimum(estimation_result.θ)

            number_non_missing = size(estimation_result.θ, 1)
            number_out_range = number_non_missing  - count(x -> (-10e-9 <= x <= 1 + 1e-8), estimation_result.θ)
            rate_out_range = round(number_out_range/number_non_missing * 100, digits = 3)
            rate_out_range = "$rate_out_range %"
            
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

            Plots.histogram!(histo_result, estimation_result, xlims = [-2, 3], xlabel = "θ",ylabel = "counts", label = label_title, fill = true, fillalpha = 0.5, bins = -2:0.1:3)

            filename_begin = "../conduct_parameter/figuretable/histogram_loglinear_loglinear_n_"
            filename_end = ".pdf"
            file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*"_"*String(estimation_method[3])filename_end

            #Plots.savefig(histo_result, file_name)
    end
    println("----------------------------------------------------------------")
end
#-----------------------------------------------------------------------------------------
# Draw the contour figure for each simulation setting
#-----------------------------------------------------------------------------------------
for t = [50, 100, 200, 1000], sigma = [0.001, 0.5, 1, 2]
    @unpack θ, γ_0 = parameter

    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end
    
    # Load the simulation data from the rds files
    filename_begin = "../conduct_parameter/output/data_loglinear_loglinear_n_"
    filename_end = ".rds"

    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end

    filename = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end

    data = load(filename)
    data = DataFrames.sort(data, [:group_id_k])
    data = data[1:t,:]

    theta_range = [-2:0.01:0.7;] 
    gamma_range = [-10:0.01:10;]

    contour_gmm = generate_contour_set_of_GMM(parameter, data, theta_range, gamma_range);
    plot_contour = Plots.plot(
        contour(theta_range, gamma_range, contour_gmm,
        xlabel="θ", ylabel="γ_0",
        title="T =$t, σ = $sigma"))
        vline!([θ], linestyle=:dash, label = "true θ")
        hline!([γ_0], linestyle=:dash, label = "true γ_0")

    filename_begin = "../conduct_parameter/figuretable/contour_loglinear_loglinear_n_"
    filename_end = ".pdf"
    file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end

    Plots.savefig(plot_contour, file_name)
end

#-----------------------------------------------------------------------------------------
# Draw the picture for each simulation setting
#-----------------------------------------------------------------------------------------
for t = [50, 100, 200, 1000], sigma = [0.001, 0.5, 1, 2]
    @unpack θ, γ_0, S = parameter

    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end
    
    # Load the simulation data from the rds files
    filename_begin = "../conduct_parameter/output/data_loglinear_loglinear_n_"
    filename_end = ".rds"

    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end

    filename = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end
    data = load(filename)
    data = DataFrames.sort(data, [:group_id_k])

    for estimation_method = estimation_methods

        difference_theta = Union{Missing, Float64}[]
        difference_gamma = Union{Missing, Float64}[]
        difference_gmm_value = Union{Missing, Float64}[]

        # Load the estimation result
        filename_estimation = "_"*String(estimation_method[1])*"_"*String(estimation_method[2])*"_"*String(estimation_method[3])
        filename_begin = "../conduct_parameter/output/parameter_hat_table_loglinear_loglinear_n_"
        filename_end = ".csv"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end
        estimation_result = DataFrame(CSV.File(file_name))

        if estimation_method[3] == :theta_constraint
            label_title = "estimation with constraint"

        else
            label_title = "estimation without constraint"
        end


        number_sample = size(estimation_result, 1)

        for s = 1:number_sample
            estimation_result_s = estimation_result[s, :]
            if estimation_result_s.θ !== missing && estimation_result_s.status_indicator == 1
                data_s = data[(s-1)*t+1:s*t, :]
                diff_theta, diff_gamma, diff_gmm_value = generate_difference_GMM_value(parameter, data_s, estimation_result_s);

                push!(difference_theta, diff_theta)
                push!(difference_gamma, diff_gamma)
                push!(difference_gmm_value, diff_gmm_value)
            end
        end

        if  estimation_method[3] == :theta_constraint
            constraint = "with constraint"
        else
            constraint = "without constraint"
        end

        filename_estimation = "_"*String(estimation_method[3])
        filename_begin = "../conduct_parameter/figuretable/diff_gmm_value_loglinear_loglinear_n_"
        filename_end = ".pdf"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end
    

        df = DataFrame(theta = difference_theta, gamma = difference_gamma, gmm_value = difference_gmm_value)
        df |> VegaLite.@vlplot(:circle, 
            x={:theta, title = "|Estimated θ - true θ|"}, 
            y={:gamma, title = "|Estimated γ_0 - true γ_0|"}, 
            color = {:gmm_value, title = "diff. GMM value", scale={scheme=:plasma}},
            title = "T = $t, σ = $sigma, $constraint", xtitle = "difference"
        )|> save(file_name)
    
    end
end