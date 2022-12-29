using LinearAlgebra, Distributions
using Statistics, Random, MultivariateStats
using JuMP, Ipopt
using DelimitedFiles, JLD, CSV, DataFrames, RData
using Plots, Combinatorics, Dates, StatsPlots
using Parameters: @unpack, @with_kw



market_parameters_log = @with_kw (
    α_0 = 10, # Demand parameter
    α_1 = 1,
    α_2 = 0.1,
    α_3 = 1,
    γ_0 = 1,  # Marginal cost parameter
    γ_1 = 1,
    γ_2 = 1,
    γ_3 = 1,
    θ = 0.3,  # Conduct paramter
    σ = 1,    # Standard deviation of the error term
    T = 50,   # Number of markets
    S = 1000, # Number of simulation
)

parameter = market_parameters_log()


estimation_methods = [(:separate,:log_constraint, :theta_constraint), (:separate,:log_constraint, :non_constraint), (:separate,:non_constraint, :theta_constraint), (:separate,:non_constraint, :non_constraint)];

# Check the summary of the eatimation result
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



        filename_estimation = "_"*String(estimation_method[1])*"_"*String(estimation_method[2])*"_"*String(estimation_method[3])
        filename_begin = "../conduct_parameter/output/parameter_hat_table_loglinear_loglinear_n_"
        filename_end   = ".csv"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end
        estimation_result = DataFrame(CSV.File(file_name))

        parameter = market_parameters_log(T = t, σ = sigma)

        @unpack θ, S = parameter

        rate_satisfy_assumption = 1 - sum(sum(1 .- θ .*(estimation_result.α_1[s] .+ estimation_result.α_2[s] .* data.z[(s-1)*t+1:s*t,:]) .<= 0)  for s = 1:S)/(S*t)

        #@show rate_satisfy_assumption

        display(describe(estimation_result))

    end
    println("------------------------------------------------------------------------------------\n")
end

# Draw histograms consisting of the estimation reuslt of θ within [-2, 3]

histogram_result_theta = []

@unpack θ = parameter

for t = [50, 100, 200, 1000], sigma =  [0.001, 0.5, 1, 2]


    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end

    histo_result = density(xlims = [-2, 3],  title = " n = $t, σ = $sigma", legend = :topright, size = (800, 600))
    
    #histo_result = histogram(xlims = [-2, 3], title = " n = $t, σ = $sigma", legend = :topright, size = (800, 600))
    
    vline!([θ], label = "true value : θ = $θ")

    for estimation_method = estimation_methods

        # Load the simulation data from the rds files
        filename_begin = "../conduct_parameter/output/data_loglinear_loglinear_n_"
        filename_end   = ".rds"



        filename_estimation = "_"*String(estimation_method[1])*"_"*String(estimation_method[2])*"_"*String(estimation_method[3])
        filename_begin = "../conduct_parameter/output/parameter_hat_table_loglinear_loglinear_n_"
        filename_end   = ".csv"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end
        estimation_result = DataFrame(CSV.File(file_name))

        # count the number of the estimation result out of [-2, 3]
        estimation_result  = dropmissing(estimation_result, :θ);
        number_non_missing = size(estimation_result.θ,1)
        number_out_range   = number_non_missing  - count(x -> (-2 <= x <= 3), estimation_result.θ)
        rate_out_range     = "$number_out_range/$number_non_missing"
        
        estimation_result = filter(x -> (-2 <= x <= 3), estimation_result.θ)

        if estimation_method[2] == :log_constraint
            if estimation_method[3] == :theta_constraint
                label_title = "both constraints ($rate_out_range are out of [-2,3])"
            else
                label_title = "only log constraint ($rate_out_range are out of [-2,3])"
            end
        else
            if estimation_method[3] == :theta_constraint
                label_title = "only theta constraint ($rate_out_range are out of [-2,3])"

            else
                label_title = "no constraint ($rate_out_range are out of [-2,3])"
            end
        end

        #histogram!(histo_result, estimation_result, xlims = [-2, 3], xlabel = "θ",ylabel = "counts", label = label_title, fill = true, fillalpha = 0.5)


        density!(histo_result, estimation_result, xlims = [-2, 3], xlabel = "θ",ylabel = "counts", label = label_title, alpha = 0.5)

    end

    push!(histogram_result_theta, histo_result)

    filename_begin = "../conduct_parameter/figuretable/histogram_loglinear_loglinear_n_"
    filename_end   = ".png"
    file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end

    savefig(histo_result, file_name)
end
   
# Put all histograms in one plot
#plot(histogram_result_theta...,  layout=(4,4), legend= true, size = (3000,3000))
