include(joinpath(dirname(@FILE),"00setting_julia.jl"))
#-----------------------------------------------------------------------------------------
# Draw histograms of the estimation reuslt of θ 

for t = [50, 100, 200, 1000], sigma =  [0.001, 0.5, 1, 2]

    @unpack θ = parameter

    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end
    
    histo_result = histogram(xlims = [0, 1], title = " n = $t, σ = $sigma", legend = :topright, size = (800, 600))
    vline!(histo_result, [θ], label = "true value : θ = $θ")

    for estimation_method = estimation_methods
        
        # Load the estimation result
        filename_estimation = "_"*String(estimation_method[1])*"_"*String(estimation_method[2])*"_"*String(estimation_method[3])
        filename_begin = "../conduct_parameter/output/parameter_hat_table_loglinear_loglinear_n_"
        filename_end   = ".csv"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end
        estimation_result = DataFrame(CSV.File(file_name))

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

        filename_begin = "../conduct_parameter/figuretable/histogram_loglinear_loglinear_n_"
        filename_end   = ".pdf"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*"_"*String(estimation_method[3])filename_end

        savefig(histo_result, file_name)
    end
end


#-----------------------------------------------------------------------------------------
# Draw the contour plot of the GMM value function with respect to the conduct parameter and the constant in the marginal cost function
function contour_set_of_GMM(parameter, data, theta_range, gamma_range)
    

    @unpack α_0, α_1, α_2, α_3,γ_0 , γ_1 ,γ_2 ,γ_3, θ, σ ,T = parameter

    γ = [γ_0 , γ_1 ,γ_2 ,γ_3]
    α = [α_0, α_1, α_2, α_3]

    Q  = data.logQ
    w  = data.w
    r  = data.r
    z  = data.z
    iv_w = data.iv_w
    iv_r = data.iv_r
    p  = data.logP
    y  = data.y

    iv = hcat(iv_w, iv_r)

    X_d = []
    X_s = []
    X   = []
    Z   = []
    Z_d = []
    Z_s = []
    P   = []

    for t = 1:T

        Z_dt = vcat(1, z[t], iv[t,:], y[t])
        Z_st = vcat(1, z[t], log(w[t]), log(r[t]), y[t])
        Z_t = [Z_dt zeros(length(Z_dt));  zeros(length(Z_st)) Z_st]'

        X_dt = vcat(1, -Q[t], -z[t].*Q[t], y[t])
        X_st = vcat(1, Q[t], log(w[t]), log(r[t]), z[t])
        X_t  = [X_dt zeros(length(X_dt));  zeros(length(X_st)) X_st]'

        push!(P, p[t])
        push!(X_d, X_dt')
        push!(X_s, X_st')
        push!(X, X_t)
        push!(Z, Z_t)
        push!(Z_d, Z_dt')
        push!(Z_s, Z_st')
    end

    Z   = reduce(vcat,(Z))
    X   = reduce(vcat,(X))
    X_d = reduce(vcat,(X_d))
    X_s = reduce(vcat,(X_s))
    Z_d = reduce(vcat,(Z_d))
    Z_s = reduce(vcat,(Z_s))

    L = size(Z, 2)
    L_d = size(Z_d,2)
    L_s = size(Z_s,2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    Ω = inv(Z_s' * Z_s)/T

    GMM_value = []


    for θ_hat = theta_range, γ_hat = gamma_range

        r = P .- γ_hat .-sum(γ[k] .* X_s[:,k] for k = 2:K_s-1) + log.(1 .- θ_hat .*(α_1 .+ α_2 .* X_s[:, end]) ) 

        g = sum(Z_s[t,:] .* r[t] for t = 1:T)

        gmm_value = sum( g[l] *Ω[l,k] * g[k] for l = 1:L_s, k = 1:L_s)

        push!(GMM_value, gmm_value)
    end

    return reshape(GMM_value, (length(gamma_range), length(theta_range)))
end

# Draw the contour figure for each simulation setting
for t = [50, 100, 200, 1000], sigma =  [0.001, 0.5, 1, 2]
    @unpack θ, γ_0 = parameter

    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end
    
    # Load the simulation data from the rds files
    filename_begin = "../conduct_parameter/output/data_loglinear_loglinear_n_"
    filename_end   = ".rds"

    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end

    filename = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end

    data = load(filename)
    data = DataFrames.sort(data, [:group_id_k])
    data = data[1:t,:]

    theta_range = [-2:0.01:0.7;] 
    gamma_range = [-10:0.01:10;]

    contour_gmm = contour_set_of_GMM(parameter, data, theta_range, gamma_range);
    plot_contour = plot(contour(theta_range, gamma_range, contour_gmm,
        xlabel="θ", ylabel="γ_0",
        title="N =$t, σ = $sigma"))
        vline!([θ], linestyle=:dash, label = "true θ")
        hline!([γ_0], linestyle=:dash, label = "true γ_0")

    filename_begin = "../conduct_parameter/figuretable/contour_loglinear_loglinear_n_"
    filename_end   = ".pdf"
    file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end

    savefig(plot_contour, file_name)
end



#-----------------------------------------------------------------------------------------

# Draw the figure of the difference of the values of the GMM objective function under the true parameters and the estimated parameter  

function value_GMM(parameter, data, estimation_result)
    
    @unpack α_0, α_1, α_2, α_3, γ_0, γ_1 ,γ_2 ,γ_3, θ, σ ,T = parameter

    γ = [γ_0 , γ_1 ,γ_2 ,γ_3]
    α = [α_0, α_1, α_2, α_3]

    Q  = data.logQ
    w  = data.w
    r  = data.r
    z  = data.z
    iv_w = data.iv_w
    iv_r = data.iv_r
    p  = data.logP
    y  = data.y

    iv = hcat(iv_w, iv_r)

    X_d = []
    X_s = []
    X   = []
    Z   = []
    Z_d = []
    Z_s = []
    P   = []

    for t = 1:T

        Z_dt = vcat(1, z[t], iv[t,:], y[t])
        Z_st = vcat(1, z[t], log(w[t]), log(r[t]), y[t])
        Z_t = [Z_dt zeros(length(Z_dt));  zeros(length(Z_st)) Z_st]'

        X_dt = vcat(1, -Q[t], -z[t].*Q[t], y[t])
        X_st = vcat(1, Q[t], log(w[t]), log(r[t]), z[t])
        X_t  = [X_dt zeros(length(X_dt));  zeros(length(X_st)) X_st]'

        push!(P, p[t])
        push!(X_d, X_dt')
        push!(X_s, X_st')
        push!(X, X_t)
        push!(Z, Z_t)
        push!(Z_d, Z_dt')
        push!(Z_s, Z_st')
    end

    Z   = reduce(vcat,(Z))
    X   = reduce(vcat,(X))
    X_d = reduce(vcat,(X_d))
    X_s = reduce(vcat,(X_s))
    Z_d = reduce(vcat,(Z_d))
    Z_s = reduce(vcat,(Z_s))

    L = size(Z, 2)
    L_d = size(Z_d,2)
    L_s = size(Z_s,2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    Ω = inv(Z_s' * Z_s)/T
    θ_hat = estimation_result.θ
    γ_hat = estimation_result.γ_0


    
    if 1 <= sum(1 .- θ_hat .*(α_1 .+ α_2 .* X_s[:, end]) .<= 0)
        return missing, missing, missing
    else
        
        r_true = P .- γ_0 .-sum(γ[k] .* X_s[:,k] for k = 2:K_s-1) + log.(1 .- θ .*(α_1 .+ α_2 .* X_s[:, end]) ) 
        r_estimated = P .- γ_hat .-sum(γ[k] .* X_s[:,k] for k = 2:K_s-1) + log.(1 .- θ_hat .*(α_1 .+ α_2 .* X_s[:, end]) ) 

        g_true      = sum(Z_s[t,:] .* r_true[t] for t = 1:T)
        g_estimated = sum(Z_s[t,:] .* r_estimated[t] for t = 1:T)

        gmm_value_true = sum( g_true[l] *Ω[l,k] * g_true[k] for l = 1:L_s, k = 1:L_s)
        gmm_value_estimated = sum( g_estimated[l] *Ω[l,k] * g_estimated[k] for l = 1:L_s, k = 1:L_s)

        difference_theta = abs(θ_hat - θ)
        difference_gamma = abs(γ_hat - γ_0)
        difference_gmm_value = abs(gmm_value_true - gmm_value_estimated)

        return difference_theta, difference_gamma, difference_gmm_value
    end
end


# Draw the picture for each simulation setting
for t = [50, 100, 200, 1000], sigma =  [0.001, 0.5, 1, 2]
    @unpack θ, γ_0, S = parameter

    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end
    
    # Load the simulation data from the rds files
    filename_begin = "../conduct_parameter/output/data_loglinear_loglinear_n_"
    filename_end   = ".rds"

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
        filename_end   = ".csv"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end
        estimation_result = DataFrame(CSV.File(file_name))

        if estimation_method[3] == :theta_constraint
            label_title = "estimation with constraint"

        else
            label_title = "estimation without constraint"
        end


        number_sample = size(estimation_result,1)

        for s = 1:number_sample
            estimation_result_s = estimation_result[s,:]
            if estimation_result_s.θ !== missing && estimation_result_s.status_indicator == 1
                data_s = data[(s-1)*t+1:s*t,:]
                diff_theta, diff_gamma, diff_gmm_value = value_GMM(parameter, data_s, estimation_result_s);

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
        filename_end   = ".pdf"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end
    

        df = DataFrame(theta = difference_theta, gamma = difference_gamma, gmm_value = difference_gmm_value)
        df |> @vlplot(:circle, 
            x={:theta, title = "|Estimated θ - true θ|"}, 
            y={:gamma, title = "|Estimated γ_0 - true γ_0|"}, 
            color = {:gmm_value, title = "diff. GMM value", scale={scheme=:plasma}},
            title = "n = $t, σ = $sigma, $constraint", xtitle = "difference"
        )|> save(file_name)
    
    end
end