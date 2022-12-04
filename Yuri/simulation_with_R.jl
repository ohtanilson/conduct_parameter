using LinearAlgebra, Distributions
using Statistics, Random, MultivariateStats
using DataFrames, RData




function TwoSLS_estimation_R(data)

    T = size(data, 1)

    Q = data.Q
    W = data.w
    R = data.r
    Z = data.z
    P = data.P
    iv_w = data.iv_w
    iv_r = data.iv_r
    Y = data.y
    
    Z_demand = hcat(ones(T), Z, iv_w, iv_r, Y)
    Q_hat = Z_demand * pinv(Z_demand' * Z_demand)* (Z_demand' * Q)
    X_d = hcat(ones(T), Q_hat, Z .* Q_hat, Y)
    α_hat = pinv(X_d' * X_d) * (X_d' * P)

    #=    
    data_1st = DataFrame(Q = Q, Z = Z, IV_w = iv_w, IV_r = iv_r, Y = Y)
    iv_1st = lm(@formula(Q ~ Z + IV_w + IV_r + Y), data_1st, dropcollinear=false)

    Q_hat = predict(iv_1st)

    data_2nd = DataFrame(P = P, Q_hat = -Q_hat, Q_Z = -Z .* Q_hat, Y = Y)
    iv_2nd= lm(@formula(P ~ Q_hat + Q_Z + Y), data_2nd, dropcollinear=false)
    α_hat_GLM = coef(iv_2nd)
    =#

    # Supply side estimation 

    Z_supply = hcat(ones(T), W, R, Z, Y)
    Q_hat = Z_supply * pinv(Z_supply' * Z_supply) * Z_supply' *  Q
    X_s = hcat(ones(T), Q_hat, W, R, (α_hat[2] .+ α_hat[3] .* Z) .* Q_hat)
    γ_hat = pinv(X_s' * X_s) * (X_s' * P)

    θ_hat = γ_hat[end]

    return α_hat, γ_hat[1:end-1], θ_hat

end

function simulation_2sls_R(filename, sample_size, sigma)


    data = load(filename)
    data = DataFrames.sort(data, [:group_id_k])

    α_est = Vector{Float64}[]
    γ_est = Vector{Float64}[]
    θ_est = Float64[]

    for s = 1:1000

        data_s = data[(s-1)*sample_size+1:s*sample_size,:]

        α_est_s, γ_est_s, θ_est_s = TwoSLS_estimation_R(data_s)

        push!(α_est, α_est_s)
        push!(γ_est, γ_est_s)
        push!(θ_est, θ_est_s)

    end

    α_est = reduce(vcat, α_est')
    γ_est = reduce(vcat, γ_est')
    θ_est = reduce(vcat, θ_est')

    result = DataFrame(
        σ = sigma,
        T = sample_size,
        α_0 = α_est[:,1],
        α_1 = α_est[:,2],
        α_2 = α_est[:,3],
        α_3 = α_est[:,end],
        γ_0 = γ_est[:,1],
        γ_1 = γ_est[:,2],
        γ_2 = γ_est[:,3],
        γ_3 = γ_est[:,4],
        θ = θ_est
    )

    return describe(result, :mean, :std, :min, :max, :median)

end



for sigma =  [0.001, 0.5, 1, 2] , sample_size = [50, 100, 200, 1000]
    filename_begin = "../conduct_parameter/R/output/data_linear_linear_n_"
    filename_end   = ".rds"

    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end

    filename = filename_begin*string(sample_size)*"_sigma_"*string(sigma)*filename_end

    result = simulation_2sls_R(filename, sample_size, sigma)

    display(result)
end