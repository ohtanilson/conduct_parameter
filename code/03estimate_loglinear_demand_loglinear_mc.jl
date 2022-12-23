using LinearAlgebra, Distributions
using Statistics, Random, MultivariateStats
using JuMP, Ipopt
using DelimitedFiles, JLD, CSV, DataFrames, RData
using Plots, Combinatorics, Dates, StatsPlots
using Parameters: @unpack, @with_kw


# for sigma =  [0.001, 0.5, 1, 2] , sample_size = [50, 100, 200, 1000]
#     filename_begin = "../conduct_parameter/output/data_loglinear_loglinear_n_"
#     filename_end   = ".rds"

#     if sigma == 1 || sigma == 2
#         sigma = Int64(sigma)
#     end

#     filename = filename_begin*string(sample_size)*"_sigma_"*string(sigma)*filename_end
#     data = load(filename)
#     data = DataFrames.sort(data, [:group_id_k])

# end




function GMM_estimation_separate(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d)
    

    QZ_hat = Z_d * inv(Z_d' * Z_d) * Z_d' * (Z_d[:,2] .* Q)
    Q_hat = Z_d * inv(Z_d' * Z_d) * Z_d' *  Q
    X_dd = hcat(ones(T), -Q_hat, -QZ_hat, X_d[:,end])
    α_hat = inv(X_dd' * X_dd) * (X_dd' * P)


    L = size(Z, 2)
    L_d = size(Z_d,2)
    L_s = size(Z_s,2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    Ω = inv(Z_s' * Z_s)/T

    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "tol", 1e-15)
    set_optimizer_attribute(model, "max_iter", 1000)
    set_optimizer_attribute(model, "acceptable_tol", 1e-12)
    set_silent(model)
    @variable(model, γ[k = 1:K_s-1])
    @variable(model, 0 <= θ <= 1)

    r = Any[];
    g = Any[];
    for t =1:T
        push!(r, @NLexpression(model, P[t]- sum(γ[k] * X_s[t,k] for k = 1:K_s-1) + log(1 - θ *(α_hat[2] + α_hat[3] * X_s[t, end]) ) ) )
    end

    for l = 1:L_s
        push!(g, @NLexpression(model, sum(Z_s[t,l] * r[t] for t = 1:T)))
    end

    @NLobjective(model, Min, sum( g[l] *Ω[l,k] * g[k] for l = 1:L_s, k = 1:L_s))

    optimize!(model)
    
    γ_hat = value.(γ)
    θ_hat = value.(θ)

    return α_hat, γ_hat, θ_hat#, termination_status_code(termination_status(model))
end


function estimation_nonlinear_2SLS_separate(parameter, data)
    @unpack T = parameter

    Q  = data.Q
    w  = data.w
    r  = data.r
    z  = data.z
    iv_w = data.iv_w
    iv_r = data.iv_r
    p  = data.P
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
        Z_st = vcat(1, z[t], w[t], r[t], y[t])
        Z_t = [Z_dt zeros(length(Z_dt));  zeros(length(Z_st)) Z_st]'

        X_dt = vcat(1, -Q[t], -z[t].*Q[t], y[t])
        X_st = vcat(1, Q[t], w[t], r[t], z[t])
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


    #α_hat, γ_hat, θ_hat, status = GMM_estimation_separate(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d)
    α_hat, γ_hat, θ_hat = GMM_estimation_separate(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d)

    return α_hat, γ_hat, θ_hat, status
end



sigma = 0.001
sample_size = 50
filename_begin = "../conduct_parameter/output/data_with_demand_hat_loglinear_loglinear_n_"
filename_end   = ".rds"

if sigma == 1 || sigma == 2
    sigma = Int64(sigma)
end

filename = filename_begin*string(sample_size)*"_sigma_"*string(sigma)*filename_end
data = load(filename)
data = DataFrames.sort(data, [:group_id_k])

estimation_nonlinear_2SLS_separate(parameter, data)



# function simulation_nonlinear_2SLS_separate(parameter, data)

#     @unpack T, S, σ = parameter 

#     α_est = Vector{Float64}[]
#     γ_est = Vector{Float64}[]
#     θ_est = Float64[]
#     status = Int64[]

#     for s = 1:S
#         data_s = data[(s-1)*T+1:s*T,:]
#         α_est_s, γ_est_s, θ_est_s, status_s = estimation_nonlinear_2SLS_separate(parameter, data_s)
    
#         push!(α_est, α_est_s)
#         push!(γ_est, γ_est_s)
#         push!(θ_est, θ_est_s)
#         push!(status,status_s)
#     end
    
#     α_est = reduce(vcat, α_est')
#     γ_est = reduce(vcat, γ_est')
#     θ_est = reduce(vcat, θ_est')
#     status = reduce(vcat, status)

#     estimation_result = DataFrame(
#     T = T,
#     σ = σ,
#     α_0 = α_est[:,1],
#     α_1 = α_est[:,2],
#     α_2 = α_est[:,3],
#     α_3 = α_est[:,4],
#     γ_0 = γ_est[:,1],
#     γ_1 = γ_est[:,2],
#     γ_2 = γ_est[:,3],
#     γ_3 = γ_est[:,4],
#     θ = θ_est,
#     status = status)

#     return estimation_result

# end



