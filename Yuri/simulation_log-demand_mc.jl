using LinearAlgebra, Distributions
using Statistics, Random, MultivariateStats
using JuMP, Ipopt
using DelimitedFiles, JLD, CSV, DataFrames
using Plots, Combinatorics, Dates, StatsPlots
using Parameters: @unpack, @with_kw


#---------------------------------------------------------------------------------------------------------

## Replication of Hyde and Perloff (1995)

#---------------------------------------------------------------------------------------------------------

market_parameters_log = @with_kw (
    θ = 0.5,
    σ = 1,
    T = 200,
    S = 1000,
    α = 1,
    β = 2,
    γ = 3,
    A = 1.2,
    α_0 = 1.8,
    α_1 = 1.2,
    α_2 = -0.5,
    α_3 = 1.0,
    γ_0 = -1/γ * log(A) - α/γ * log(α) - β/γ * log(β),
    γ_1 = (1- γ)/γ,
    γ_2 = α/γ,
    γ_3 = β/γ,

)

mutable struct market_data_log
    Q::Vector{Float64}
    W::Vector{Float64}
    R::Vector{Float64}
    Z::Vector{Float64}
    P::Vector{Float64}
    Y::Vector{Float64}
    MC::Vector{Float64}
    IV::Matrix{Float64}
end


function simulation_data_log(parameter,s)


    Random.seed!(s*1234)

    @unpack α_0, α_1, α_2, α_3 ,γ_0, γ_1, γ_2, γ_3, θ, σ, T, α, β, γ, A = parameter
    
    Q = Float64[];
    P = Float64[];
    MC = Float64[];
    W  = Float64[];
    R  = Float64[];
    Z = Float64[];   
    Y = Float64[];
    IV = zeros(T,2);

    for t = 1:T
        r_t = rand(Uniform(0,1))
        z_t = rand(Uniform(5,10))
        w_t = rand(Uniform(1,3))
        iv_w_t = w_t + randn()
        iv_r_t = r_t + randn()

        y_t = randn()
        ε_c = rand(Normal(0,σ))
        ε_d = rand(Normal(0,σ))


        log_Q_t = (α_0 + log(1 - θ * (α_1 + α_2 * z_t))  + α_3 * y_t - ( γ_0 + γ_2 * log(w_t) + γ_3 * log(r_t)) + ε_d- ε_c)/ (γ_1 + α_1 + α_2 * z_t) # Equilibrium total quantity

        log_p_t = α_0 - (α_1 + α_2 * z_t )* log_Q_t  + α_3 * y_t + ε_d # The demand function

        push!(Q, log_Q_t)
        push!(P,log_p_t)
        push!(W,log(w_t))
        push!(R,log(r_t))
        push!(Z,z_t)
        push!(Y,y_t)
        IV[t,:] .= (iv_w_t, iv_r_t)

    end

    data = market_data_log(Q, W, R, Z, P,Y, MC, IV)
    return data
end

function GMM_estimation(T, P, Z, X, X_s, X_d, Ω)
    
    L = size(Z, 2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "tol", 1e-15)
    set_optimizer_attribute(model, "max_iter", 1000)
    set_optimizer_attribute(model, "acceptable_tol", 1e-12)
    set_silent(model)
    @variable(model, β[k = 1:K_d+K_s-1])
    @variable(model, θ)

    r = Any[];
    g = Any[];
    for t =1:T
        push!(r, @NLexpression(model, P[t]- sum(β[k] * X[2*t-1,k] for k = 1:K_d) ))
        push!(r, @NLexpression(model, P[t]- sum(β[k] * X[2*t,k] for k = K_d+1:K_d+K_s-1) + log(1 - θ *(β[2] + β[3] * X[2*t, end]) ) ) )
    end

    for l = 1:L
        push!(g, @NLexpression(model, sum(Z[t,l] * r[t] for t = 1:2*T)))
    end

    @NLobjective(model, Min, sum( g[l] *Ω[l,k] * g[k] for l = 1:L, k = 1:L))

    optimize!(model)
    
    α_hat = value.(β)[1:K_d]
    γ_hat = value.(β)[K_d+1:end]
    θ_hat = value.(θ)

    return α_hat, γ_hat, θ_hat
end


function nonlinear_TwoSLS(parameter, data)
    @unpack T = parameter

    Q  = data.Q
    w  = data.W
    r  = data.R
    z  = data.Z
    iv = data.IV
    p  = data.P
    y  = data.Y

    X_d = []
    X_s = []
    X   = []
    Z   = []
    Z_d = []
    Z_s = []
    P   = []

    for t = 1:T

        Z_dt = vcat(1, z[t], w[t], r[t], y[t])
        Z_st = vcat(1, z[t], iv[t,:], y[t])
        Z_t = [Z_dt zeros(length(Z_dt));  zeros(length(Z_st)) Z_st]'

        X_dt = vcat(1, -Q[t], -z[t].*Q[t], y[t])
        X_st = vcat(1, Q[t], w[t], r[t], z[t])
        X_t  = [X_dt zeros(length(X_dt));  zeros(length(X_st)) X_st]'

        push!(P, p[t])
        push!(X_d, X_dt')
        push!(X_s, X_st')
        push!(X, X_t)
        push!(Z, Z_t)
        push!(Z_d, Z_dt)
        push!(Z_s, Z_st)
    end

    Z   = reduce(vcat,(Z))
    X   = reduce(vcat,(X))
    X_d = reduce(vcat,(X_d))
    X_s = reduce(vcat,(X_s))
    Z_d = reduce(vcat,(Z_d))
    Z_s = reduce(vcat,(Z_s))

    Ω = inv(Z'Z/T)

    α_hat, γ_hat, θ_hat = GMM_estimation(T, P, Z, X, X_s, X_d, Ω)

    return α_hat, γ_hat, θ_hat
end


function simulation_nonlinear_2sls()

    Results = []

    for t in [50, 100, 200, 1000]
        
        parameter = market_parameters_log(T = t)

        @unpack S = parameter 

        α_est = Vector{Float64}[]
        γ_est = Vector{Float64}[]
        θ_est = Float64[]

        for s = 1:S
            data_log_s = simulation_data_log(parameter, s);
            α_est_s, γ_est_s, θ_est_s = nonlinear_TwoSLS(parameter, data_log_s)
        
            push!(α_est, α_est_s)
            push!(γ_est, γ_est_s)
            push!(θ_est, θ_est_s)
        end
        
        α_est = reduce(vcat, α_est')
        γ_est = reduce(vcat, γ_est')
        θ_est = reduce(vcat, θ_est')

        result = DataFrame(
        T = t,
        α_0 = α_est[:,1],
        α_1 = α_est[:,2],
        α_2 = α_est[:,3],
        α_3 = α_est[:,4],
        γ_0 = γ_est[:,1],
        γ_1 = γ_est[:,2],
        γ_2 = γ_est[:,3],
        γ_3 = γ_est[:,4],
        θ = θ_est)

        push!(Results, describe(result, :mean, :std, :min, :max, :median))

    end

    return Results
end

nonlinear_2sls_result = simulation_nonlinear_2sls()



#-------------------------------------------------------------------------------------------


function nonlinear_ThreeSLS(parameter, data)

    @unpack T = parameter

    Q  = data.Q
    w  = data.W
    r  = data.R
    z  = data.Z
    iv = data.IV
    p  = data.P
    y  = data.Y

    X_d = []
    X_s = []
    X   = []
    Z   = []
    Z_d = []
    Z_s = []
    P   = []

    for t = 1:T

        Z_st = vcat(1, z[t], iv[t,:], y[t])
        Z_dt = vcat(1, z[t], w[t], r[t], y[t])

        Z_t = [Z_dt zeros(length(Z_dt));  zeros(length(Z_st)) Z_st]'

        X_dt = vcat(1, -Q[t], -z[t].*Q[t], y[t])
        X_st = vcat(1, Q[t], w[t], r[t], z[t])
        X_t  = [X_dt zeros(length(X_dt));  zeros(length(X_st)) X_st]'

        push!(P, p[t])
        push!(X_d, X_dt')
        push!(X_s, X_st')
        push!(X, X_t)
        push!(Z, Z_t)
        push!(Z_d, Z_dt)
        push!(Z_s, Z_st)
    end

    Z   = reduce(vcat,(Z))
    X   = reduce(vcat,(X))
    X_d = reduce(vcat,(X_d))
    X_s = reduce(vcat,(X_s))
    Z_d = reduce(vcat,(Z_d))
    Z_s = reduce(vcat,(Z_s))
    

    L = size(Z, 2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    Ω = inv(Z'Z/T)

    α_2sls, γ_2sls, θ_2sls = GMM_estimation(T, P, Z, X, X_s, X_d, Ω)



    u_d = P .- sum(α_2sls[k] .* X_d[:,k] for k = 1:K_d)
    u_s = P .- sum(γ_2sls[k] .* X_s[:,k] for k = 1:K_s-1) .+ log.(1 .- θ_2sls*(α_2sls[2] .+ α_2sls[3] * X_s[:,end]))

    u_hat = reduce(vcat, [[u_d[t], u_s[t]] for t = 1:T])

    Ω = zeros(L,L)

    for t = 1:T
        Ω .+= Z[2*t-1:2*t,:]' * u_hat[2*t - 1:2*t] * u_hat[2*t - 1:2*t]' * Z[2*t-1:2*t,:]
    end

    Ω  = inv(Ω/T)

    α_hat, γ_hat, θ_hat =  GMM_estimation(T, P, Z, X, X_s, X_d, Ω)
    return α_hat, γ_hat, θ_hat
end


function simulation_nonlinear_3sls()

    Results = []

    for t in [50, 100, 200, 1000]
        
        parameter = market_parameters_log(T = t)

        @unpack S = parameter 

        α_est = Vector{Float64}[]
        γ_est = Vector{Float64}[]
        θ_est = Float64[]

        for s = 1:S
            data_log_s = simulation_data_log(parameter, s);
            α_est_s, γ_est_s, θ_est_s = nonlinear_ThreeSLS(parameter, data_log_s)
        
            push!(α_est, α_est_s)
            push!(γ_est, γ_est_s)
            push!(θ_est, θ_est_s)
        end
        
        α_est = reduce(vcat, α_est')
        γ_est = reduce(vcat, γ_est')
        θ_est = reduce(vcat, θ_est')

        result = DataFrame(
        T = t,
        α_0 = α_est[:,1],
        α_1 = α_est[:,2],
        α_2 = α_est[:,3],
        α_3 = α_est[:,4],
        γ_0 = γ_est[:,1],
        γ_1 = γ_est[:,2],
        γ_2 = γ_est[:,3],
        γ_3 = γ_est[:,4],
        θ = θ_est)

        push!(Results, describe(result, :mean, :std, :min, :max, :median))

    end

    return Results
end

nonlinear_3sls_result = simulation_nonlinear_3sls()






#nonlinear_3sls_result = nonlinear_ThreeSLS(parameter, data_log, α_2sls, γ_2sls, θ_2sls)






