using LinearAlgebra, Distributions
using Statistics, Random, MultivariateStats
using JuMP, Ipopt
using DelimitedFiles, JLD, CSV, DataFrames
using Plots, Combinatorics, Dates, StatsPlots
using Parameters: @unpack, @with_kw
#---------------------------------------------------------------------------------------------------------

market_parameters = @with_kw (
    α_0 = 10,
    α_1 = 1,
    α_2 = 1,
    α_3 = 1,
    γ_0 = 1,
    γ_1 = 1,
    γ_2 = 1,
    γ_3 = 1,
    θ = 0.5,
    σ = 1,
    T = 50,
    S = 1000,
)


mutable struct market_data
    Q::Vector{Float64}
    W::Vector{Float64}
    R::Vector{Float64}
    Z::Vector{Float64}
    IV::Matrix{Float64}
    P::Vector{Float64}
    Y::Vector{Float64}
end



function simulation_data(parameter, s)

    @unpack α_0, α_1, α_2, α_3,γ_0 , γ_1 ,γ_2 ,γ_3, θ,σ ,T = parameter
    
    Random.seed!(s * 12345)

    Q = Float64[];
    P = Float64[];
    W = Float64[];
    R = Float64[];
    Z = Float64[];
    Y = Float64[];

    IV = zeros(T,2);
    for t = 1:T
        z_t = rand(Normal(10,1))
        w_t = rand(Normal(3,1))
        r_t = rand(Normal(0,1))
        iv_w_t = w_t + randn()
        iv_r_t = r_t + randn()
        Y_t = randn()

        ε_d = rand(Normal(0,σ))
        ε_c = rand(Normal(0,σ))

        Q_t = (α_0 + α_3 * Y_t - γ_0 - γ_2 * w_t - γ_3 * r_t + ε_d - ε_c)/((1+θ) * (α_1 + α_2 *z_t) + γ_1)

        p_t = α_0 - α_1 * Q_t  - α_2 * z_t * Q_t + α_3 * Y_t + ε_d 

        push!(Q, Q_t)
        push!(P, p_t)
        push!(W, w_t)
        push!(R, r_t)
        push!(Z, z_t)
        push!(Y, Y_t)
        IV[t,:] .= (iv_w_t, iv_r_t)
    end

    data = market_data(Q, W, R, Z, IV, P, Y)
    return data
end


function TwoSLS_estimation(parameter, data)

    Q = data.Q
    W = data.W
    R = data.R
    Z = data.Z
    P = data.P
    IV = data.IV
    Y  = data.Y

    @unpack  T = parameter

    Z_demand = hcat(ones(T), Z, IV, Y)
    Q_hat = Z_demand * inv(Z_demand' * Z_demand) * Z_demand' *  Q
    X_d = hcat(ones(T), -Q_hat, -Z .* Q_hat, Y)
    α_hat = inv(X_d' * X_d) * (X_d' * P)


    #=
    data_1st = DataFrame(Q = Q, Z = Z, IV_w = IV[:,1], IV_r = IV[:,2], Y = Y)
    iv_1st = lm(@formula(Q ~ Z + IV_w + IV_r + Y), data_1st)

    Q_hat = predict(iv_1st)

    data_2nd = DataFrame(P = P, Q_hat = -Q_hat, Q_Z = -Z .* Q_hat, Y = Y)
    iv_2nd= lm(@formula(P ~ Q_hat + Q_Z + Y), data_2nd)
    α_hat = coef(iv_2nd)
    =#
    
    # Supply side estimation 

    Z_supply = hcat(ones(T), W, R, (α_hat[2] .+ α_hat[3] .* Z), Y)
    Q_hat = Z_supply * inv(Z_supply' * Z_supply) * Z_supply' *  Q
    X_s = hcat(ones(T), Q_hat, W, R, (α_hat[2] .+ α_hat[3] .* Z) .* Q_hat)
    γ_hat = inv(X_s' * X_s) * (X_s' * P)


    θ_hat = γ_hat[end]

    return α_hat, γ_hat[1:end-1], θ_hat
end


function simulation_2sls()

    parameter = market_parameters()
    α_est = Vector{Float64}[]
    γ_est = Vector{Float64}[]
    θ_est = Float64[]

    for s = 1:1000
        data_s = simulation_data(parameter, s);
        α_est_s, γ_est_s, θ_est_s = TwoSLS_estimation(parameter, data_s)

        push!(α_est, α_est_s)
        push!(γ_est, γ_est_s)
        push!(θ_est, θ_est_s)

    end

    α_est = reduce(vcat, α_est')
    γ_est = reduce(vcat, γ_est')
    θ_est = reduce(vcat, θ_est')

    @show mean(α_est, dims = 1)
    @show mean(γ_est, dims = 1)
    @show mean(θ_est, dims = 1)
end


simulation_2sls()

Results = []

for sigma in [0.001, 0.5, 1, 2], t in [50, 100, 200, 1000]
    

    parameter = market_parameters(σ = sigma, T = t)
    α_est = Vector{Float64}[]
    γ_est = Vector{Float64}[]
    θ_est = Float64[]

    for s = 1:1000
        data_s = simulation_data(parameter, s);
        α_est_s, γ_est_s, θ_est_s = TwoSLS_estimation(parameter, data_s)

        push!(α_est, α_est_s)
        push!(γ_est, γ_est_s)
        push!(θ_est, θ_est_s)
    end

    α_est = reduce(vcat, α_est')
    γ_est = reduce(vcat, γ_est')
    θ_est = reduce(vcat, θ_est')

    result = DataFrame(
    σ = sigma,
    T = t,
    α_0 = α_est[:,1],
    α_1 = α_est[:,2],
    α_2 = α_est[:,3],
    α_3 = α_est[:,end],
    γ_0 = γ_est[:,1],
    γ_1 = γ_est[:,2],
    γ_2 = γ_est[:,3],
    γ_3 = γ_est[:,4],
    θ = θ_est)

    push!(Results, describe(result, :mean, :std, :min, :max, :median))

end













function TwoSLS_estimation(parameter, data)


    @unpack T = parameter

    Q  = data.Q
    w  = data.W
    r  = data.R
    z  = data.Z
    iv = data.IV
    p  = data.P
    Y  = data.Y

    P = Vector{Float64}[]
    Z = Matrix{Float64}[]
    X = Matrix{Float64}[]

    for t = 1:T
        Z_td = vcat(z[t], w[t], r[t], Y[t], iv[t,:])
        Z_ts = vcat(z[t], w[t], r[t], Y[t])

        Z_t = [Z_td zeros(length(Z_td));  zeros(length(Z_ts)) Z_ts]'

        X_td = vcat(1, -Q[t],-z[t].*Q[t], Y[t])
        X_ts = vcat(1, Q[t], w[t], r[t], z[t].*Q[t])

        X_t = [X_td zeros(length(X_td)); zeros(length(X_ts)) X_ts]'

        push!(P, vcat(p[t], p[t]))
        push!(X, X_t)
        push!(Z, Z_t)
    end

    Z = reduce(vcat,(Z))
    X = reduce(vcat,(X))
    P = reshape(reduce(vcat,transpose.(P))', (T * 2))

    β_hat = inv(X' * Z * inv(Z'Z) * Z' * X) * X' * Z * inv(Z'Z) * Z' * P

    α_0_hat, α_1_hat, α_2_hat  = β_hat[1], β_hat[2], β_hat[3]
    γ_0_hat, γ_2_hat, γ_3_hat = β_hat[4], β_hat[7],   β_hat[8]

    θ_hat = β_hat[6]/α_2_hat
    γ_1_hat = β_hat[5] - θ_hat * α_1_hat

    #return α_hat, γ_hat[1:end-1], θ_hat
    return β_hat, (α_0_hat, α_1_hat, α_2_hat ), (γ_0_hat, γ_1_hat, γ_2_hat, γ_3_hat), θ_hat
end


function ThreeSLS_estimation(parameter, data)

    @unpack T = parameter


    Q  = data.Q
    w  = data.W
    r  = data.R
    z  = data.Z
    iv = data.IV
    P  = data.P

    Z = []
    Y = []
    X = []

    for t = 1:T
       
        Z_ts = vcat(z[t], z[t] * Q[t], w[t], r[t], iv[t,:])
        Z_td = vcat(z[t], z[t] * Q[t], w[t], r[t], iv[t,:])

        Z_t = [Z_td zeros(length(Z_ts));  zeros(length(Z_td)) Z_ts]'

        X_td = vcat(1, -Q[t], -z[t].*Q[t])
        X_ts = vcat(1, Q[t], w[t], r[t], z[t].*Q[t])

        X_t = [X_td zeros(length(X_td)); zeros(length(X_ts)) X_ts]'

        push!(Y, vcat(P[t], P[t]))
        push!(X, X_t)
        push!(Z, Z_t)
    end


    Z = reduce(vcat,(Z))
    X = reduce(vcat,(X))
    Y = reshape(reduce(vcat,transpose.(Y))', (T * 2))

    β_hat = inv(X' * Z * inv(Z'Z) * Z' * X) * X' * Z * inv(Z'Z) * Z' * Y
    u = Y - X * β_hat
    Ω_hat = zeros(2,2)
    
    for t = 1:T
        Ω_hat .+= u[2*t-1:2*t] * u[2*t-1:2*t]'./T
    end

    W_hat = inv(Z' * ( kron(diagm(ones(T)), Ω_hat) )  * Z )

    β_hat = inv(X' * Z * W_hat * Z' * X) * X' * Z * W_hat * Z' * Y

    α_0_hat, α_1_hat, α_2_hat  = β_hat[1], β_hat[2], β_hat[3]
    γ_0_hat, γ_2_hat, γ_3_hat = β_hat[4], β_hat[7],   β_hat[8]

    θ_hat = β_hat[6]/α_2_hat
    γ_1_hat = β_hat[5] - θ_hat * α_1_hat

    return β_hat, (α_0_hat, α_1_hat, α_2_hat ), (γ_0_hat, γ_1_hat, γ_2_hat, γ_3_hat), θ_hat

end

β_hat, α_hat, γ_hat, θ_hat = TwoSLS_estimation(parameter, data)

β_hat, α_hat, γ_hat, θ_hat = ThreeSLS_estimation(parameter, data)






function Tsls_model(s)

    α = [1, 2, 3, 4]
    γ = [1, 2, 3, 4]


    X = Vector{Float64}[]
    Z = Vector{Float64}[]
    Y = Float64[]

    Random.seed!(s * 100)
    for t = 1:1000

        ε = randn()
        η = randn()

        Z_t = vcat(1, randn(3))
        X_exo = randn(2)
        X_end = γ' * Z_t + η
        X_t = vcat(1, X_exo, X_end)
        Y_t = α' * X_t + ε

        push!(Z, Z_t)
        push!(Y, Y_t)
        push!(X, X_t)
    end

    Z = reduce(vcat,(Z'))
    X = reduce(vcat,(X'))
    Y = reduce(vcat,(Y))

    X_end_hat = Z * (inv(Z' * Z) * Z' * X[:,end])

    X = hcat(X[:,1:end-1], X_end_hat)
    α_IV = inv(X' * X) * (X' * Y)

    #display(α_IV)

    #α_IV = inv(X' * Z * inv(Z' * Z) * Z' * X) * (X' * Z * inv(Z' * Z) * Z' * Y) 

    return α_IV
end

α_IV = Vector{Float64}[]
for s = 1:1000
    α_IV_s =  Tsls_model(s)

    push!(α_IV, α_IV_s)
end

α_IV = reduce(vcat, α_IV')

mean(α_IV, dims = 1)




function ols_model()

    α = [1, 2, 3, 4]

    X = Vector{Float64}[]
    Y = Float64[]

    for t = 1:100
        ε = randn()

        X_t = vcat(1, randn(3))
        Y_t = α' * X_t + ε

        push!(Y, Y_t)
        push!(X, X_t)
    end

    X = reduce(vcat,(X'))
    Y = reduce(vcat,(Y))

    α_ols = inv(X' * X) * (X' * Y)

    return α_ols
end

α_ols = Vector{Float64}[]
for s = 1:1000
    α_ols_s =  ols_model()

    push!(α_ols, α_ols_s)
end

α_ols = reduce(vcat, α_ols')

mean(α_ols, dims = 1)

#---------------------------------------------------------------------------------------------------------

## Replication of Hyde and Perloff (2015)

#---------------------------------------------------------------------------------------------------------

market_parameters_log = @with_kw (
    θ = 0.5,
    σ = 1,
    T = 50,
    α = 1,
    β = 2,
    γ = 3,
    A = 1.2,
    α_0 = 1.8,
    α_1 = 1.2,
    α_2 = -0.5,
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
    MC::Vector{Float64}
    IV::Matrix{Float64}
end

parameter  = market_parameters_log()

function simulation_data_log(parameter)


    Random.seed!(1234)

    @unpack α_1, α_2, γ_0, γ_1, γ_2, γ_3, α_0, θ, σ, T, α, β, γ, A = parameter
    
    Q = Float64[];
    P = Float64[];
    MC = Float64[];
    W  = Float64[];
    R  = Float64[];
    Z = Float64[];
    
    IV = zeros(T,2);

    for t = 1:T
        r_t = rand(Uniform(0,1))
        z_t = rand(Uniform(5,10))
        w_t = rand(Uniform(1,3))
        iv_w_t = w_t + randn()
        iv_r_t = r_t + randn()


        ε_c = rand(Normal(0,σ))
        ε_d = rand(Normal(0,σ))
        log_Q_star = (α_0 + log(1 - θ * (α_1 + α_2 * z_t)) - ( γ_0 + γ_2 * log(w_t) + γ_3 * log(r_t)) )/ (γ_1 + α_1 + α_2 * z_t) # Equilibrium total quantity
        log_Q_t = log_Q_star +  (ε_d- ε_c)/(γ_1 + α_1 + α_2 * z_t)

        mc_t = γ_0 + γ_1 *log_Q_t + γ_2 * log(w_t) + γ_3 * log(r_t)+ ε_c # Marginal Cost function

        log_p_t = α_0 - (α_1 + α_2 * z_t )* log_Q_t + ε_d # The demand function

        push!(Q, log_Q_t)
        push!(MC,mc_t)
        push!(P,log_p_t)
        push!(W,log(w_t))
        push!(R,log(r_t))
        push!(Z,z_t)
        IV[t,:] .= (iv_w_t, iv_r_t)

    end

    data = market_data_log(Q, W, R, Z, P, MC, IV)
    return data
end

data_log = simulation_data_log(parameter);

function GMM_estimation(T, Y, Z, X_s, X_d, Ω)
    
    L = size(Z, 2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "tol", 1e-15)
    set_optimizer_attribute(model, "max_iter", 1000)
    set_optimizer_attribute(model, "acceptable_tol", 1e-12)
    @variable(model, α[1:K_d])
    @variable(model, γ[1:K_s-1])
    @variable(model, θ)
    @variable(model, g[1:L])
    
    @expression(model, r_d[t = 1:T], Y[t] - sum(α[k] * X_d[t,k] for k = 1:K_d) )
    @NLexpression(model, r_s[t = 1:T], Y[t] - sum(γ[k] * X_s[t,k] for k = 1:K_s-1) + log(1 - θ*(α[2] + α[3] * X_s[t,end])))
    
    @constraint(model, [l = 1:Int64(L/2)], g[l] - sum(Z[t,l]*r_d[t] for t = 1:T) == 0)  
    @NLconstraint(model, [l = Int64(L/2)+1:L], g[l] - sum(Z[t,l] * r_s[t] for t = 1:T) ==0)

    @NLobjective(model, Min, sum( g[l] *Ω[l,k] * g[k] for l = 1:L, k = 1:L))
    optimize!(model)
    
    α_hat = value.(α)
    γ_hat = value.(γ)
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
    P  = data.P

    Y = []
    X_d = []
    X_s = []
    Z   = []

    for t = 1:T

        Z_ts = vcat(z[t], z[t] * Q[t], w[t], r[t], iv[t,:])
        Z_td = vcat(z[t], z[t] * Q[t], w[t], r[t], iv[t,:])

        Z_t = [Z_td zeros(length(Z_ts));  zeros(length(Z_td)) Z_ts]'


        X_dt = hcat(1, -Q[t], -z[t].*Q[t])
        X_st = hcat(1, Q[t], w[t], r[t], z[t])

        push!(Y, P[t])
        push!(X_d, X_dt)
        push!(X_s, X_st)
        push!(Z, Z_t)

    end

    Z = reduce(vcat,(Z))
    X_d = reduce(vcat,(X_d))
    X_s = reduce(vcat,(X_s))

    Ω = inv(Z'Z/T)


    display(Ω)
    α_hat, γ_hat, θ_hat = GMM_estimation(T, Y, Z, X_s, X_d, Ω)
    

    return α_hat, γ_hat, θ_hat
end


α_2sls, γ_2sls, θ_2sls = nonlinear_TwoSLS(parameter, data_log)


function nonlinear_ThreeSLS(parameter, data, α_2sls, γ_2sls, θ_2sls)

    @unpack T = parameter


    Q  = data.Q
    w  = data.W
    r  = data.R
    z  = data.Z
    iv = data.IV
    P  = data.P

    Y = []
    X_d = []
    X_s = []
    Z   = []

    for t = 1:T

        Z_ts = vcat(z[t], z[t] * Q[t], w[t], r[t], iv[t,:])
        Z_td = vcat(z[t], z[t] * Q[t], w[t], r[t], iv[t,:])

        Z_t = [Z_td zeros(length(Z_ts));  zeros(length(Z_td)) Z_ts]'


        X_dt = hcat(1, -Q[t], -z[t].*Q[t])
        X_st = hcat(1, Q[t], w[t], r[t], z[t])

        push!(Y, P[t])
        push!(X_d, X_dt)
        push!(X_s, X_st)
        push!(Z, Z_t)

    end

    Z = reduce(vcat,(Z))
    X_d = reduce(vcat,(X_d))
    X_s = reduce(vcat,(X_s))
    
    L = size(Z, 2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    u_d = Y .- sum(α_2sls[k] .* X_d[:,k] for k = 1:K_d)
    u_s = Y .- sum(γ_2sls[k] .* X_s[:,k] for k = 1:K_s-1) .+ log.(1 .- θ_2sls*(α_2sls[2] .+ α_2sls[3] * X_s[:,end]))


    u_hat = reduce(vcat, [[u_d[t], u_s[t]] for t = 1:T])

    Ω = zeros(L,L)

    for t = 1:T
        Ω .+= Z[2*t-1:2*t,:]' * u_hat[2*t - 1:2*t] * u_hat[2*t - 1:2*t]' * Z[2*t-1:2*t,:]
    end

    Ω  = inv(Ω/T)

    α_hat, γ_hat, θ_hat =  GMM_estimation(T, Y, Z, X_s, X_d, Ω)
    return α_hat, γ_hat, θ_hat
end

nonlinear_ThreeSLS(parameter, data_log, α_2sls, γ_2sls, θ_2sls )






