
using LinearAlgebra, Distributions
using Statistics, Random, MultivariateStats
using JuMP, Ipopt
using DelimitedFiles, JLD, CSV, DataFrames
using Plots, Combinatorics, Dates
using Parameters: @unpack, @with_kw
using Logging

#---------------------------------------------------------------------------------------------------------

market_parameters = @with_kw (
    α_1 = 1,
    α_2 = 1,
    γ_0 = 1,
    γ_1 = 1,
    γ_2 = 1,
    γ_3 = 1,
    α_3 = 0,
    α_0 = 10,
    θ = 0.5,
    σ = 2,
    T = 50,
)


mutable struct market_data
    Q::Vector{Float64}
    W::Vector{Float64}
    R::Vector{Float64}
    Z::Vector{Float64}
    IV::Matrix{Float64}
    P::Vector{Float64}
    MC::Vector{Float64}
end



function simulation_data(parameter)

    Random.seed!(1234)

    @unpack α_1, α_2 ,γ_0 , γ_1 ,γ_2 ,γ_3 , α_0, θ,σ ,T = parameter
    
    Q = Float64[];
    P = Float64[];
    MC = Float64[];
    W  = Float64[];
    R  = Float64[];
    Z = Float64[];
    
    IV = zeros(T,2);

    for t = 1:T
        r_t = rand(Normal(0,1))
        z_t = rand(Normal(10,1))
        w_t = rand(Normal(3,1))
        h_t = w_t + randn()
        k_t = r_t + randn()

        ε_c = rand(Normal(0,σ))
        ε_d = rand(Normal(0,σ))

        Q_star = (α_0 - γ_0 - γ_2 * w_t - γ_3 *r_t)/ ((1+θ) * (α_1 + α_2 *z_t) + γ_1)

        Q_t = Q_star + (ε_d - ε_c)/((1+θ) * (α_1 + α_2 *z_t) + γ_1) # The aggregate quantity, Q

        mc_t = γ_0 + γ_1 *Q_t + γ_2 * w_t + γ_3 * r_t + ε_c # Marginal cost

        p_t = α_0 - (α_1 + α_2 * z_t )* Q_t + ε_d # The demand relation

        push!(Q, Q_t)
        push!(MC,mc_t)
        push!(P,p_t)
        push!(W,w_t)
        push!(R,r_t)
        push!(Z,z_t)
        IV[t,:] .= (h_t, k_t)

    end


    data = market_data(Q, W, R, Z, IV, P, MC)
    return data
end



function TwoSLS_estimation(parameter, data)


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

        X_td = vcat(1, -Q[t],-z[t].*Q[t])
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



function Independent_TwoSLS_estimation(parameter, data)

    Q = data.Q
    W = data.W
    R = data.R
    Z = data.Z
    P = data.P
    IV = data.IV

    @unpack  α_0, α_1, α_2, γ_0, γ_1, γ_2, γ_3, θ, T = parameter

    Z_cost = hcat(W, R, Z, Z.* Q, IV)

    Q_hat = Z_cost * (inv(Z_cost' * Z_cost) * (Z_cost' * Q))
    
    X_d = hcat(ones(T), -Q_hat, -Q_hat .* Z)

    α_hat = inv(X_d' * X_d) * (X_d' * P);

    X_s = hcat(ones(T), Q_hat, W, R, (α_hat[2] .+ α_hat[3] .* Z) .* Q_hat)

    γ_hat =  inv(X_s' * X_s) * (X_s' * P); 

    θ_hat = γ_hat[end]

    return α_hat, γ_hat[1:end-1], θ_hat
end

α_hat, γ_hat, θ_hat = Independent_TwoSLS_estimation(parameter, data)


# Two additional instruments are created by adding an additional random variable drawn from N(0,1) to w and to r
# p = η + α*w + β*r +[θ * (α_1 + α_2 *z ) + γ] *Q + ε_c # Supply relationship
parameter = market_parameters()
data = simulation_data(parameter);
β_hat, α_hat, γ_hat, θ_hat = TwoSLS_estimation(parameter, data)

β_hat, α_hat, γ_hat, θ_hat = ThreeSLS_estimation(parameter, data)







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
        h_t = w_t + randn()
        k_t = r_t + randn()


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
        IV[t,:] .= (h_t, k_t)

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






