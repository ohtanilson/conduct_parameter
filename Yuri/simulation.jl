
using LinearAlgebra, Distributions
using Statistics, Random, MultivariateStats
using JuMP, Ipopt
using DelimitedFiles, JLD, CSV, DataFrames
using Plots, Combinatorics, Dates
using Parameters: @unpack, @with_kw
using Logging


market_parameters = @with_kw (
    α_1 = 1,
    α_2 = 1,
    γ_0 = 1,
    γ_1 = 1,
    γ_2 = 1,
    γ_3 = 1,
    α_3 = 0,
    α_0 = 10,
    θ = 0/5,
    σ = 1,
    Y = 1,
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




# Two additional instruments are created by adding an additional random variable drawn from N(0,1) to w and to r


# p = η + α*w + β*r +[θ * (α_1 + α_2 *z ) + γ] *Q + ε_c # Supply relationship


parameter = market_parameters()

function simulation_data(parameter)

    @unpack α_1, α_2 ,γ_0 , γ_1 ,γ_2 ,γ_3 ,α_3,α_0, θ,σ ,Y, T = parameter
    
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
        h_t = w + randn()
        k_t = r + randn()

        ε_c = rand(Normal(0,σ))
        ε_d = rand(Normal(0,σ))

        Q_t = (α_0 + α_3 * Y - γ_0 - γ_2 * w_t - γ_3 *r_t + (ε_d - ε_c) )/ ((1+θ) * (α_1 + α_2 *z_t) + γ_1)  # Equilibrium total quantity

        mc_t = γ_0 + γ_1 *Q_t + γ_2 * w_t + γ_3 * r_t + ε_c # Marginal Cost function

        p_t = α_0 - (α_1 + α_2*z_t )* Q_t + α_3*Y + ε_d # The demand function

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

data = simulation_data(parameter);

function TSLS_estimation(parameter, data)

    Q = data.Q
    W = data.W
    R = data.R
    Z = data.Z
    IV = data.IV
    P = data.P

    @unpack T = parameter


    Z_cost = hcat(ones(T), W, R, Z)

    display(Z_cost)

    Q_hat = Z_cost * (inv(Z_cost' * Z_cost) * (Z_cost' * Q))
    
    X_d = hcat(ones(T), -Q_hat, -Q_hat .* Z)

    α_hat = inv(X_d' * X_d) * (X_d' * P);


    Q = hcat(ones(T), )



    X_s = hcat(ones(T), Q, W, R, (α_hat[1] .+ α_hat[2] .* Z) .* Q)

    γ_hat =  inv(X_s' * X_s) * (X_s' * P); 

    θ_hat = γ_hat[end]

    return α_hat, γ_hat[1:end-1], θ_hat
end

TSLS_estimation(parameter, data)