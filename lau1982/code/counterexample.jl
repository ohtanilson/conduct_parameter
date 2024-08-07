
# Assume that the reader already downloaded the above packages
using LinearAlgebra
using Random
using Distributions
using CSV
using DataFrames
using Plots
using StatsPlots


@kwdef mutable struct simulation_setting
    S :: Int64 = 1000
    T :: Int64 = 1000
    N :: Int64 = 3
    K_d :: Int64 = 3
    K_s :: Int64 = 2
    α :: Vector{Float64} = [-1, 0.5, 0.2, -0.1] # first element is for Q
    β :: Vector{Float64} = [1, -0.2, 0.1]       # first element is for Q
    θ :: Float64 = 1/N
    σ :: Float64 = 0.1
end

function simulaiton_data(setting::simulation_setting)

    (; S, T, K_d, K_s, α, β, θ, σ) = setting

    data = []
    for s = 1:S
        Random.seed!(s)

        for t = 1:T

            X_d = vcat(exp(1), rand(LogNormal(0,0.5), K_d-1))
            X_s = rand(LogNormal(0,0.5), K_s)
            ε_d = rand(Normal(0, σ))
            ε_s = rand(Normal(0, σ))
    
            logQ = (sum(β[k+1] .* log.(X_s[k]) for k = 1:K_s) .+ ε_s .- sum(α[k+1] .* log.(X_d[k]) for k = 1:K_d) .- ε_d .- log(1 + θ * α[1]) )./(α[1] - β[1])

            logP = α[1] .* logQ .+ sum(α[k+1] .* log.(X_d[k]) for k = 1:K_d) .+ ε_d
            
            P = exp.(logP)
            Q = exp.(logQ)

            data_s = DataFrame("P" => P, "Q" => Q)

            for k = 1:K_d
                data_s = hcat(data_s, DataFrame("X_d$k" => X_d[k]))
            end

            for k = 1:K_s
                data_s = hcat(data_s, DataFrame("X_s$k" => X_s[k]))
            end
            
            data_s.simulation_id .= s
            data_s.market_id .= t
    
            push!(data, data_s)
        end
    end
    data = vcat(data...)

    return data
end

setting = simulation_setting()

data = simulaiton_data(setting)

##

α_simulation = []
β_simulation = []
θ_simulation = []

for s = 1:setting.S
    data_s = data[data.simulation_id .== s, :]

    P = data_s.P
    Q = data_s.Q

    X_d = hcat([data_s[:,"X_d$k"] for k = 1:setting.K_d]...)
    X_s = hcat([data_s[:,"X_s$k"] for k = 1:setting.K_s]...)

    X = hcat(X_d, X_s)

    # demand estiamtion, first-stage
    Q_hat = exp.(log.(X) * inv(log.(X)' * log.(X)) * (log.(X)' * log.(Q)))

    α_hat = inv(log.(hcat(Q_hat, X_d))' * log.(hcat(Q_hat, X_d))) * (log.(hcat(Q_hat, X_d))' * log.(P))

    Z = hcat(exp.(ones(setting.T)), Q_hat, X_s)

    γ_hat = inv(log.(Z)' * log.(Z)) * (log.(Z)' * log.(P))
    
    β_hat = γ_hat[2:end]

    θ_hat = (exp(-γ_hat[1]) - 1)/α_hat[1]

    push!(α_simulation, α_hat')
    push!(β_simulation, β_hat')
    push!(θ_simulation, θ_hat)

end


α_simulation = vcat(α_simulation...)
β_simulation = vcat(β_simulation...)
θ_simulation = vcat(θ_simulation...)


Bias_α = mean(abs.(α_simulation .- setting.α'), dims = 1)
Bias_β = mean(abs.(β_simulation .- setting.β'), dims = 1)
Bias_θ = mean(abs.(θ_simulation .- setting.θ'), dims = 1)

RMSE_α = sqrt.(mean((α_simulation .- setting.α').^2, dims = 1))
RMSE_β = sqrt.(mean((β_simulation .- setting.β').^2, dims = 1))
RMSE_θ = sqrt.(mean((θ_simulation .- setting.θ').^2, dims = 1))