using LinearAlgebra, Distributed
using Random
using Distributions
using CSV
using DataFrames
using Plots
using VegaLite
using Parameters: @with_kw, @unpack
using JuMP, Ipopt, Optim

# Set parameters for linear model
market_parameters = @with_kw (
    α_0 = 10, # Demand parameter
    α_1 = 1,
    α_2 = 1,
    α_3 = 1,
    γ_0 = 1,  # Marginal cost parameter
    γ_1 = 1,
    γ_2 = 1,
    γ_3 = 1,
    θ_0 = 0.5,  # Conduct paramter
    σ = 1,    # Standard deviation of the error term
    T = 50,   # Number of markets
    S = 100,   # Number of simulation
    start_θ = 0.0,
    start_γ = [0.0, 0.0, 0.0, 0.0]
)

mutable struct SIMULATION_SETTING
    α_0::Float64
    α_1::Float64
    α_2::Float64 
    α_3::Float64
    γ_0::Float64
    γ_1::Float64
    γ_2::Float64
    γ_3::Float64
    θ_0::Float64
    start_θ::Float64
    start_γ::Vector{Float64}
    T::Int64
    σ::Float64
    data::DataFrame
    estimation_method::Symbol
    starting_value::Symbol
    tol_level::Symbol
    simulation_index::Int64
end
