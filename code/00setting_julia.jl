# Uncomment the followings to download packages.
#using Pkg
# Pkg.add("LinearAlgebra")
# Pkg.add("Random")
# Pkg.add("Distributions")
# Pkg.add("CSV")
# Pkg.add("DataFrames")
# Pkg.add("Plots")
#Pkg.add("VegaLite")
# Pkg.add("Parameters")
# Pkg.add("JuMP")
# Pkg.add("Ipopt")
# Pkg.add("Optim")
# Pkg.add("MathOptInterface")


number_workers = addprocs(5)

#----------------------------------------------------------------------------------
# We assume that the reader already downloaded the above packages
@everywhere using LinearAlgebra
@everywhere using Random
@everywhere using Distributions
@everywhere using CSV
@everywhere using DataFrames
@everywhere using Plots
@everywhere using MathOptInterface
@everywhere using VegaLite
@everywhere using Parameters: @with_kw, @unpack
@everywhere using JuMP, Ipopt, Optim

# Set parameters for log-linear model
@everywhere market_parameters_log = @with_kw (
    α_0 = 20, # Demand parameter
    α_1 = 1,
    α_2 = 0.1,
    α_3 = 1,
    γ_0 = 5,  # Marginal cost parameter
    γ_1 = 1,
    γ_2 = 1,    
    γ_3 = 1,
    θ_0 = 0.5,  # Conduct paramter
    σ = 1,    # Standard deviation of the error term
    T = 50,   # Number of markets
    S = 1000,   # Number of simulation
    start_θ = 0.0,
    start_γ = [0.0, 0.0, 0.0, 0.0]
)

@everywhere mutable struct SIMULATION_SETTING
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
    estimation_method::Tuple{Symbol, Symbol, Symbol,Symbol}
    starting_value::Symbol
    tol_level::Symbol
    simulation_index::Int64
end