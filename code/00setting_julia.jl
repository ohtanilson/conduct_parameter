# Uncomment the followings to download packages.
# using Pkg
# Pkg.add("LinearAlgebra")
# Pkg.add("Random")
# Pkg.add("Distributions")
# Pkg.add("CSV")
# Pkg.add("DataFrames")
# Pkg.add("Plots")
# Pkg.add("VegaLite")
# Pkg.add("Parameters")
# Pkg.add("JuMP")
# Pkg.add("Ipopt")

#----------------------------------------------------------------------------------
# We assume that the reader already downloaded the above packages
using LinearAlgebra
using Random
using Distributions
using CSV
using DataFrames
using Plots
using VegaLite
using Parameters: @with_kw, @unpack
using JuMP, Ipopt, Optim

# Set parameters for log-linear model
market_parameters_log = @with_kw (
    α_0 = 10, # Demand parameter
    α_1 = 1,
    α_2 = 0.1,
    α_3 = 1,
    γ_0 = 5,  # Marginal cost parameter
    γ_1 = 1,
    γ_2 = 1,
    γ_3 = 1,
    θ = 0.2,  # Conduct paramter
    σ = 1,    # Standard deviation of the error term
    T = 50,   # Number of markets
    S = 1000, # Number of simulation
    start_θ = 0.0,
    start_γ = [0.0, 0.0, 0.0, 0.0]
)