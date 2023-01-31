module ReplicationMatsumuraOtani


# Third party packages

import LinearAlgebra
import Distributions
import CSV
import DataFrames
import Plots
import VegaLite
import Parameters
import JuMP, Ipopt


export JuMP, Ipopt, CSV


using LinearAlgebra
using Distributions
using DataFrames
using Plots
using VegaLite
using Parameters

# LinearAlgebra
export rank

# Distributions
export Uniform

# DataFrame
export DataFrame, File, dropmissing, load, sort

# Plots
export histogram, plot, vline!, hline!, histogram!, savefig, contour

# Parameters
export @unpack, @with_kw
# VegaLite
export @vlplot, save

end # module

using .ReplicationMatsumuraOtani



# Set parameters for log-linear model

market_parameters_log = @with_kw (
    α_0 = 10, # Demand parameter
    α_1 = 1,
    α_2 = 0.1,
    α_3 = 1,
    γ_0 = 1,  # Marginal cost parameter
    γ_1 = 1,
    γ_2 = 1,
    γ_3 = 1,
    θ = 0.3,  # Conduct paramter
    σ = 1,    # Standard deviation of the error term
    T = 50,   # Number of markets
    S = 1000, # Number of simulation
    start_θ = 0.0,
    start_γ = [0.0, 0.0, 0.0, 0.0]
)

parameter = market_parameters_log()

estimation_methods = [(:separate,:non_constraint, :non_constraint), (:separate,:non_constraint, :theta_constraint)];