using Pkg;
 Pkg.add("LinearAlgebra")
 Pkg.add("Random")
 Pkg.add("Distributions")
 Pkg.add("CSV")
 Pkg.add("DataFrames")
 Pkg.add("Plots")
 Pkg.add("VegaLite")
 Pkg.add("Parameters")
 Pkg.add("JuMP")
 Pkg.add("Ipopt")
 Pkg.add("Optim")
 Pkg.add("RData")
#-----------------------------------------------------------------
#-----------------------------------------------------------------
# This code is based on https://github.com/magerton/julia-slurm-example

const IN_SLURM = "SLURM_JOBID" in keys(ENV)

#load packages
using Distributed
IN_SLURM && using ClusterManagers

# Adding processors
# If we are not using a cluster computing, IN_SLURM has true boolian variable. Then ClusterManagers becomes active and add the processors in the cluster. parse() changes strings into numbers. 
# If not, ClusterManagers becomes innactive and just add processors in your local machine

if IN_SLURM
    number_workers = addprocs_slurm(parse(Int, ENV["SLURM_NTASKS"]))
    print("\n")
else
    number_workers = addprocs(Sys.CPU_THREADS-1)
end

println("Pallalerization is ready \n")
println("We have ",length(workers()), " workers\n")

#----------------------------------------------------------------------------
# Run your code which needs parallelization

include("../code/03_2estimate_loglinear_loglinear_mpec.jl")

#-----------------------------------------------------------------------------
rmprocs(number_workers)
println("Processores are removed")