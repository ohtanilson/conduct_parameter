#!/bin/bash
#SBATCH --job-name=conduct_parameter_testing    # create a short name for your job
#SBATCH --nodes=1                 # node count
#SBATCH --ntasks=20               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1         # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=5G          # memory per cpu-core (4G is default)
#SBATCH --time=24:00:00           # total run time limit (HH:MM:SS)
#SBATCH --mail-type=All           # send email when job begins and ends
#SBATCH --mail-user=ym23@rice.edu,so19@rice.edu
#SBATCH --output=conduct_parameter_testing.out

module purge
module load Julia/1.8.2

julia "../conduct_parameter/code/10_2implement_linear_GMM_optimal_instrument_on_NOTS.jl"