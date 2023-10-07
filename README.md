# Replication files for conduct parameter papers written by Matsumura and Otani: 

This replication package contains the data and code to generate all results reported in 
(1) Resolving the Conflict on Conduct Parameter Estimation in homogeneous goods Markets between Breshnahan (1982) and Perloff and Shen (2012). Economics Letters (2023): 111193.
(2) Conduct Parameter Estimation in Homogeneous Goods Markets with Equilibrium Existence and Uniqueness Conditions: The Case of Log-linear Specification. (2023)
(3) Finite Sample Performance of a Conduct Parameter Test in Homogenous Goods Markets (2023)

For replicating results, you need to run the following R and julia scripts in /code in order.
  Data generation step
  - 01simulate_data.R # for (1) and (2)
  - 01_2_simulate_data_linear_for_test.R # for (3)
  Estimation step
  - 02estimate_linear_linear.R # for (1) 
  - 02_2_test_linear_linear.R # for (3)
  - 03_1_estimate_loglinear_loglinear.jl # for (2)
  - 03_2_estimate_loglinear_loglinear_mpec.jl # for (2)
  - 03_3_estimate_linear_linear_mpec.jl # for (2)
  - 03_4_estimate_linear_linear_test.jl # for (2)
The above processes save the results in /output and figures and tables in /figuretable.
The following HTML files correspond with .rmd scripts with the same index reports all Monte Carlo simulation results in formatted tabs and generate latex outputs for our drafts. Interested readers can check results in .html files.
  - 04_1construct_figuretable_linear_linear.Rmd # for (1) 
  - 04_1construct_figuretable_linear_linear.html
  - 04_2_construct_figuretable_loglinear_loglinear.Rmd # for (2)
  - 04_2_construct_figuretable_loglinear_loglinear.html
  - 04_3_construct_figuretable_test_linear_linear.Rmd # for (3)
  - 04_3_construct_figuretable_test_linear_linear.html
  - 04_4_construct_figuretable_test_linear_linear_iv_polynomial.Rmd # for (3)
  - 04_4_construct_figuretable_test_linear_linear_iv_polynomial.html
  - 04_5_construct_figuretable_test_linear_linear_optimal_iv.Rmd # for (3)
  - 04_5_construct_figuretable_test_linear_linear_optimal_iv.html

The remaining files are used for cluster computing and checking numerical errors by Julia.

Note that we did not push output/ folder here because the folder contains large files.
