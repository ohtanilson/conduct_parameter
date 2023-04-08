# Replication files for: Resolving the Conflict on Conduct Parameter Estimation in homogeneous goods Markets between Breshnahan (1982) and Perloff and Shen (2012). 

This replication package contains the data and the code to generate the results reported in "Matsumura, Y., & Otani, S. (2023). Resolving the Conflict on Conduct Parameter Estimation in homogeneous goods Markets between Breshnahan (1982) and Perloff and Shen (2012). arXiv preprint arXiv:2301.06665."

For replicating results, you need to run the following R scripts 
  - 01simulate_data.R
  - 01simulate_data.Rmd
  - 01simulate_data.html
  - 02estimate_linear_linear.R
  - 02estimate_linear_linear.Rmd
  - 02estimate_linear_linear.html
  - 04_1construct_figuretable.Rmd
  - 04_1construct_figuretable.html
  - 04_2construct_figuretable_bias_rmse.Rmd
  - 04_2construct_figuretable_bias_rmse.html

in code/ and save the results in output/ and figures and tables in figuretable/. 
Each HTML file corresponds with .rmd script with the same index.
Note that we did not push output/ folder here because the folder contains large files.

The rest of files for conducting the simulation of nonlinear model, but their results are not in the paper.
We are planning to write a new paper about the nonlinear model.
