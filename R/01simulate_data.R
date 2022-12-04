rm(list = ls())
library(magrittr)

# set function ----
generate_data <-
  function(
    target_n_observation,
    target_k_observation,
    target_sigma
  ){
    n_observation <-
      target_n_observation
    k_observation <-
      target_k_observation
    nk <-
      n_observation * k_observation
    ### exogenous variable ----
    w <-
      rnorm(nk, mean = 3, sd = 1)
    r <-
      rnorm(nk, mean = 0, sd = 1)
    y <-
      rnorm(nk, mean = 0, sd = 1)
    z <-
      rnorm(nk, mean = 10, sd = 1)
    ### instrumental variable ----
    iv_w <-
      w + rnorm(nk, mean = 0, sd = 1)
    iv_r <-
      r + rnorm(nk, mean = 0, sd = 1)
    
    ## set parameter ----
    theta <-
      0.5
    alpha0 <-
      10
    alpha1 <-
      1
    alpha2 <-
      1
    alpha3 <-
      1
    gamma0 <-
      1
    gamma1 <-
      1
    gamma2 <-
      1
    gamma3 <-
      1
    
    ## set error ----
    sigma <-
      target_sigma
    epsilon_c <-
      rnorm(nk, mean = 0, sd = sigma)
    epsilon_d <-
      rnorm(nk, mean = 0, sd = sigma)
    
    
    # generate data ----
    # aggregate quantity ----
    Q <- 
      (alpha0 + alpha3 * y - gamma0 - gamma2 * w - 
         gamma3 * r + (epsilon_d - epsilon_c))/
      ((1 + theta) * (alpha1 + alpha2 * z) + gamma1)
    # aggregate price ----  
    P <- 
      alpha0 - (alpha1 + alpha2 * z) * Q + 
      alpha3 * y + epsilon_d
    
    group_id_k <-
      rep(
        c(1:k_observation), 
        n_observation
      )
    
    data <-
      cbind(group_id_k,
            Q,
            P,
            w,
            r,
            y,
            z,
            iv_w,
            iv_r) %>% 
      tibble::as_tibble()
    return(data)
  }

generate_data_without_demand_shifter_y <-
  function(
    target_n_observation,
    target_k_observation,
    target_sigma
  ){
    n_observation <-
      target_n_observation
    k_observation <-
      target_k_observation
    nk <-
      n_observation * k_observation
    ### exogenous variable ----
    w <-
      rnorm(nk, mean = 3, sd = 1)
    r <-
      rnorm(nk, mean = 0, sd = 1)
    # y <-
    #   rnorm(nk, mean = 0, sd = 1)
    z <-
      rnorm(nk, mean = 10, sd = 1)
    ### instrumental variable ----
    iv_w <-
      w + rnorm(nk, mean = 0, sd = 1)
    iv_r <-
      r + rnorm(nk, mean = 0, sd = 1)
    
    ## set parameter ----
    theta <-
      0.5
    alpha0 <-
      10
    alpha1 <-
      1
    alpha2 <-
      1
    alpha3 <-
      1
    gamma0 <-
      1
    gamma1 <-
      1
    gamma2 <-
      1
    gamma3 <-
      1
    
    ## set error ----
    sigma <-
      target_sigma
    epsilon_c <-
      rnorm(nk, mean = 0, sd = sigma)
    epsilon_d <-
      rnorm(nk, mean = 0, sd = sigma)
    
    
    # generate data ----
    # aggregate quantity ----
    Q <- 
      (alpha0 - gamma0 - gamma2 * w - 
         gamma3 * r + (epsilon_d - epsilon_c))/
      ((1 + theta) * (alpha1 + alpha2 * z) + gamma1)
    # aggregate price ----  
    P <- 
      alpha0 - (alpha1 + alpha2 * z) * Q + 
      + epsilon_d
    
    group_id_k <-
      rep(
        c(1:k_observation), 
        n_observation
      )
    
    data <-
      cbind(group_id_k,
            Q,
            P,
            w,
            r,
            # y,
            z,
            iv_w,
            iv_r) %>% 
      tibble::as_tibble()
    return(data)
  }
# set constant ----
# one thousand replications of experiments with 50 observations each
## set list ----
set.seed(1)
n_observation_list <-
  c(50, 100, 200, 1000)
sigma_list <-
  c(0.001, 0.5, 1.0, 2.0)
# generate and save data ----
## linear demand and linear cost ----
for(nn in 1:length(n_observation_list)){
  for(ss in 1:length(sigma_list)){
    temp_nn <-
      n_observation_list[nn]
    temp_sigma <-
      sigma_list[ss]
    data_linear_linear <-
      generate_data(
        target_n_observation = temp_nn,
        target_k_observation = 1000,
        target_sigma = temp_sigma
      )
    filename <-
      paste(
        "R/output/data_linear_linear_",
        "n_",
        temp_nn,
        "_sigma_",
        temp_sigma,
        ".rds",
        sep = ""
      )
    cat(filename,"\n")
    saveRDS(data_linear_linear,
            file = filename)
  }
}

for(nn in 1:length(n_observation_list)){
  for(ss in 1:length(sigma_list)){
    temp_nn <-
      n_observation_list[nn]
    temp_sigma <-
      sigma_list[ss]
    data_linear_linear_without_demand_shifter_y <-
      generate_data_without_demand_shifter_y(
        target_n_observation = temp_nn,
        target_k_observation = 1000,
        target_sigma = temp_sigma
      )
    filename <-
      paste(
        "R/output/data_linear_linear_",
        "n_",
        temp_nn,
        "_sigma_",
        temp_sigma,
        "_without_demand_shifter_y",
        ".rds",
        sep = ""
      )
    cat(filename,"\n")
    saveRDS(data_linear_linear_without_demand_shifter_y,
            file = filename)
  }
}

## log-linear demand and linear cost ----

## log-linear demand and log-linear cost ----

## linear demand and log-linear cost ----

