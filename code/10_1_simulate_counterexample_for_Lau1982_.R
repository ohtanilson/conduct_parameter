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

    # exogenous demand shifter
    logX_d1 <- log(runif(nk, min = 1, max = 3))
    logX_d2 <- log(runif(nk, min = 1, max = 3))

    # exogenous cost shifter
    logX_s1 <- log(runif(nk, min = 1, max = 3))
    logX_s2 <- log(runif(nk, min = 1, max = 3))
  

    # instrumental variables for cost shifter
    iv_s1 <- logX_s1 + rnorm(nk, mean = 0, sd = 1)
    iv_s2 <- logX_s2 + rnorm(nk, mean = 0, sd = 1)

    
    # ---- set parameter ----

    # conduct parameter
    theta <- 0.5

    # for demand
    alpha0 <- 5
    alpha1 <- 1
    alpha2 <- 1

    # for supply
    beta0 <- 1
    beta1 <- 1
    beta2 <- 1

    # ---- set error ----
    sigma <- target_sigma
    epsilon_d <- rnorm(nk, mean = 0, sd = target_sigma)
    epsilon_s <- rnorm(nk, mean = 0, sd = target_sigma)

    group_id_k <- rep(1:k_observation, n_observation)

    #  ---- generate data ----
    # aggregate quantity
    logQ <- 
      (beta1 * logX_s1 +
          beta2 * logX_s2 -
          log(1 + theta * alpha0) +
          epsilon_s -
          alpha1 * logX_d1 -
          alpha2 * logX_d2 -
          epsilon_d)/
      (alpha0 - beta0)

    # aggregate price
    logP <- alpha0 * logQ +
                alpha1 * logX_d1 +
                alpha2 * logX_d2 +
                epsilon_d

    data <- cbind(
      group_id_k,
      logQ,
      logP,
      logX_d1,
      logX_d2,
      logX_s1,
      logX_s2,
      iv_s1,
      iv_s2
    ) |>
      as.data.frame()
    
    return(data)
  }

# set constant ----
# one thousand replications of experiments with 50 observations each
## set list ----
set.seed(1)
n_observation_list <-
  c(
    50,
    100,
    200, 
    1000
    )
sigma_list <-
  c(
    0.001,
    0.5, 
    1.0, 
    2.0
    )

# generate and save data ----
## linear demand and linear cost ----
### with demand shifter ----
for(nn in 1:length(n_observation_list)){
  for(ss in 1:length(sigma_list)){
    temp_nn <-
      n_observation_list[nn]
    temp_sigma <-
      sigma_list[ss]
    data <-
      generate_data(
        target_n_observation = temp_nn,
        target_k_observation = 1000,
        target_sigma = temp_sigma
      )
    filename <-
      paste(
        "output/data_counterexample_for_Lau1982_",
        "n_",
        temp_nn,
        "_sigma_",
        temp_sigma,
        ".rds",
        sep = ""
      )
    cat(filename,"\n")
    saveRDS(data,
            file = filename)
  }
}

