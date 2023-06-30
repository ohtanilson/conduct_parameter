rm(list = ls())
library(magrittr)

# set function ----
generate_data <-
  function(
    target_n_observation,
    target_k_observation,
    target_sigma,
    specification,
    demand_shifter_dummy,
    target_alpha2
  ){
    n_observation <-
      target_n_observation
    k_observation <-
      target_k_observation
    nk <-
      n_observation * k_observation
    ### exogenous variable ----
    if(specification == "linear_linear"){
      # linear specification
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
        target_alpha2#1
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
    }else{
      # log-linear specification
      w <-
        runif(nk, min = 1, max = 3)
      r <-
        #runif(nk, min = 0, max = 1)
        runif(nk, min = 1, max = 3)
      y <-
        runif(nk, min = 1, max = 3)
        #rnorm(nk, mean = 0, sd = 1)
      z <-
        runif(nk, min = 0, max = 1)
      ### instrumental variable ----
      iv_w <-
        w + rnorm(nk, mean = 0, sd = 1)
      iv_r <-
        r + rnorm(nk, mean = 0, sd = 1)
      ## set parameter ----
      theta <-
        0.2#0.5
      alpha0 <-
        20
      alpha1 <-
        1
      alpha2 <-
        target_alpha2#0.1#1
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
    }
    ## set error ----
    sigma <-
      target_sigma
    epsilon_c <-
      rnorm(nk, mean = 0, sd = sigma)
    epsilon_d <-
      rnorm(nk, mean = 0, sd = sigma)
    
    group_id_k <-
      rep(
        c(1:k_observation), 
        n_observation
      )
    # generate data ----
    if(demand_shifter_dummy == TRUE){
      if(specification == "linear_linear"){
        # aggregate quantity
        Q <- 
          (alpha0 +
             alpha3 * y - 
             gamma0 -
             gamma2 * w - 
             gamma3 * r +
             (epsilon_d -
                epsilon_c))/
          ((1 + theta) * 
             (alpha1 +
                alpha2 * z) +
             gamma1)
        # aggregate price
        P <- 
          alpha0 -
          (alpha1 +
             alpha2 * z) * Q + 
          alpha3 * y +
          epsilon_d
        data <-
          cbind(
            group_id_k,
            Q,
            P,
            w,
            r,
            y,
            z,
            iv_w,
            iv_r) %>% 
          tibble::as_tibble()
      }else{
        # aggregate quantity
        #log_Q_t = (α_0 + log(1 - θ * (α_1 + α_2 * z_t))  
        #+ α_3 * log(y_t) - ( γ_0 + γ_2 * log(w_t) + γ_3 * log(r_t)) + ε_d- ε_c)
        # / (γ_1 + α_1 + α_2 * z_t) # Equilibrium total quantity
        logQ <- 
          (alpha0 + 
             log(1 - 
                   theta *
                   (alpha1 +
                      alpha2 * z)) +
             alpha3 * log(y) - 
             (gamma0 + 
                gamma2 *log(w) +
                gamma3 *log(r)) +
             (epsilon_d -
                epsilon_c)
           )/
          (gamma1 + alpha1 + alpha2 * z)

        # aggregate price
        # log_p_t = α_0 - (α_1 + α_2 * z_t )* log_Q_t  +
        # α_3 * log_y_t + ε_d # The demand function
        logP <- 
          alpha0 - 
          (alpha1 +
             alpha2 * z) * logQ + 
              alpha3 * log(y) +
          epsilon_d
        logw <-
          log(w)
        logr <-
          log(r)
        logy <-
          log(y)
        data <-
          cbind(
            group_id_k,
            logQ,
            logP,
            logw,
            logr,
            logy,
            z,
            iv_w,
            iv_r) %>% 
          tibble::as_tibble()
      }
      
    }else{
      if(specification == "linear_linear"){
        # aggregate quantity
        Q <- 
          (alpha0 -
             0 * y -
             gamma0 -
             gamma2 * w - 
             gamma3 * r +
             (epsilon_d -
                epsilon_c))/
          ((1 + theta) *
             (alpha1 +
                alpha2 * z) +
             gamma1)
        # aggregate price
        P <- 
          alpha0 - 
          (alpha1 + alpha2 * z) * Q + 
          0 * y + 
          epsilon_d
        data <-
          cbind(
            group_id_k,
            Q,
            P,
            w,
            r,
            # drop alpha3 y
            z,
            iv_w,
            iv_r) %>% 
          tibble::as_tibble()
      }else{
        # aggregate quantity
        #log_Q_t = (α_0 + log(1 - θ * (α_1 + α_2 * z_t))  
        #+ α_3 * y_t - ( γ_0 + γ_2 * log(w_t) + γ_3 * log(r_t)) + ε_d- ε_c)
        # / (γ_1 + α_1 + α_2 * z_t) # Equilibrium total quantity
        logQ <- 
          (alpha0 + 
             log(1 - 
                   theta *
                   (alpha1 +
                      alpha2 * z)) +
             0 * log(y) - 
             (gamma0 +
                gamma2 * log(w) +
                gamma3 *log(r)) +
             (epsilon_d - 
                epsilon_c))/
          (gamma1 + 
             alpha1 + 
             alpha2 * z)
        # aggregate price
        # log_p_t = α_0 - (α_1 + α_2 * z_t )* log_Q_t  +
        # α_3 * y_t + ε_d # The demand function
        logP <-
          alpha0 - 
          (alpha1 + 
             alpha2 * z) * logQ + 
          0 * log(y) + 
          epsilon_d
        logw <-
          log(w)
        logr <-
          log(r)
        logy <-
          log(y)
        data <-
          cbind(group_id_k,
                logQ,
                logP,
                logw,
                logr,
                logy,
                z,
                iv_w,
                iv_r) %>% 
          tibble::as_tibble()
      }
      
    }
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
        target_sigma = temp_sigma,
        specification = "linear_linear",
        demand_shifter_dummy = TRUE,
        target_alpha2 = 1
      )
    filename <-
      paste(
        "output/data_linear_linear_",
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

### without demand shifter ----
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
        target_sigma = temp_sigma,
        specification = "linear_linear",
        demand_shifter_dummy = FALSE,
        target_alpha2 = 1
      )
    filename <-
      paste(
        "output/data_linear_linear_",
        "n_",
        temp_nn,
        "_sigma_",
        temp_sigma,
        "_without_demand_shifter_y",
        ".rds",
        sep = ""
      )
    cat(filename,"\n")
    saveRDS(
      data,
      file = filename)
  }
}

## log-linear demand and log-linear cost ----
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
        target_sigma = temp_sigma,
        specification = "loglinear_loglinear",
        demand_shifter_dummy = TRUE,
        target_alpha2 = 0.1#- 0.2
      )
    filename <-
      paste(
        "output/data_loglinear_loglinear_",
        "n_",
        temp_nn,
        "_sigma_",
        temp_sigma,
        ".rds",
        sep = ""
      )
    cat(filename,"\n")
    saveRDS(
      data,
      file = filename)
  }
}

### without demand shifter ----

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
        target_sigma = temp_sigma,
        specification = "loglinear_loglinear",
        demand_shifter_dummy = FALSE,
        target_alpha2 = 0.1
      )
    filename <-
      paste(
        "output/data_loglinear_loglinear_",
        "n_",
        temp_nn,
        "_sigma_",
        temp_sigma,
        "_without_demand_shifter_y",
        ".rds",
        sep = ""
      )
    cat(filename,"\n")
    saveRDS(
      data,
      file = filename)
  }
}


# Appendix ---- 
## linear with zero sigma ----
temp_sigma_list <-
  0
for(nn in 1:length(n_observation_list)){
  for(ss in 1:length(temp_sigma_list)){
    temp_nn <-
      n_observation_list[nn]
    temp_sigma <-
      temp_sigma_list[ss]
    data <-
      generate_data(
        target_n_observation = temp_nn,
        target_k_observation = 1000,
        target_sigma = temp_sigma,
        specification = "linear_linear",
        demand_shifter_dummy = TRUE,
        target_alpha2 = 1
      )
    filename <-
      paste(
        "output/data_linear_linear_",
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

## linear with zero alpha2 ----
temp_target_alpha2 <-
  0
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
        target_sigma = temp_sigma,
        specification = "linear_linear",
        demand_shifter_dummy = TRUE,
        target_alpha2 = temp_target_alpha2
      )
    filename <-
      paste(
        "output/data_linear_linear_",
        "n_",
        temp_nn,
        "_sigma_",
        temp_sigma,
        "alpha2_0",
        ".rds",
        sep = ""
      )
    cat(filename,"\n")
    saveRDS(data,
            file = filename)
  }
}
