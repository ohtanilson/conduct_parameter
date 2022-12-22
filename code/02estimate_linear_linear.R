rm(list = ls())
library(magrittr)

# set function ----
estimate_demand <-
  function(target_data,
           target_demand_formula,
           demand_shifter_dummy){
    data <-
      target_data
    res_demand <-
      data %>% 
      split(
        .$group_id_k
      ) %>% 
      purrr::map(
        ~ AER::ivreg(
          formula = target_demand_formula,
          data = .x)
      ) 
    R2_demand <-
      res_demand %>% 
      purrr::map_dbl(~summary(.)$r.squared)
    res_demand <-
      res_demand %>% 
      purrr::map(summary)
    alpha0_hat <-
      res_demand %>% 
      purrr::map_dbl(~coef(.)[1]) 
    alpha1_hat <-
      res_demand %>% 
      purrr::map_dbl(~coef(.)[2]) 
    # note the sign needs to be 
    # transformed into negative
    alpha1_hat <-
      - alpha1_hat
    if(demand_shifter_dummy == TRUE){
      # note the sign needs to be 
      # transformed into negative
      alpha2_hat <-
        res_demand %>% 
        purrr::map_dbl(~coef(.)[4]) 
      alpha2_hat <-
        - alpha2_hat
      alpha3_hat <-
        res_demand %>% 
        purrr::map_dbl(~coef(.)[3]) 
      demand_hat <-
        cbind(
          alpha0_hat,
          alpha1_hat,
          alpha2_hat,
          alpha3_hat,
          R2_demand
        ) %>% 
        tibble::as_tibble() %>% 
        dplyr::mutate(
          group_id_k =
            dplyr::row_number()
        ) 
    }else{
      alpha2_hat <-
        res_demand %>% 
        purrr::map_dbl(~coef(.)[3]) 
      # note the sign needs to be 
      # transformed into negative
      alpha2_hat <-
        - alpha2_hat
      demand_hat <-
        cbind(
          alpha0_hat,
          alpha1_hat,
          alpha2_hat,
          R2_demand
        ) %>% 
        tibble::as_tibble() %>% 
        dplyr::mutate(
          group_id_k =
            dplyr::row_number()
        ) 
    }
    data_with_demand_hat <-
      data %>% 
      dplyr::left_join(
        demand_hat,
        by = c("group_id_k" = 
                 "group_id_k")
      ) %>% 
      dplyr::mutate(
        composite_z =
          (alpha1_hat + 
               alpha2_hat * z)
      )
    return(data_with_demand_hat)
  }



estimate_supply <-
  function(target_data_with_demand_hat,
           target_supply_formula){
    data <-
      target_data_with_demand_hat
    res_supply <-
      data %>% 
      split(
        .$group_id_k
      ) %>% 
      purrr::map(
        ~ AER::ivreg(
          formula = target_supply_formula,
          data = .x)
      ) 
    R2_supply <-
      res_supply %>% 
      purrr::map_dbl(~summary(.)$r.squared)
    res_supply <-
      res_supply %>% 
      purrr::map(summary)
    gamma0_hat <-
      res_supply %>% 
      purrr::map_dbl(~coef(.)[1]) 
    gamma1_hat <-
      res_supply %>% 
      purrr::map_dbl(~coef(.)[2]) 
    gamma2_hat <-
      res_supply %>% 
      purrr::map_dbl(~coef(.)[3]) 
    gamma3_hat <-
      res_supply %>% 
      purrr::map_dbl(~coef(.)[4]) 
    theta_hat <-
      res_supply %>% 
      purrr::map_dbl(~coef(.)[5]) 
    supply_hat <-
      cbind(
        gamma0_hat,
        gamma1_hat,
        gamma2_hat,
        gamma3_hat,
        theta_hat,
        R2_supply
      ) %>% 
      tibble::as_tibble() %>% 
      dplyr::mutate(
        group_id_k =
          dplyr::row_number()
      ) 
    data_with_demand_hat_and_supply_hat <-
      data %>% 
      dplyr::left_join(
        supply_hat,
        by = c("group_id_k" = 
                 "group_id_k")
      ) 
    return(data_with_demand_hat_and_supply_hat)
  }

estimate_supply_loglinear <-
  function(target_data_with_demand_hat){
    data <-
      target_data_with_demand_hat
    g <-
      function(theta,
               data){
        constant <- 
          data[,"constant"]
        y <- 
          data[,"y"]
        logQ <-
          data[,"logQ"]
        logP <-
          data[,"logP"]
        w <-
          data[,"w"]
        r <-
          data[,"r"]
        z <-
          data[,"z"]
        iv_w <-
          data[,"iv_w"]
        iv_r <-
          data[,"iv_r"]
        composite_z <-
          data[,"composite_z"]
        linear_terms <-
          theta[1] * constant + 
          theta[2] * logQ +
          theta[3] * log(w) +
          theta[4] * log(r) 
        nonlinear_term <-
          - log(1 - theta[5] *
                  composite_z)
        u <- 
          c(logP - 
              (nonlinear_term +
              linear_terms)
            )
        return(u * y)
      }
    model_gmm_lm <-
      momentfit::momentModel(
        g = g,
        x = data %>% 
          dplyr::mutate(constant = 1) %>% 
          #dplyr::select(y, intercept, x1, x2, x3) %>% 
          as.matrix(),
        theta0 = runif(5)
        )
    fit_gmm_lm <-
      momentfit::gmmFit(
        model_gmm_lm
        )
    summary(fit_gmm_lm)
    
    res_supply <-
      data %>% 
      split(
        .$group_id_k
      ) %>% 
      purrr::map(
        ~ momentfit::gmmFit(
          formula = target_supply_formula,
          data = .x)
      ) 
  }


estimate_demand_and_supply <-
  function(target_data,
           target_demand_formula,
           target_supply_formula,
           demand_shifter_dummy){
    ## demand ----
    data_with_demand_hat <-
      estimate_demand(
        target_data = 
          target_data,
        target_demand_formula =
          target_demand_formula,
        demand_shifter_dummy =
          demand_shifter_dummy)
    ## supply ----
    data_with_demand_hat_and_supply_hat <-
      estimate_supply(
        target_data_with_demand_hat =
          data_with_demand_hat,
        target_supply_formula =
          target_supply_formula)
    ## pick up estimated parameters ----
    if(demand_shifter_dummy == T){
      parameter_hat_table <-
        data_with_demand_hat_and_supply_hat %>% 
        dplyr::distinct(
          group_id_k,
          .keep_all = T
        ) %>% 
        dplyr::select(
          group_id_k,
          alpha0_hat:alpha3_hat,
          gamma0_hat:theta_hat,
          R2_demand,
          R2_supply
        )
    }else{
      parameter_hat_table <-
        data_with_demand_hat_and_supply_hat %>% 
        dplyr::distinct(
          group_id_k,
          .keep_all = T
        ) %>% 
        dplyr::select(
          group_id_k,
          alpha0_hat:alpha2_hat, # drop alpha3
          gamma0_hat:theta_hat,
          R2_demand,
          R2_supply
        )
    }
    return(parameter_hat_table)
  }


# set constant ----
n_observation_list <-
  c(50, 100, 200, 1000)
sigma_list <-
  c(0.001, 0.5, 1.0, 2.0)

# load, estimate, and save data ----
## linear demand and linear cost ----
### with demand_shifter_y ----
for(nn in 1:length(n_observation_list)){
  for(ss in 1:length(sigma_list)){
    temp_nn <-
      n_observation_list[nn]
    temp_sigma <-
      sigma_list[ss]
    filename <-
      paste(
        "linear_linear_",
        "n_",
        temp_nn,
        "_sigma_",
        temp_sigma,
        sep = ""
      )
    cat(filename,"\n")
    # load 
    target_data <-
      readRDS(file = 
                here::here(
                  paste(
                    "output/data_",
                    filename,
                    ".rds",
                    sep = ""
                  )
                )
      )
    # assign(filename,
    #        temp_data)
    # estimate 
    linear_demand_formula <-
      "P ~ Q + Q:z + y|y + z + iv_w + iv_r"
    linear_demand_linear_supply_formula <-
      paste("P ~ composite_z:Q + Q + w + r|",
            "composite_z + w + r + y")
    parameter_hat_table <-
      estimate_demand_and_supply(
        target_data =
          target_data,
        target_demand_formula = 
          linear_demand_formula,
        target_supply_formula =
          linear_demand_linear_supply_formula,
        demand_shifter_dummy = TRUE)
    # save 
    saveRDS(parameter_hat_table,
            file = paste(
              "output/",
              "parameter_hat_table",
              filename,
              ".rds",
              sep = ""
              )
            )
  }
}
modelsummary::datasummary_skim(
  fmt = 3,
  parameter_hat_table)
### without demand_shifter_y ----
for(nn in 1:length(n_observation_list)){
  for(ss in 1:length(sigma_list)){
    temp_nn <-
      n_observation_list[nn]
    temp_sigma <-
      sigma_list[ss]
    filename <-
      paste(
        "linear_linear_",
        "n_",
        temp_nn,
        "_sigma_",
        temp_sigma,
        "_without_demand_shifter_y",
        sep = ""
      )
    cat(filename,"\n")
    # load 
    target_data <-
      readRDS(file = 
                here::here(
                  paste(
                    "output/data_",
                    filename,
                    ".rds",
                    sep = ""
                  )
                )
      )
    # assign(filename,
    #        temp_data)
    # estimate 
    linear_demand_formula <-
      "P ~ Q + Q:z|z + iv_w + iv_r"
    linear_demand_linear_supply_formula <-
      paste("P ~ composite_z:Q + Q + w + r|",
            "composite_z + w + r + iv_w + iv_r")
    parameter_hat_table <-
      estimate_demand_and_supply(
        target_data =
          target_data,
        target_demand_formula = 
          linear_demand_formula,
        target_supply_formula =
          linear_demand_linear_supply_formula,
        demand_shifter_dummy = F)
    # save 
    saveRDS(parameter_hat_table,
            file = paste(
              "output/",
              "parameter_hat_table",
              filename,
              ".rds",
              sep = ""
            )
    )
  }
}

modelsummary::datasummary_skim(
  fmt = 3,
  parameter_hat_table)

## log-linear demand and log-linear cost ----

### with demand_shifter_y ----
temp_nn <-
  n_observation_list[1]
temp_sigma <-
  sigma_list[1]
filename <-
  paste(
    "loglinear_loglinear_",
    "n_",
    temp_nn,
    "_sigma_",
    temp_sigma,
    sep = ""
  )
cat(filename,"\n")
# load 
target_data <-
  readRDS(file = 
            here::here(
              paste(
                "output/data_",
                filename,
                ".rds",
                sep = ""
              )
            )
  )
# assign(filename,
#        temp_data)
# estimate 
## demand ----
target_demand_formula =
  "logP ~ logQ + logQ:z + y|y + z + iv_w + iv_r"
data_with_demand_hat <-
  estimate_demand(
    target_data = 
      target_data,
    target_demand_formula =
      target_demand_formula,
    demand_shifter_dummy =
      TRUE)
modelsummary::datasummary_skim(
  fmt = 3,
  data_with_demand_hat %>% 
    dplyr::select(
      group_id_k,
      alpha0_hat:R2_demand
    ) %>% 
    dplyr::distinct(group_id_k,
                    .keep_all = T)
  ) 


for(nn in 1:length(n_observation_list)){
  for(ss in 1:length(sigma_list)){
    temp_nn <-
      n_observation_list[nn]
    temp_sigma <-
      sigma_list[ss]
    filename <-
      paste(
        "loglinear_loglinear_",
        "n_",
        temp_nn,
        "_sigma_",
        temp_sigma,
        sep = ""
      )
    cat(filename,"\n")
    # load
    target_data <-
      readRDS(file =
                here::here(
                  paste(
                    "output/data_",
                    filename,
                    ".rds",
                    sep = ""
                  )
                )
      )
    # assign(filename,
    #        temp_data)
    # estimate
    target_demand_formula <-
      "logP ~ logQ + logQ:z + y|y + z + iv_w + iv_r"
    data_with_demand_hat <-
      estimate_demand(
        target_data = 
          target_data,
        target_demand_formula =
          target_demand_formula,
        demand_shifter_dummy =
          TRUE)
    # save
    saveRDS(data_with_demand_hat,
            file = paste(
              "output/",
              "data_with_demand_hat_",
              filename,
              ".rds",
              sep = ""
            )
    )
  }
}



## supply ----
data <-
  data_with_demand_hat %>% 
  dplyr::mutate(
    constant =
      1
  )
fun <-
  function(theta){
  constant <- 
    data[,"constant"]
  y <- 
    data[,"y"]
  logQ <-
    data[,"logQ"]
  logP <-
    data[,"logP"]
  w <-
    data[,"w"]
  r <-
    data[,"r"]
  z <-
    data[,"z"]
  iv_w <-
    data[,"iv_w"]
  iv_r <-
    data[,"iv_r"]
  composite_z <-
    data[,"composite_z"]
  linear_terms <-
    theta[1] * constant + 
    theta[2] * logQ +
    theta[3] * log(w) +
    theta[4] * log(r) 
  nonlinear_term <-
    #- log(1 - theta[5] *
    - log(1 - 0.3 *        
            composite_z)
  u <- 
    c(logP - 
        (nonlinear_term +
           linear_terms)
    )
  return(sum(u * y))
}
# fun = function(x) 
#   x[1]*exp(-(x[1]^2 + x[2]^2)) + (x[1]^2 + x[2]^2)/20
# pracma::fminunc(x0 = c(1, 1, 1, 1), 
#                 fun)



# find numerical error point ----
constant <- 
  data[,"constant"]
y <- 
  data[,"y"]
logQ <-
  data[,"logQ"]
logP <-
  data[,"logP"]
w <-
  data[,"w"]
r <-
  data[,"r"]
z <-
  data[,"z"]
iv_w <-
  data[,"iv_w"]
iv_r <-
  data[,"iv_r"]
composite_z <-
  data[,"composite_z"]
linear_terms <-
  theta[1] * constant + 
  theta[2] * logQ +
  theta[3] * log(w) +
  theta[4] * log(r) 
nonlinear_term <-
  #- log(1 - theta[5] *
  - log(1 - 0.3 *        
          composite_z)

sum(is.na(nonlinear_term))




# appendix: test AER::ivreg, lm(2sls) ----
filename <-
  paste(
    "linear_linear_",
    "n_",
    50,
    "_sigma_",
    1,
    sep = ""
  )
target_data <-
  readRDS(file = 
            here::here(
              paste(
                "output/data_",
                filename,
                ".rds",
                sep = ""
              )
            )
  )
target_data_k <-
  target_data %>% 
  dplyr::filter(
    group_id_k == 1
  ) %>% 
  dplyr::mutate(
    Qz =
      Q * z
  )
## interaction ----
res_ivreg <-
  AER::ivreg(
    formula = "P ~ Q + Q:z + y|y + z + iv_w + iv_r",
    data = target_data_k)
res_lm_first_stage <-
  lm(
    formula = "Q ~ y + z + iv_w + iv_r",
    data = target_data_k)
target_data_k$Q_hat <-
  fitted.values(
    res_lm_first_stage
    )
res_lm_first_stage_interaction <-
  lm(
    formula = "Qz ~ y + z + iv_w + iv_r",
    data = target_data_k)
target_data_k$Q_hat_z <-
  fitted.values(
    res_lm_first_stage_interaction
  )
res_lm_second_stage <-
  lm(
    formula = "P ~ Q_hat + Q_hat_z + y",
    data = target_data_k)
## no interaction ----
res_ivreg_no_interaction <-
  AER::ivreg(
    formula = "P ~ Q + y|y + z + iv_w + iv_r",
    data = target_data_k)
res_lm_first_stage_no_interaction <-
  lm(
    formula = "Q ~ y + z + iv_w + iv_r",
    data = target_data_k)
target_data_k$Q_hat <-
  fitted.values(
    res_lm_first_stage_no_interaction
  )
res_lm_second_stage_no_interaction <-
  lm(
    formula = "P ~ Q_hat + y",
    data = target_data_k)
## compare ----
res_ivreg$coefficients
res_lm_second_stage$coefficients
res_ivreg_no_interaction$coefficients
res_lm_second_stage_no_interaction$coefficients
