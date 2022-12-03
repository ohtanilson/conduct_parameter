rm(list = ls())
library(magrittr)

# set function ----
estimate_demand <-
  function(target_data,
           target_demand_formula){
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
      ) %>% 
      purrr::map(summary)
    alpha0_hat <-
      res_demand %>% 
      purrr::map_dbl(~coef(.)[1]) 
    alpha1_hat <-
      res_demand %>% 
      purrr::map_dbl(~coef(.)[2]) 
    alpha3_hat <-
      res_demand %>% 
      purrr::map_dbl(~coef(.)[3]) 
    alpha2_hat <-
      res_demand %>% 
      purrr::map_dbl(~coef(.)[4]) 
    demand_hat <-
      cbind(
        alpha0_hat,
        alpha1_hat,
        alpha2_hat,
        alpha3_hat
      ) %>% 
      tibble::as_tibble() %>% 
      dplyr::mutate(
        group_id_k =
          dplyr::row_number()
      ) 
    data_with_demand_hat <-
      data %>% 
      dplyr::left_join(
        demand_hat,
        by = c("group_id_k" = 
                 "group_id_k")
      ) %>% 
      dplyr::mutate(
        composite_z =
          # note the sign needs to be 
          # transformed into negative
          - (alpha1_hat + 
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
      ) %>% 
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
        theta_hat
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

estimate_demand_and_supply <-
  function(target_data,
           target_demand_formula,
           target_supply_formula){
    ## demand ----
    data_with_demand_hat <-
      estimate_demand(
        target_data = 
          target_data,
        target_demand_formula =
          target_demand_formula)
    ## supply ----
    data_with_demand_hat_and_supply_hat <-
      estimate_supply(
        target_data_with_demand_hat =
          data_with_demand_hat,
        target_supply_formula =
          target_supply_formula)
    ## pick up estimated parameters ----
    parameter_hat_table <-
      data_with_demand_hat_and_supply_hat %>% 
      dplyr::distinct(
        group_id_k,
        .keep_all = T
      ) %>% 
      dplyr::select(
        group_id_k,
        alpha0_hat:alpha3_hat,
        gamma0_hat:theta_hat
      )
    return(parameter_hat_table)
  }

# set constant ----
n_observation_list <-
  c(50, 100, 200, 1000)
sigma_list <-
  c(0.001, 0.5, 1.0, 2.0)

# load, estimate, and save data ----
## linear demand and linear cost ----
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
                    "R/output/data_",
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
          linear_demand_linear_supply_formula)
    # save 
    saveRDS(parameter_hat_table,
            file = paste(
              "R/output/",
              "parameter_hat_table",
              filename,
              ".rds",
              sep = ""
              )
            )
  }
}

## log-linear demand and linear cost ----

## log-linear demand and log-linear cost ----

## linear demand and log-linear cost ----


