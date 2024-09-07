rm(list = ls())
library(magrittr)

# set function ----
estimate_demand <-
  function(
    target_data,
    target_demand_formula
    ){

    data <- target_data

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

    alpha2_hat <-
        res_demand %>% 
        purrr::map_dbl(~coef(.)[3]) 

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

    data_with_demand_hat <-
      data %>% 
      dplyr::left_join(
        demand_hat,
        by = 
          c(
            "group_id_k" = 
              "group_id_k"
            )
        )

    return(data_with_demand_hat)
  }



estimate_supply <-
  function(
    target_data_with_demand_hat,
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

    intercept_hat <-
        res_supply %>% 
        purrr::map_dbl(~coef(.)[1]) 

    beta0_hat <-
        res_supply %>% 
        purrr::map_dbl(~coef(.)[2]) 
    beta1_hat <-
        res_supply %>% 
        purrr::map_dbl(~coef(.)[3]) 
    beta2_hat <-
        res_supply %>% 
        purrr::map_dbl(~coef(.)[4]) 

    supply_hat <-
      cbind(
        intercept_hat,
        beta0_hat,
        beta1_hat,
        beta2_hat,
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
      )  |>
      dplyr::mutate(
        # we compute theta_hat from:
        #   -intercept_hat = log(1 + theta_hat * alpha0_hat)
        theta_hat = (exp(-intercept_hat) - 1) / alpha0_hat
      )

    return(data_with_demand_hat_and_supply_hat)
  }

estimate_demand_and_supply <-
  function(
    target_data,
    target_demand_formula,
    target_supply_formula
    ){
    ## demand ----
    data_with_demand_hat <-
      estimate_demand(
        target_data = target_data,
        target_demand_formula = target_demand_formula
        )
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
            alpha0_hat:alpha2_hat,
            beta0_hat:beta2_hat,
            theta_hat,
            R2_demand,
            R2_supply
        )   
  
    return(parameter_hat_table)
  }


# set constant ----
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
        "counterexample_for_Lau1982_",
        "n_",
        temp_nn,
        "_sigma_",
        temp_sigma,
        sep = ""
      )
    cat(filename,"\n")
    # load 
    target_data <-
      readRDS(
        file = 
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
    demand_formula <-
      "logP ~ -1 + logQ + logX_d1 + logX_d2 | logX_d1 + logX_d2 + iv_s1 + iv_s2"
    supply_formula <-
        "logP ~ logQ + logX_s1 + logX_s2 | logX_s1 + logX_s2 + logX_d1 + logX_d2"

    parameter_hat_table <-
      estimate_demand_and_supply(
        target_data = target_data,
        target_demand_formula = demand_formula,
        target_supply_formula = supply_formula
        )
    # save 
    saveRDS(
      parameter_hat_table,
      file = 
        paste(
          "output/",
          "parameter_hat_table_",
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

colnames(parameter_hat_table)
