rm(list = ls())
library(magrittr)

# set function ----
estimate_demand <-
  function(
    target_data,
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
        by = 
          c(
            "group_id_k" = 
              "group_id_k"
            )
        ) %>% 
      dplyr::mutate(
        composite_z =
          (alpha1_hat + 
               alpha2_hat * z)
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
    theta_hat_t_stats <-
      res_supply %>% 
      purrr::map_dbl(~coef(.)[5,"t value"]) 
    supply_hat <-
      cbind(
        gamma0_hat,
        gamma1_hat,
        gamma2_hat,
        gamma3_hat,
        theta_hat,
        theta_hat_t_stats,
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

estimate_demand_and_supply <-
  function(
    target_data,
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
          gamma0_hat:theta_hat_t_stats,
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
          gamma0_hat:theta_hat_t_stats,
          R2_demand,
          R2_supply
        )
    }
    return(parameter_hat_table)
  }


# set constant ----
n_observation_list <-
  c(
    50,
    100,
    200, 
    1000,
    2000,
    5000,
    10000
  )
sigma_list <-
  c(
    #0.001,
    #0.5, 
    1.0#, 
    #2.0
  )
alpha2_list <-
  c(
    0.1,
    0.5,
    1.0,
    5.0,
    20.0
  )
theta_list <-
  c(
    0.05,
    0.1,
    0.2,
    0.33,
    0.5,
    1.0
  )
# load, estimate, and save data ----
## linear demand and linear cost ----
### with demand_shifter_y ----
#### benchmark ----
for(nn in 1:length(n_observation_list)){
  for(ss in 1:length(sigma_list)){
    for(aa in 1:length(alpha2_list)){
      for(tt in 1:length(theta_list)){
        temp_nn <-
          n_observation_list[nn]
        temp_sigma <-
          sigma_list[ss]
        temp_theta <-
          theta_list[tt]
        temp_alpha2 <-
          alpha2_list[aa]
        filename <-
          paste(
            "data_linear_linear_",
            "n_",
            temp_nn,
            "_theta_",
            temp_theta,
            "_alpha2_",
            temp_alpha2,
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
                  "output/testing_project/",
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
            demand_shifter_dummy = TRUE
            )
        # save 
        saveRDS(
          parameter_hat_table,
          file = 
            paste(
              "output/testing_project/",
              "parameter_hat_table",
              filename,
              ".rds",
              sep = ""
            )
        )
      }
    }
  }
}
modelsummary::datasummary_skim(
  fmt = 3,
  parameter_hat_table
  )

#### 1st approximation (a second-order polynomial) ----
for(nn in 1:length(n_observation_list)){
  for(ss in 1:length(sigma_list)){
    for(aa in 1:length(alpha2_list)){
      for(tt in 1:length(theta_list)){
        temp_nn <-
          n_observation_list[nn]
        temp_sigma <-
          sigma_list[ss]
        temp_theta <-
          theta_list[tt]
        temp_alpha2 <-
          alpha2_list[aa]
        filename <-
          paste(
            "data_linear_linear_",
            "n_",
            temp_nn,
            "_theta_",
            temp_theta,
            "_alpha2_",
            temp_alpha2,
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
                  "output/testing_project/",
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
                # iv benchmark
                "composite_z + w + r + y",
                # iv polynomial 
                "+ composite_z^2 + w^2 + r^2 + y^2",
                "+ composite_z:w + composite_z:r + composite_z:y + w:r + w:y + r:y"#,
                #"+ iv_w + iv_r + iv_w^2 + iv_r^2"
                )
        parameter_hat_table <-
          estimate_demand_and_supply(
            target_data =
              target_data,
            target_demand_formula = 
              linear_demand_formula,
            target_supply_formula =
              linear_demand_linear_supply_formula,
            demand_shifter_dummy = TRUE
          )
        # save 
        saveRDS(
          parameter_hat_table,
          file = 
            paste(
              "output/testing_project/",
              "parameter_hat_table_iv_polynomial_",
              filename,
              ".rds",
              sep = ""
            )
        )
      }
    }
  }
}
modelsummary::datasummary_skim(
  fmt = 3,
  parameter_hat_table
)


# check standard error difference between julia and R ----
nn = 1
ss = 1
aa = 1
tt = 1
temp_nn <-
  n_observation_list[nn]
temp_sigma <-
  sigma_list[ss]
temp_theta <-
  theta_list[tt]
temp_alpha2 <-
  alpha2_list[aa]
filename <-
  paste(
    "data_linear_linear_",
    "n_",
    temp_nn,
    "_theta_",
    temp_theta,
    "_alpha2_",
    temp_alpha2,
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
          "output/testing_project/",
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
## demand ----
data_with_demand_hat <-
  estimate_demand(
    target_data = 
      target_data,
    target_demand_formula = 
      linear_demand_formula,
    demand_shifter_dummy = 
      TRUE)
## supply ----
# data_with_demand_hat_and_supply_hat <-
#   estimate_supply(
#     target_data_with_demand_hat =
#       data_with_demand_hat,
#     target_supply_formula =
#       linear_demand_linear_supply_formula
#     )
res_supply <-
  data_with_demand_hat %>% 
  split(
    .$group_id_k
  ) %>% 
  purrr::map(
    ~ AER::ivreg(
      formula = linear_demand_linear_supply_formula,
      data = .x)
  ) 
cat(filename,"\n")
summary(res_supply[[1]])
lmtest::coeftest(
  res_supply[[1]],
  vcov = sandwich::NeweyWest(res_supply[[1]], lag=0, prewhite=FALSE, adjust=TRUE, verbose=TRUE)
)

