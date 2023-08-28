rm(list = ls())
library(magrittr)
# set constant ----
sigma_list <-
  c(
    #0.001,
     0.5,
     1.0,
     2.0
  )

# load, estimate, and save data ----
# n_observation_list <-
#   c(
#     50,
#     100,
#     200,
#     1000,
#     1500
#   )
# loglinear loglinear
# for(nn in 1:length(n_observation_list)){
#   for(ss in 1:length(sigma_list)){
#     temp_nn <-
#       n_observation_list[nn]
#     temp_sigma <-
#       sigma_list[ss]
#     filename <-
#       paste(
#         "loglinear_loglinear_",
#         "n_",
#         temp_nn,
#         "_sigma_",
#         temp_sigma,
#         sep = ""
#       )
#     cat(filename,"\n")
#     # load 
#     target_data <-
#       readRDS(
#         file = 
#           here::here(
#             paste(
#               "output/data_",
#               filename,
#               ".rds",
#               sep = ""
#             )
#           )
#       )
#     assign(
#       filename,
#       target_data)
#   }
# }
n_observation_list <-
  c(
    #50,
    100, 
    200,
    1000,
    1500 #addd
  )

get_gmm_obj_value <-
  function(
    target_data,
    target_estimated_parameter
  ){
    GMM_obj_value <-
      rep(0, 1000)
    for(kk in 1:1000){
      estimated_parameters <-
        target_estimated_parameter[kk,1:9]
      target_data_k <-
        target_data %>% 
        dplyr::filter(
          group_id_k == kk
        )
      # compute residual
      residuals <-
        target_data_k %>% 
        dplyr::mutate(
          residual_demand =
            logP - 
            estimated_parameters[1] +
            (estimated_parameters[2] + 
               estimated_parameters[3] * z) * logQ -
            estimated_parameters[4] * logy,
          residual_supply =
            logP + 
            log(1 - estimated_parameters[9] * 
                  (estimated_parameters[2] + 
                     estimated_parameters[3] * z)
            ) -
            estimated_parameters[5] -
            estimated_parameters[6] * logQ -
            estimated_parameters[7] * logw -
            estimated_parameters[8] * logr
        ) %>% 
        dplyr::select(
          residual_demand,
          residual_supply
        )
      demand_iv <-
        target_data_k %>% 
        dplyr::mutate(
          constant = 1
        ) %>% 
        dplyr::select(
          constant,
          z,
          logy,
          iv_w,
          iv_r
        )
      
      supply_iv <-
        target_data_k %>% 
        dplyr::mutate(
          constant = 1
        ) %>% 
        dplyr::select(
          constant,
          z,
          logw,
          logr,
          logy
        )
      individual_moment <-
        cbind(
          as.matrix(demand_iv)[,1] * 
            residuals$residual_demand,
          as.matrix(demand_iv)[,2] * 
            residuals$residual_demand,
          as.matrix(demand_iv)[,3] * 
            residuals$residual_demand,
          as.matrix(demand_iv)[,4] * 
            residuals$residual_demand,
          as.matrix(demand_iv)[,5] * 
            residuals$residual_demand,
          as.matrix(supply_iv)[,1] * 
            residuals$residual_supply,
          as.matrix(supply_iv)[,2] * 
            residuals$residual_supply,
          as.matrix(supply_iv)[,3] * 
            residuals$residual_supply,
          as.matrix(supply_iv)[,4] * 
            residuals$residual_supply,
          as.matrix(supply_iv)[,5] * 
            residuals$residual_supply
        ) 
      moment <-
        cbind(
          mean(as.matrix(individual_moment)[,1]),
          mean(as.matrix(individual_moment)[,2]),
          mean(as.matrix(individual_moment)[,3]),
          mean(as.matrix(individual_moment)[,4]),
          mean(as.matrix(individual_moment)[,5]),
          mean(as.matrix(individual_moment)[,6]),
          mean(as.matrix(individual_moment)[,7]),
          mean(as.matrix(individual_moment)[,8]),
          mean(as.matrix(individual_moment)[,9]),
          mean(as.matrix(individual_moment)[,10])
        )
      temp_ZZ <-
        matrix(
          rep(0, 100),
          10,
          10
        )
      for(t in 1:dim(demand_iv)[1]){
        Z_t <-
          matrix(
            c(
              as.matrix(demand_iv)[t,1],
              as.matrix(demand_iv)[t,2],
              as.matrix(demand_iv)[t,3],
              as.matrix(demand_iv)[t,4],
              as.matrix(demand_iv)[t,5],
              rep(0, 5),
              rep(0, 5),
              as.matrix(supply_iv)[t,1],
              as.matrix(supply_iv)[t,2],
              as.matrix(supply_iv)[t,3],
              as.matrix(supply_iv)[t,4],
              as.matrix(supply_iv)[t,5]
            ),
            2,
            10,
            byrow = TRUE
          )
        temp_ZZ <-
          temp_ZZ +
          t(Z_t) %*% Z_t
      }
      
      weight_matrix <-
        solve(
          temp_ZZ/dim(demand_iv)[1]
        ) 
      
      GMM_obj_value[kk] <-
        moment %*%
        weight_matrix %*%
        t(moment)
      if(kk %% 100 == 0){
        cat(kk)
      }
    }
    # restore data
    return(GMM_obj_value)
  }

get_gmm_obj_value_at_true_values <-
  function(
    target_data,
    target_estimated_parameter
  ){
    GMM_obj_value <-
      rep(0, 1000)
    for(kk in 1:1000){
      true_parameter_loglinear <-
        c(
          20,#10,
          1,
          0.1,#target_alpha2
          1,
          5, # gamma0
          1,
          1,
          1,
          0.5#0.2
        )
      estimated_parameters <-
        true_parameter_loglinear
      target_data_k <-
        target_data %>% 
        dplyr::filter(
          group_id_k == kk
        )
      # compute residual
      residuals <-
        target_data_k %>% 
        dplyr::mutate(
          residual_demand =
            logP - 
            estimated_parameters[1] +
            (estimated_parameters[2] + 
               estimated_parameters[3] * z) * logQ -
            estimated_parameters[4] * logy,
          residual_supply =
            logP + 
            log(1 - estimated_parameters[9] * 
                  (estimated_parameters[2] + 
                     estimated_parameters[3] * z)
            ) -
            estimated_parameters[5] -
            estimated_parameters[6] * logQ -
            estimated_parameters[7] * logw -
            estimated_parameters[8] * logr
        ) %>% 
        dplyr::select(
          residual_demand,
          residual_supply
        )
      demand_iv <-
        target_data_k %>% 
        dplyr::mutate(
          constant = 1
        ) %>% 
        dplyr::select(
          constant,
          z,
          logy,
          iv_w,
          iv_r
        )
      
      supply_iv <-
        target_data_k %>% 
        dplyr::mutate(
          constant = 1
        ) %>% 
        dplyr::select(
          constant,
          z,
          logw,
          logr,
          logy
        )
      individual_moment <-
        cbind(
          as.matrix(demand_iv)[,1] * 
            residuals$residual_demand,
          as.matrix(demand_iv)[,2] * 
            residuals$residual_demand,
          as.matrix(demand_iv)[,3] * 
            residuals$residual_demand,
          as.matrix(demand_iv)[,4] * 
            residuals$residual_demand,
          as.matrix(demand_iv)[,5] * 
            residuals$residual_demand,
          as.matrix(supply_iv)[,1] * 
            residuals$residual_supply,
          as.matrix(supply_iv)[,2] * 
            residuals$residual_supply,
          as.matrix(supply_iv)[,3] * 
            residuals$residual_supply,
          as.matrix(supply_iv)[,4] * 
            residuals$residual_supply,
          as.matrix(supply_iv)[,5] * 
            residuals$residual_supply
        ) 
      moment <-
        cbind(
          mean(as.matrix(individual_moment)[,1]),
          mean(as.matrix(individual_moment)[,2]),
          mean(as.matrix(individual_moment)[,3]),
          mean(as.matrix(individual_moment)[,4]),
          mean(as.matrix(individual_moment)[,5]),
          mean(as.matrix(individual_moment)[,6]),
          mean(as.matrix(individual_moment)[,7]),
          mean(as.matrix(individual_moment)[,8]),
          mean(as.matrix(individual_moment)[,9]),
          mean(as.matrix(individual_moment)[,10])
        )
      temp_ZZ <-
        matrix(
          rep(0, 100),
          10,
          10
        )
      for(t in 1:dim(demand_iv)[1]){
        Z_t <-
          matrix(
            c(
              as.matrix(demand_iv)[t,1],
              as.matrix(demand_iv)[t,2],
              as.matrix(demand_iv)[t,3],
              as.matrix(demand_iv)[t,4],
              as.matrix(demand_iv)[t,5],
              rep(0, 5),
              rep(0, 5),
              as.matrix(supply_iv)[t,1],
              as.matrix(supply_iv)[t,2],
              as.matrix(supply_iv)[t,3],
              as.matrix(supply_iv)[t,4],
              as.matrix(supply_iv)[t,5]
            ),
            2,
            10,
            byrow = TRUE
          )
        temp_ZZ <-
          temp_ZZ +
          t(Z_t) %*% Z_t
      }
      
      weight_matrix <-
        solve(
          temp_ZZ/dim(demand_iv)[1]
        ) 
      
      GMM_obj_value[kk] <-
        moment %*%
        weight_matrix %*%
        t(moment)
      if(kk %% 100 == 0){
        cat(kk)
      }
    }
    # restore data
    return(GMM_obj_value)
  }


constraint_list <-
  c(
    #"_separate_non_constraint_non_constraint",
    #"_separate_theta_constraint_non_constraint",
    #"_separate_theta_constraint_slope_constraint",
    #"_separate_log_constraint_theta_constraint",
    #"_simultaneous_non_constraint_non_constraint",
    #"_simultaneous_theta_constraint_non_constraint",
    "_simultaneous_theta_constraint_slope_constraint"#,
    #"_simultaneous_log_constraint_theta_constraint",
    #"_optim_nelder_mead_theta_constraint_slope_constraint"
  )
# gmm value loglinear with constraint ----
for(nn in 1:length(n_observation_list)){
  for(ss in 1:length(sigma_list)){
    for(cc in 1:length(constraint_list)){
      ## data import ----
      temp_nn <-
        n_observation_list[nn]
      temp_sigma <-
        sigma_list[ss]
      temp_constraint <-
        constraint_list[cc]
      filename <-
        paste(
          "loglinear_loglinear_",
          "n_",
          temp_nn,
          "_sigma_",
          temp_sigma,
          sep = ""
        )
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
      assign(
        filename,
        target_data
        )
      ## estimate import ----
      filename <-
        paste(
          "loglinear_loglinear_",
          "n_",
          temp_nn,
          "_sigma_",
          temp_sigma,
          temp_constraint,
          sep = ""
        )
      cat(filename,"\n")
      temp_data <-
        readr::read_csv(
          file =
            here::here(
              paste(
                "output/",
                "parameter_hat_table_",
                filename,
                ".csv",
                sep = ""
              )
            ),
          show_col_types = FALSE
        )
      # drop T, sigma, status
      temp_data <-
        temp_data[
          ,-c(
            1, # T
            2, # sigma
            12 # status
          )]
      colnames(temp_data) <-
        c(
          "alpha0_hat",
          "alpha1_hat",
          "alpha2_hat",
          "alpha3_hat",
          "gamma0_hat",
          "gamma1_hat",
          "gamma2_hat",
          "gamma3_hat",
          "theta_hat",
          "status_indicator"
        )
      temp_data <-
        temp_data %>% 
        dplyr::rename(
          `$\\alpha_{0}$` =
            alpha0_hat,
          `$\\alpha_{1}$` =
            alpha1_hat,
          `$\\alpha_{2}$` =
            alpha2_hat,
          `$\\alpha_{3}$` =
            alpha3_hat,
          `$\\gamma_{0}$` =
            gamma0_hat,
          `$\\gamma_{1}$` =
            gamma1_hat,
          `$\\gamma_{2}$` =
            gamma2_hat,
          `$\\gamma_{3}$` =
            gamma3_hat,
          `$\\theta$` =
            theta_hat#,
          # `$R^{2}$ (demand)` =
          #   R2_demand,
          # `$R^{2}$ (supply)` =
          #   R2_supply
        ) %>% 
        dplyr::mutate(
          filename =
            paste(
              "(",
              nn,
              ") ",
              "$n=",
              temp_nn,
              "$",
              sep = ""
            ),
          locally_solved_percent =
            format(
              round(
                sum(
                  status_indicator,
                  na.rm = TRUE)/
                  dim(temp_data)[1],
                digits = 3) * 100,
              nsmall = 3)
        )
      assign(
        filename,
        temp_data
        ) 
      parameter_hat_temp_data <-
        temp_data
      
      ## compute gmm obj ----
      system.time(
        GMM_obj_value <-
          get_gmm_obj_value(
            target_data = 
              target_data,
            target_estimated_parameter =
              temp_data
          )
      )
      write.csv(
        # for minimization
        - GMM_obj_value,
        file =
          here::here(
            paste(
              "output/",
              "gmm_objective_value_",
              filename,
              ".csv",
              sep = ""
            )
          )
      )
    }
  }
}
# mpec loglinear with demand shifter
constraint_list <-
  c(#"_mpec_separate_theta_constraint_slope_constraint",
    #"_mpec_non_constraint_non_constraint",
    #"_mpec_non_constraint_theta_constraint",
    "_mpec_theta_constraint_slope_constraint"
  )
for(nn in 1:length(n_observation_list)){
  for(ss in 1:length(sigma_list)){
    for(cc in 1:length(constraint_list)){
      temp_nn <-
        n_observation_list[nn]
      temp_sigma <-
        sigma_list[ss]
      temp_constraint <-
        constraint_list[cc]
      ## data import ----
      temp_nn <-
        n_observation_list[nn]
      temp_sigma <-
        sigma_list[ss]
      temp_constraint <-
        constraint_list[cc]
      filename <-
        paste(
          "loglinear_loglinear_",
          "n_",
          temp_nn,
          "_sigma_",
          temp_sigma,
          sep = ""
        )
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
      assign(
        filename,
        target_data
      )
      ## estimate import ----
      filename <-
        paste(
          "loglinear_loglinear_",
          "n_",
          temp_nn,
          "_sigma_",
          temp_sigma,
          temp_constraint,
          sep = ""
        )
      cat(filename,"\n")
      temp_data <-
        readr::read_csv(
          file =
            here::here(
              paste(
                "output/",
                "parameter_hat_table_",
                filename,
                ".csv",
                sep = ""
              )
            ),
          show_col_types = FALSE
        )
      # drop T, sigma, status
      temp_data <-
        temp_data[
          ,-c(
            1, # T
            2, # sigma
            12 # status
          )]
      colnames(temp_data) <-
        c(
          "alpha0_hat",
          "alpha1_hat",
          "alpha2_hat",
          "alpha3_hat",
          "gamma0_hat",
          "gamma1_hat",
          "gamma2_hat",
          "gamma3_hat",
          "theta_hat",
          "status_indicator"
        )
      temp_data <-
        temp_data %>% 
        dplyr::rename(
          `$\\alpha_{0}$` =
            alpha0_hat,
          `$\\alpha_{1}$` =
            alpha1_hat,
          `$\\alpha_{2}$` =
            alpha2_hat,
          `$\\alpha_{3}$` =
            alpha3_hat,
          `$\\gamma_{0}$` =
            gamma0_hat,
          `$\\gamma_{1}$` =
            gamma1_hat,
          `$\\gamma_{2}$` =
            gamma2_hat,
          `$\\gamma_{3}$` =
            gamma3_hat,
          `$\\theta$` =
            theta_hat#,
          # `$R^{2}$ (demand)` =
          #   R2_demand,
          # `$R^{2}$ (supply)` =
          #   R2_supply
        ) %>% 
        dplyr::mutate(
          filename =
            paste(
              "(",
              nn,
              ") ",
              "$n=",
              temp_nn,
              "$",
              sep = ""
            ),
          locally_solved_percent =
            format(
              round(
                sum(
                  status_indicator,
                  na.rm = TRUE)/
                  dim(temp_data)[1],
                digits = 3) * 100,
              nsmall = 3)
        )
      assign(
        filename,
        temp_data
        ) 
      ## compute gmm obj ----
      system.time(
        GMM_obj_value <-
          get_gmm_obj_value(
            target_data = 
              target_data,
            target_estimated_parameter =
              temp_data
          )
      )
      write.csv(
        # for minimization
        - GMM_obj_value,
        file =
          here::here(
            paste(
              "output/",
              "gmm_objective_value_",
              filename,
              ".csv",
              sep = ""
            )
          )
      )
      
    }
  }
}

# gmm value loglinear with constraint at true value ----
for(nn in 1:length(n_observation_list)){
  for(ss in 1:length(sigma_list)){
    for(cc in 1:length(constraint_list)){
      ## data import ----
      temp_nn <-
        n_observation_list[nn]
      temp_sigma <-
        sigma_list[ss]
      temp_constraint <-
        constraint_list[cc]
      filename <-
        paste(
          "loglinear_loglinear_",
          "n_",
          temp_nn,
          "_sigma_",
          temp_sigma,
          sep = ""
        )
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
      assign(
        filename,
        target_data
      )
      ## estimate import ----
      filename <-
        paste(
          "loglinear_loglinear_",
          "n_",
          temp_nn,
          "_sigma_",
          temp_sigma,
          temp_constraint,
          sep = ""
        )
      cat(filename,"\n")
      temp_data <-
        readr::read_csv(
          file =
            here::here(
              paste(
                "output/",
                "parameter_hat_table_",
                filename,
                ".csv",
                sep = ""
              )
            ),
          show_col_types = FALSE
        )
      # drop T, sigma, status
      temp_data <-
        temp_data[
          ,-c(
            1, # T
            2, # sigma
            12 # status
          )]
      colnames(temp_data) <-
        c(
          "alpha0_hat",
          "alpha1_hat",
          "alpha2_hat",
          "alpha3_hat",
          "gamma0_hat",
          "gamma1_hat",
          "gamma2_hat",
          "gamma3_hat",
          "theta_hat",
          "status_indicator"
        )
      temp_data <-
        temp_data %>% 
        dplyr::rename(
          `$\\alpha_{0}$` =
            alpha0_hat,
          `$\\alpha_{1}$` =
            alpha1_hat,
          `$\\alpha_{2}$` =
            alpha2_hat,
          `$\\alpha_{3}$` =
            alpha3_hat,
          `$\\gamma_{0}$` =
            gamma0_hat,
          `$\\gamma_{1}$` =
            gamma1_hat,
          `$\\gamma_{2}$` =
            gamma2_hat,
          `$\\gamma_{3}$` =
            gamma3_hat,
          `$\\theta$` =
            theta_hat#,
          # `$R^{2}$ (demand)` =
          #   R2_demand,
          # `$R^{2}$ (supply)` =
          #   R2_supply
        ) %>% 
        dplyr::mutate(
          filename =
            paste(
              "(",
              nn,
              ") ",
              "$n=",
              temp_nn,
              "$",
              sep = ""
            ),
          locally_solved_percent =
            format(
              round(
                sum(
                  status_indicator,
                  na.rm = TRUE)/
                  dim(temp_data)[1],
                digits = 3) * 100,
              nsmall = 3)
        )
      assign(
        filename,
        temp_data
      ) 
      parameter_hat_temp_data <-
        temp_data
      
      ## compute gmm obj at true values----
      system.time(
        GMM_obj_value <-
          get_gmm_obj_value_at_true_values(
            target_data = 
              target_data,
            target_estimated_parameter =
              temp_data
          )
      )
      write.csv(
        # for minimization
        - GMM_obj_value,
        file =
          here::here(
            paste(
              "output/",
              "gmm_objective_value_at_true_values_",
              filename,
              ".csv",
              sep = ""
            )
          )
      )
    }
  }
}
