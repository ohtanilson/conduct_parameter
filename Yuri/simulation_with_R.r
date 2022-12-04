library(magrittr)
library(tidyverse)
library(matlib)


R = 1000

n_observation_list <-
  c(50, 100, 200, 1000)
sigma_list <-
  c(0.001, 0.5, 1.0, 2.0)

data_variable <- 
  c("group_id_k", "Q", "P", "w", "r", "y", "z", "iv_w", "iv_r")
result_variable <- 
  c("alpha_0", "alpha_1", "alpha_2", "alpha_3", "gamma_0" , "gamma_1" , "gamma_2" , "gamma_3", "theta")


tsls_demand_estimation <- 
  function(dta, R, temp_nn, temp_sigma){
    
    result <- matrix(, ncol = 9, nrow = R)
    
    for(rr in 1:R){
      
      dta_s <- as.data.frame(dta[rr])
      names(dta_s) <- data_variable
      
      Z_d <- cbind(rep(1, temp_nn), dta_s$z, dta_s$y, dta_s$iv_w, dta_s$iv_r)
      
      Q_hat <- Z_d %*% (solve(t(Z_d) %*% Z_d) %*% t(Z_d) %*% dta_s$Q)
      
      X_d <- cbind(rep(1, temp_nn), Q_hat, dta_s$z * Q_hat, dta_s$y)
      
      alpha_hat <- solve(t(X_d) %*% X_d) %*% t(X_d) %*% dta_s$P
      
      Z_s <- cbind(rep(1, temp_nn), dta_s$z, dta_s$y, dta_s$w, dta_s$r)
      
      Q_hat <- Z_s %*% solve(t(Z_s) %*% Z_s) %*% t(Z_s) %*% dta_s$Q
      
      X_s <- cbind(rep(1, temp_nn), Q_hat, dta_s$w, dta_s$r, (alpha_hat[2] + alpha_hat[3] * dta_s$z ) * Q_hat)
  
      gamma_hat <- solve(t(X_s) %*% X_s) %*% t(X_s) %*% dta_s$P
      
      result_s <- rbind(alpha_hat, gamma_hat)
      
      result[rr, ] <- result_s
      
    }
    
    result <- as.data.frame(result)
    names(result) <- result_variable
    return( result)
  }




for(nn in 1:length(n_observation_list)){
  for(ss in 1:length(sigma_list)){
    temp_nn <-
      n_observation_list[nn]
    temp_sigma <-
      sigma_list[ss]
    
    dta <- readRDS(file = 
                     here::here(
                       paste(
                         "R/output/data_linear_linear_n_",
                         temp_nn,
                         "_sigma_",
                         temp_sigma,
                         ".rds",
                         sep = ""
                       )
                     )
    )%>% 
      split(
        .$group_id_k
      )
    

    result <- tsls_demand_estimation(dta, R, temp_nn, temp_sigma)
    print(summary(result))
  }
}


