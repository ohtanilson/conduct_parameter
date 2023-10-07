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
  paste("P ~ -1 + one + composite_z:Q + Q + w + r|",
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
res_demand <-
  target_data %>% 
  split(
    .$group_id_k
  ) %>% 
  purrr::map(
    ~ AER::ivreg(
      formula = linear_demand_formula,
      data = .x)
  ) 
summary(res_demand[[100]])
vcov(res_demand[[100]])
sum(res_demand[[100]]$residuals^2)/(50-4)
Xd <-
  as.matrix(
    cbind(rep(1, 50),
          data_with_demand_hat %>% 
            dplyr::filter(group_id_k == 100) %>%
            dplyr::mutate(Qz = z*Q) %>% 
            dplyr::select(Q,Qz,y)
    )
  )
P <-
  as.matrix(
    data_with_demand_hat %>% 
      dplyr::filter(group_id_k == 100) %>%
      dplyr::select(P)
  )
Q <-
  as.matrix(
    data_with_demand_hat %>% 
      dplyr::filter(group_id_k == 100) %>%
      dplyr::select(Q)
  )
Qz <-
  as.matrix(
    data_with_demand_hat %>% 
      dplyr::filter(group_id_k == 100) %>%
      dplyr::mutate(Qz = z*Q) %>% 
      dplyr::select(Qz)
  )
Zd <-
  as.matrix(
    cbind(rep(1, 50),
          data_with_demand_hat %>% 
            dplyr::filter(group_id_k == 100) %>%
            dplyr::select(y, z, iv_w, iv_r)
    )
  )
Qhat <-
  Zd %*% 
  solve(crossprod(Zd, Zd)) %*%
  crossprod(Zd, Q)
Qzhat <-
  Zd %*% 
  solve(crossprod(Zd, Zd)) %*%
  crossprod(Zd, Qz)
Xdhat <-
  as.matrix(
    cbind(rep(1, 50),
          Qhat,
          Qzhat,
          data_with_demand_hat %>% 
            dplyr::filter(group_id_k == 100) %>%
            dplyr::select(y)
    )
  )
alphahat <-
  solve(t(Xdhat) %*% Xdhat) %*% (t(Xdhat) %*% P)
residual_d <-
  P - Xdhat %*% alphahat
P - Xd %*% alphahat ==
  res_demand[[100]]$residuals
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
sigmahat_supply <-
  res_supply %>% 
  purrr::map_dbl(~summary(.)$sigma)
group_id_k <-
  c(1:length(unique(data_with_demand_hat$group_id_k)))
sigmahat_supply <-
  cbind(
    group_id_k,
    sigmahat_supply
  ) %>% 
  tibble::as_tibble()
# recompute using optimal iv
linear_demand_linear_supply_formula_optimal <-
  paste("P ~ -1 + one + composite_z:fitted_values_of_quantity_on_z + fitted_values_of_quantity_on_z + w + r|",
        "composite_z + fitted_values_of_quantity_on_z + w + r + y")
data_with_sigmahat <-
  data_with_demand_hat %>% 
  dplyr::left_join(
    sigmahat_supply,
    by = c("group_id_k" = "group_id_k")
  ) %>% 
  # construct optimal iv
  dplyr::mutate(
    one = one,
    Q_optimal = Q/(sigmahat_supply)^2,
    w_optimal = w/(sigmahat_supply)^2,
    r_optimal = r/(sigmahat_supply)^2,
    y_optimal = y/(sigmahat_supply)^2
  )
res_supply_optimal <-
  data_with_sigmahat %>% 
  split(
    .$group_id_k
  ) %>% 
  purrr::map(
    ~ AER::ivreg(
      formula = linear_demand_linear_supply_formula_optimal,
      data = .x)
  ) 

# res_supply <-
#   res_supply %>% 
#   purrr::map(summary)
# gamma0_hat <-
#   res_supply %>% 
#   purrr::map_dbl(~coef(.)[1]) 
# gamma1_hat <-
#   res_supply %>% 
#   purrr::map_dbl(~coef(.)[2]) 
# gamma2_hat <-
#   res_supply %>% 
#   purrr::map_dbl(~coef(.)[3]) 
# gamma3_hat <-
#   res_supply %>% 
#   purrr::map_dbl(~coef(.)[4]) 
# theta_hat <-
#   res_supply %>% 
#   purrr::map_dbl(~coef(.)[5]) 
# theta_hat_t_stats <-
#   res_supply_optimal %>% 
#   purrr::map_dbl(~coef(.)[5,"t value"]) 
summary(res_supply_optimal[[99]])
summary(res_supply[[99]])
## bread ----
Xs <-
  as.matrix(
    cbind(rep(1, 50),
          data_with_demand_hat %>% 
            dplyr::filter(group_id_k == 100) %>%
            dplyr::mutate(composite_zQ = composite_z*Q) %>% 
            dplyr::select(Q,w,r,composite_zQ)
    )
  )
Zs <-
  as.matrix(
    cbind(rep(1, 50),
          data_with_demand_hat %>% 
            dplyr::filter(group_id_k == 100) %>%
            dplyr::select(w, r, y, composite_z)
    )
  )
xhat <-
  Zs %*% 
  solve(crossprod(Zs, Zs)) %*%
  crossprod(Zs, Xs)
sandwich::bread(res_supply[[100]]) ==
  solve(crossprod(cbind(xhat), cbind(xhat))) *
  50
## meat ----
sandwich::meat(res_supply[[100]]) ==
  crossprod(res_supply[[100]]$residuals * cbind(xhat))/50
## vcov
sandwich::sandwich(res_supply[[100]]) ==
  sandwich::bread(res_supply[[100]]) %*% 
  sandwich::meat(res_supply[[100]]) %*%
  sandwich::bread(res_supply[[100]])/50 
stats::vcov(res_supply[[100]]) ==
  sandwich::vcovHC(res_supply[[100]],
                   type = "const")
stats::vcov(res_supply[[100]]) ==
  solve(crossprod(cbind(xhat), cbind(xhat))) *
  sum(res_supply[[100]]$residuals^2)/(50-5)
## se ----
summary(res_supply[[100]])
sqrt(diag(
  solve(crossprod(cbind(xhat), cbind(xhat))) *
    sum(res_supply[[100]]$residuals^2)/(50-5)
)
)


vcov <-
  solve(t(Zs) %*% Xs) %*% t(Zs) %*% 
  Zs %*% solve(t(Zs) %*% Xs) * 
  (sum(res_supply[[100]]$residuals^2)/(50 - 5))

# vcov <-
#   sum(res_supply[[100]]$residuals^2)/(50 - 5) *invXX
sqrt(diag(vcov))
vcov(res_supply[[99]])
lmtest::coeftest(
  res_supply[[1]],
  vcov = sandwich::NeweyWest(res_supply[[1]], lag=0, prewhite=FALSE, adjust=TRUE, verbose=TRUE)
)

# use example to check bread and meat ----
x <- sin(1:10)
y <- rnorm(10)
z <- rnorm(10)
fm <- AER::ivreg(y ~ x|z)
## bread ----
#             (Intercept)         x
# (Intercept)    1.616317 -4.367364
# x             -4.367364 30.948133
xhat <-
  cbind(1, z) %*% 
  solve(crossprod(cbind(1, z), cbind(1, z))) %*%
  crossprod(cbind(1, z), cbind(x))
sandwich::bread(fm) ==
  solve(crossprod(cbind(1, xhat), cbind(1, xhat))) *
  10
## meat ----
#             (Intercept)          x
# (Intercept)   3.2656098 0.11709824
# x             0.1170982 0.06016913
sandwich::meat(fm) ==
  crossprod(fm$residuals * cbind(1, xhat))/10
psi <-
  sandwich::estfun(fm)
crossprod(as.matrix(psi))/10
## vcov ----
#             (Intercept)         x
# (Intercept)   0.8025801 -2.309365
# x            -2.3093655  8.826262
sandwich::sandwich(fm) ==
  sandwich::bread(fm) %*% 
  sandwich::meat(fm) %*%
  sandwich::bread(fm)/10
stats::vcov(fm) ==
  sandwich::vcovHC(fm, type = "const")
stats::vcov(fm) ==
  solve(crossprod(cbind(1, xhat), cbind(1, xhat))) *
  sum(fm$residuals^2)/(10-2)
sandwich::sandwich(fm) ==
  sandwich::vcovHC(fm, type = "HC")
## se ----
summary(fm)
sqrt(diag(sandwich::vcovHC(fm, type = "const")))

# check se prcedure ----
df <- data.frame(ones = c(1, 1, 1, 1, 1, 1, 1),
                 rating=c(67, 75, 79, 85, 90, 96, 97),
                 points=c(8, 12, 16, 15, 22, 28, 24),
                 assists=c(4, 6, 6, 5, 3, 8, 7),
                 rebounds=c(1, 4, 3, 3, 2, 6, 7))

#fit multiple linear regression model
# some data (taken from Roland's example)
x = c(1, 2, 3, 4)
y = c(2.1, 3.9, 6.3, 7.8)
xx = cbind(c(1,1,1,1),x)
# fitting a linear model
fit = lm(y ~ x)

# get vector of all standard errors of the coefficients
summary(fit)
sqrt(diag(vcov(fit)))
vcov(fit)
sum(fit$residual^2)/(length(fit$residual)-2)*solve(t(xx) %*% xx)
