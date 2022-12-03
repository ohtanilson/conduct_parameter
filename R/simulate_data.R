library(magrittr)
# set constant ----
# one thousand replications of experiments with 50 observations each
## set variable ----
set.seed(1)
n_observation <-
  50
k_observation <-
  1000
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
  1
epsilon_c <-
  rnorm(nk, mean = 0, sd = sigma)
epsilon_d <-
  rnorm(nk, mean = 0, sd = sigma)



# generate data ----
# aggregate quantity ----
Q <- 
  (alpha0 + alpha3 * y - gamma0 - gamma2 * w - gamma3 * r + (epsilon_d - epsilon_c))/
  ((1 + theta) * (alpha1 + alpha2 * z) + gamma1)
# aggregate price ----  
P <- 
  alpha0 - (alpha1 + alpha2 * z) * Q + alpha3 * y + epsilon_d

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
plot(data$Q,
     data$P)

# estimation ----
## demand ----
demand_formula <-
  "P ~ Q + Q:z + y|y + z + iv_w + iv_r"
supply_formula <-
  paste("P ~ composite_z:Q + Q + w + r|",
        "composite_z + w + r + y")
res_demand <-
  AER::ivreg(
    formula = demand_formula,
     data = data)
res_demand <-
  data %>% 
  split(
    .$group_id_k
  ) %>% 
  purrr::map(
    ~ AER::ivreg(
      formula = demand_formula,
      data = .x)
  ) %>% 
  purrr::map(summary)# %>%
  #purrr::map_dbl("coefficients")
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

data <-
  data %>% 
  dplyr::left_join(
    demand_hat,
    by = c("group_id_k" = 
             "group_id_k")
  ) %>% 
  dplyr::mutate(
    composite_z =
      alpha1_hat + 
      alpha2_hat * z
  )
## supply ----
res_supply <-
  AER::ivreg(
    formula = supply_formula,
     data = data)
res_supply <-
  data %>% 
  split(
    .$group_id_k
  ) %>% 
  purrr::map(
    ~ AER::ivreg(
      formula = supply_formula,
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

data <-
  data %>% 
  dplyr::left_join(
    supply_hat,
    by = c("group_id_k" = 
             "group_id_k")
  ) 
# merge ----
grouped_data <-
  data %>% 
  dplyr::distinct(
    group_id_k,
    .keep_all = T
  ) %>% 
  dplyr::select(
    group_id_k,
    alpha0_hat:alpha3_hat,
    gamma0_hat:theta_hat
  )
modelsummary::datasummary_skim(grouped_data)
# save ----
saveRDS(data,
        file = "R/data.rds")
