library(magrittr)
# set constant ----
#one thousand replications of experiments with 50 observations each
## set variable ----
set.seed(1)
n <-
  50
### exogenous variable ----
w <-
  rnorm(n, mean = 3, sd = 1)
r <-
  rnorm(n, mean = 0, sd = 1)
z <-
  rnorm(n, mean = 10, sd = 1)
### instrumental variable ----
iv_w <-
  w + rnorm(n, mean = 0, sd = 1)
iv_r <-
  r + rnorm(n, mean = 0, sd = 1)

## set parameter ----
theta <-
  0.5
alpha0 <-
  10
alpha1 <-
  1
alpha2 <-
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
  rnorm(n, mean = 0, sd = sigma)
epsilon_d <-
  rnorm(n, mean = 0, sd = sigma)



# generate data ----
# aggregate quantity ----
Q <- 
  (alpha0 - gamma0 - gamma2 * w - gamma3 * r + (epsilon_d - epsilon_c))/
  ((1 + theta) * (alpha1 + alpha2 * z) + gamma1)
# aggregate price ----  
P <- 
  alpha0 - (alpha1 + alpha2 * z) * Q + epsilon_d
  
data <-
  cbind(Q,
        P,
        w,
        r,
        z,
        iv_w,
        iv_r) %>% 
  tibble::as_tibble()
plot(data$Q,
     data$P)

res_demand <-
  AER::ivreg(
    formula = "P ~ Q + Q:z|z + iv_w + iv_r",
     data = data)

data <-
  data %>% 
  dplyr::mutate(
    alpha1_hat = res_demand$coefficients["Q"],
    alpha2_hat = res_demand$coefficients["Q:z"]
  ) %>% 
  dplyr::mutate(
    composite_z =
      alpha1_hat + 
      alpha2_hat * z
  )

res_supply <-
  AER::ivreg(
    formula = paste("P ~ composite_z:Q + Q + w + r|",
                    "composite_z + w + r + iv_w + iv_r"),
     data = data)

data <-
  data %>% 
  dplyr::mutate(
    gamma0_hat = res_supply$coefficients["(Intercept)"],
    gamma1_hat = res_supply$coefficients["Q"],
    gamma2_hat = res_supply$coefficients["w"],
    gamma3_hat = res_supply$coefficients["r"],
    theta_hat = res_supply$coefficients["composite_z:Q"]
  )

estimate_table <-
  list("Demand" = 
         res_demand,
       "Supply" = 
         res_supply
  )
modelsummary::modelsummary(
  estimate_table
)