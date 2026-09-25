# Simple simulation study for the manuscript
# Ullmann T., Heinze G., Kappenberg F., Henrion M., Sauerbrei W., Collins G., Leonhardt C.-S., 
# Nold M., and Dunkler D. for TG2 of the STRATOS initiative (2026). The Problem with Univariable 
# Selection in Regression Modelling — and What To Do Instead.

library(MASS)
library(dplyr)
library(abe)
library(xtable)

set.seed(123)

alpha = 0.05
beta1 = 0.5
beta2 = 0.3
beta = c(beta1, beta2)
n = 40
rho = -0.35
r2 = 0.9

# function for determining the sigma^2 that is needed to attain a given R^2 value 
calc_sigma2 = function(beta, rho, r2) {
  D_P = mvrnorm(n = 100000, mu = c(0,0), Sigma = matrix(c(1,rho,rho,1), ncol=2))
  Xb <- as.matrix(D_P) %*% beta
  sigma2 = var(Xb)*(1-r2)/r2
  # check:
  #y <- Xb + rnorm(100000, 0, sqrt(sigma2))
  #summary(lm(y ~ Xb))$r.squared
}

sim_func <- function(nsim, n, beta, rho, sigma2) {
  results = data.frame(X1_coef_univar = numeric(), X1_p_univar = numeric(), X2_coef_univar = numeric(), X2_p_univar = numeric(),
                       X1_coef_multivar = numeric(), X1_p_multivar = numeric(), X2_p_multivar = numeric(), X2_coef_multivar = numeric())
  for (i in 1:nsim) {
    data_sim = mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(1,rho,rho,1), ncol=2))
    Xb <- as.matrix(data_sim) %*% beta
    y <- Xb + rnorm(n, 0, sqrt(sigma2))
    x1 = data_sim[,1]; x2 = data_sim[,2]
    data = data.frame(x1 = x1, x2 = x2, y = y)
    mod_univar_x1 = lm(y ~ x1, data = data)
    mod_univar_x2 = lm(y ~ x2, data = data)
    mod_multivar = lm(y ~ x1 + x2, data = data, x = T, y = T)
    mod_BE_005 = abe(mod_multivar, data = data, alpha = 0.05, tau = Inf, type.test = "F", verbose = F)
    
    results = rbind(results, data.frame(corr = cor(x1, x2), 
                                        X1_coef_univar = coef(mod_univar_x1)[2], X1_se_univar = summary(mod_univar_x1)$coefficients[2, "Std. Error"], X1_p_univar = summary(mod_univar_x1)$coefficients[2, "Pr(>|t|)"],
                                        X2_coef_univar = coef(mod_univar_x2)[2], X2_se_univar = summary(mod_univar_x2)$coefficients[2, "Std. Error"], X2_p_univar = summary(mod_univar_x2)$coefficients[2, "Pr(>|t|)"],
                                        X1_coef_multivar = coef(mod_multivar)[2], X1_se_multivar = summary(mod_multivar)$coefficients[2, "Std. Error"], X1_p_multivar = summary(mod_multivar)$coefficients[2, "Pr(>|t|)"],
                                        X2_coef_multivar = coef(mod_multivar)[3], X2_se_multivar = summary(mod_multivar)$coefficients[3, "Std. Error"], X2_p_multivar = summary(mod_multivar)$coefficients[3, "Pr(>|t|)"],
                                        BE_correct = (length(coef(mod_BE_005)) == 3)))
                                        
  }
  return(results)
}

# main scenario
set.seed(123)
sigma2 = calc_sigma2(beta = c(0.5, 0.3), rho = -0.35, r2 = 0.9)
results = sim_func(nsim = 2000, n = 40, beta = c(0.5, 0.3), rho = -0.35, sigma2)
apply(results, MARGIN = 2, FUN = mean)
sum(results$X1_p_univar < 0.05)/2000 # alpha = 0.05
sum(results$X2_p_univar < 0.05)/2000 # alpha = 0.05
sum(results$X2_p_univar < 0.2)/2000 # alpha = 0.25
sum(results$X2_p_univar < 0.5)/2000 # alpha = 0.5

mean(results$X1_coef_univar[results$X2_p_univar >= 0.05])
mean(results$X1_coef_multivar[results$X2_p_univar < 0.05])
mean(results$X2_coef_multivar[results$X2_p_univar < 0.05])

# generate table with first 10 sim. repetitions for paper
res10 = results[1:10,1:13]
rownames(res10) = 1:10
res10$X1_p_univar_alt = ifelse(res10$X1_p_univar < 0.01, "p < 0.01", paste0("p = ", round(res10$X1_p_univar,2)))
res10$X2_p_univar_alt = ifelse(res10$X2_p_univar < 0.01, "p < 0.01", paste0("p = ", round(res10$X2_p_univar,2)))

res10$univar_X1 = paste0(round(res10$X1_coef_univar, 2), ", ",res10$X1_p_univar_alt)
res10$univar_X2 = paste0(round(res10$X2_coef_univar, 2), ", ", res10$X2_p_univar_alt)

res10$X1_p_multivar_alt = ifelse(res10$X1_p_multivar < 0.01, "p < 0.01", paste0("p = ", round(res10$X1_p_multivar,2)))
res10$X2_p_multivar_alt = ifelse(res10$X2_p_multivar < 0.01, "p < 0.01", paste0("p = ", round(res10$X2_p_multivar,2)))

res10$multivar_X1 = paste0(round(res10$X1_coef_multivar, 2), ", ", res10$X1_p_multivar_alt)
res10$multivar_X2 = paste0(round(res10$X2_coef_multivar, 2), ", ", res10$X2_p_multivar_alt)

res10 = res10 %>% dplyr::select(corr, univar_X1, univar_X2, multivar_X1, multivar_X2)
res10$corr = round(res10$corr, 2)
xtable(res10)

# try different parameter values: 

# flip signs of beta_1 and rho
set.seed(123)
sigma2 = calc_sigma2(beta = c(-0.5, 0.3), rho = 0.35, r2 = 0.9)
results = sim_func(nsim = 2000, n = 40, beta = c(-0.5, 0.3), rho = 0.35, sigma2)
apply(results, MARGIN = 2, FUN = mean)
sum(results$X1_p_univar < 0.05)/2000 # alpha = 0.05
sum(results$X2_p_univar < 0.05)/2000 # alpha = 0.05

# increase sample size
set.seed(123)
sigma2 = calc_sigma2(beta = c(0.5, 0.3), rho = -0.35, r2 = 0.9)
results = sim_func(nsim = 2000, n = 80, beta = c(0.5, 0.3), rho = -0.35, sigma2)
apply(results, MARGIN = 2, FUN = mean)
sum(results$X1_p_univar < 0.05)/2000 # alpha = 0.05
sum(results$X2_p_univar < 0.05)/2000 # alpha = 0.05

# larger effect sizes
set.seed(123)
sigma2 = calc_sigma2(beta = c(0.5, 0.4), rho = -0.35, r2 = 0.9)
results = sim_func(nsim = 2000, n = 40, beta = c(0.5, 0.4), rho = -0.35, sigma2)
apply(results, MARGIN = 2, FUN = mean)
sum(results$X1_p_univar < 0.05)/2000 # alpha = 0.05
sum(results$X2_p_univar < 0.05)/2000 # alpha = 0.05
sum(results$X1_p_univar < 0.05 & results$X2_p_univar < 0.05)/2000 # alpha = 0.05

set.seed(123)
sigma2 = calc_sigma2(beta = c(0.5, 0.5), rho = -0.35, r2 = 0.9)
results = sim_func(nsim = 2000, n = 40, beta = c(0.5, 0.5), rho = -0.35, sigma2)
apply(results, MARGIN = 2, FUN = mean)
sum(results$X1_p_univar < 0.05)/2000 # alpha = 0.05
sum(results$X2_p_univar < 0.05)/2000 # alpha = 0.05
sum(results$X1_p_univar < 0.05 & results$X2_p_univar < 0.05)/2000 # alpha = 0.05

# smaller R2
set.seed(123)
sigma2 = calc_sigma2(beta = c(0.5, 0.3), rho = -0.35, r2 = 0.4)
results = sim_func(nsim = 2000, n = 40, beta = c(0.5, 0.3), rho = -0.35, sigma2)
apply(results, MARGIN = 2, FUN = mean)
sum(results$X1_p_univar < 0.05)/2000 # alpha = 0.05
sum(results$X2_p_univar < 0.05)/2000 # alpha = 0.05
sum(results$X1_p_univar < 0.05 & results$X2_p_univar < 0.05)/2000 # alpha = 0.05

# rho = -0.2
set.seed(123)
sigma2 = calc_sigma2(beta = c(0.5, 0.3), rho = -0.2, r2 = 0.9)
results = sim_func(nsim = 2000, n = 40, beta = c(0.5, 0.3), rho = -0.2, sigma2)
apply(results, MARGIN = 2, FUN = mean)
sum(results$X1_p_univar < 0.05)/2000 # alpha = 0.05
sum(results$X2_p_univar < 0.05)/2000 # alpha = 0.05

# rho = -0.5
set.seed(123)
sigma2 = calc_sigma2(beta = c(0.5, 0.3), rho = -0.5, r2 = 0.9)
results = sim_func(nsim = 2000, n = 40, beta = c(0.5, 0.3), rho = -0.5, sigma2)
apply(results, MARGIN = 2, FUN = mean)
sum(results$X1_p_univar < 0.05)/2000 # alpha = 0.05
sum(results$X2_p_univar < 0.05)/2000 # alpha = 0.05

# rho = 0.35
set.seed(123)
sigma2 = calc_sigma2(beta = c(0.5, 0.3), rho = 0.35, r2 = 0.9)
results = sim_func(nsim = 2000, n = 40, beta = c(0.5, 0.3), rho = 0.35, sigma2)
apply(results, MARGIN = 2, FUN = mean)
sum(results$X1_p_univar < 0.05)/2000 # alpha = 0.05
sum(results$X2_p_univar < 0.05)/2000

# rho = 0
set.seed(123)
sigma2 = calc_sigma2(beta = c(0.5, 0.3), rho = 0, r2 = 0.9)
results = sim_func(nsim = 2000, n = 40, beta = c(0.5, 0.3), rho = 0, sigma2)
apply(results, MARGIN = 2, FUN = mean)
sum(results$X1_p_univar < 0.05)/2000 # alpha = 0.05
sum(results$X2_p_univar < 0.05)/2000
mean(results[results$X2_p_univar >= 0.05,]$corr)

set.seed(123)
sigma2 = calc_sigma2(beta = c(1, 0.3), rho = 0, r2 = 0.9)
results = sim_func(nsim = 2000, n = 40, beta = c(1, 0.3), rho = 0, sigma2)
apply(results, MARGIN = 2, FUN = mean)
sum(results$X1_p_univar < 0.05)/2000 # alpha = 0.05
sum(results$X2_p_univar < 0.05)/2000
sum(results$X1_p_multivar < 0.05)/2000 # alpha = 0.05
sum(results$X2_p_multivar < 0.05)/2000


