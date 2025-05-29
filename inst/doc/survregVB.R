## ----pre,echo=FALSE,results='hide'--------------------------------------------
library(knitr)
opts_chunk$set(warning = FALSE, message = FALSE, cache = FALSE)

## ----loadLibrary--------------------------------------------------------------
library(survregVB)
library(survival)

## -----------------------------------------------------------------------------
fit <- survregVB(
  formula = Surv(time, infect) ~ trt + fev,
  data = dnase,
  alpha_0 = 501,
  omega_0 = 500,
  mu_0 = c(4.4, 0.25, 0.04),
  v_0 = 1,
  max_iteration = 10000,
  threshold = 0.0005,
  na.action = na.omit
)
print(fit)
summary(fit)

## -----------------------------------------------------------------------------
fit_frailty <- survregVB(
  formula = Surv(Time.15, delta.15) ~ x1 + x2,
  data = simulation_frailty,
  alpha_0 = 3,
  omega_0 = 2,
  mu_0 = c(0, 0, 0),
  v_0 = 0.1,
  lambda_0 = 3,
  eta_0 = 2,
  cluster = cluster,
  max_iteration = 100,
  threshold = 0.01
)
print(fit_frailty)
summary(fit_frailty)

## ----echo=FALSE---------------------------------------------------------------
sessionInfo()

