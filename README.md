# survregVB

## Overview

`survregVB` is an R package that provides Bayesian inference for log-logistic accelerated failure time (AFT) models used in survival analysis as a faster alternative to Markov chain Monte Carlo (MCMC) methods. The details of the Variational Bayes algorithms with and without shared frailty can be found in [Xian et al., (2024a)](https://doi.org/10.1007/s11222-023-10365-6) and [Xian et al., (2024b)](https://doi.org/10.48550/ARXIV.2408.00177) respectively.

## Installation

To install `survregVB`, use the following command:

``` r
remotes::install_github("https://github.com/chengqianxian/survregVB")
```

## Usage

### Loading the Package

``` r
library(survregVB)
library(survival) 
```

### Fitting a Basic Model

``` r
# Example using dataset included in the package
data(dnase)

# Fit a survival model
fit <- survregVB(formula = Surv(time, infect) ~ trt + fev, data = dnase,
                 alpha_0 = 501, omega_0 = 500, mu_0 = c(4.4, 0.25, 0.04), v_0 = 1)

# Print summary
summary(fit)
```

### Fitting a Model with Frailty

``` r
# Using dataset included in the package
data(simulation_frailty)

# Fit a survival model with shared frailty 
fit_frailty <- survregVB(formula = Surv(Time.15, delta.15) ~ x1 + x2, data = simulation_frailty,
                         alpha_0 = 3, omega_0 = 2, mu_0 = c(0, 0, 0), v_0 = 0.1,
                         lambda_0 = 3, eta_0 = 2, cluster = cluster)

# Print summary
summary(fit_frailty)
```

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/chengqianxian/survregVB/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chengqianxian/survregVB/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->
