#' Variational Bayesian Analysis of Correlated Survival Data Using a
#' Log-Logistic Accelerated Failure Time Model
#'
#' Called by \code{survregVB} to do the actual parameter and ELBO computations
#' for correlated survival data with shared frailty (a random intercept).
#' This routine does no checking that the arguments are the proper length
#' or type.
#'
#' @name survregVB.frailty.fit
#'
#' @inheritParams survregVB
#' @param Y A `Surv` object containing 2 columns: time and event.
#' @param X A design matrix including covariates with first column of ones
#'  to represent the intercept.
#'
#' @returns A list containing results of the fit.
#'
#' @details
#' Implements the Variational Bayes algorithm with random intercepts proposed
#' in the paper "Fast variational bayesian inference for correlated survival
#' data: An application to invasive mechanical ventilation duration analysis".
#'
#' For right-censored survival time \eqn{T_{ij}} of the \eqn{j_{th}} subject
#' from the \eqn{i_{th}} cluster in the sample, in a sample, \eqn{i=1,...,K}
#' and \eqn{j=1,...,n_i}, the shared-frailty log-logistic AFT model is specified
#' as follows:
#'
#' \eqn{\log(T_{ij})=\gamma_i+X_{ij}^T\beta+b\epsilon_{ij}}, where
#' - \eqn{X_{ij}} is a column vector of length \eqn{p, p\ge2} containing
#'    \eqn{p-1} covariates and a constant one to incorporate the intercept
#' (i.e., \eqn{X_i=(1,x_{ij1},...,x_{ij(p-1)})^T}),
#' - \eqn{\beta} is the corresponding vector of coefficients for the fixed effects,
#' - \eqn{\gamma_i} is a random intercept for the \eqn{i_{th}} cluster,
#' - \eqn{\epsilon_{ij}} is a random variable following a standard logistic
#' distribution, and
#' - \emph{b} is a scale parameter.
#'
#' @export
#' @examples
#' fit <- survregVB.frailty.fit(
#'   X = matrix(c(rep(1, 75), simulation_frailty$x1, simulation_frailty$x2), nrow = 75),
#'   Y = survival::Surv(simulation_frailty$Time, simulation_frailty$delta),
#'   alpha_0 = 3,
#'   omega_0 = 2,
#'   mu_0 = c(0, 0, 0),
#'   v_0 = 0.1,
#'   lambda_0 = 3,
#'   eta_0 = 2,
#'   cluster = simulation_frailty$cluster
#' )
#'
#' @references Xian, C., Souza, C. P. E. de, He, W., Rodrigues, F. F.,
#'   & Tian, R. (2024). "Fast variational bayesian inference for correlated
#'   survival data: An application to invasive mechanical ventilation
#'   duration analysis." https://doi.org/10.48550/ARXIV.2408.00177
#' @seealso \code{\link{survregVB}}
survregVB.frailty.fit <- function(Y, X, alpha_0, omega_0, mu_0, v_0, lambda_0,
                                  eta_0, cluster, max_iteration = 100,
                                  threshold = 0.0001) {
  y <- log(Y[, 1])
  delta <- Y[, 2]
  cluster_names <- unique(cluster)
  cluster <- as.numeric(as.factor(cluster))
  n <- nrow(X)
  K <- length(unique(cluster))

  alpha <- alpha_star(alpha_0, delta) # fixed always
  omega <- omega_0
  lambda <- lambda_star(lambda_0, K)

  eta <- eta_0

  converged <- FALSE
  iteration <- 0

  expectation_b <- omega_0 / (alpha_0 - 1)
  curr_mu <- mu_0
  curr_tau <- rep(0, K)
  curr_elbo <- 0

  while (!converged && iteration <= max_iteration) {
    iteration <- iteration + 1
    Sigma <- Sigma_star_cluster(
      y, X, delta, v_0, alpha, omega, curr_mu,
      curr_tau, expectation_b, cluster
    )
    mu <- mu_star_cluster(
      y, X, delta, mu_0, v_0, alpha, omega, curr_mu,
      Sigma, curr_tau, expectation_b, cluster
    )
    sigma <- sigma_squared_star(
      y, X, delta, alpha, omega, mu, curr_tau,
      lambda, eta, expectation_b, cluster
    )
    tau <- tau_star(
      y, X, delta, alpha, omega, mu, curr_tau, sigma,
      expectation_b, cluster
    )
    omega <- omega_star_cluster(
      y, X, delta, omega_0, mu, tau,
      expectation_b, cluster
    )
    eta <- eta_star(eta_0, tau, sigma)

    elbo <- elbo_cluster(
      y, X, delta, alpha_0, omega_0, mu_0, v_0,
      lambda_0, eta_0, alpha, omega, mu, Sigma, tau,
      sigma, lambda, eta, expectation_b, cluster
    )

    if (abs(elbo - curr_elbo) > threshold &&
      sum(abs(mu - curr_mu)) > threshold) {
      converged <- FALSE
    } else {
      converged <- TRUE
    }

    expectation_b <- omega / (alpha - 1)
    curr_elbo <- elbo
    curr_mu <- mu
    curr_tau <- tau
  }

  mu <- c(mu)
  names(mu) <- colnames(X)
  dimnames(Sigma) <- list(colnames(X), colnames(X))
  names(tau) <- cluster_names
  names(sigma) <- cluster_names

  return_list <- list(
    clustered = TRUE,
    ELBO = unname(elbo),
    alpha = alpha,
    omega = unname(omega),
    mu = mu,
    Sigma = Sigma,
    tau = tau,
    sigma = sigma,
    lambda = lambda,
    eta = eta,
    iterations = iteration,
    n = n
  )

  if (converged == FALSE) {
    warning(
      "The max iteration has been achieved and the algorithm has not converged\n"
    )
    return_list$not_converged <- TRUE
  }

  return_list
}
