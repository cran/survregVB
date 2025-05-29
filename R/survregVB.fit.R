#' Variational Bayesian Analysis of Survival Data Using a Log-Logistic
#' Accelerated Failure Time Model
#'
#' Called by \code{survregVB} to do the actual parameter and ELBO
#' computations. This routine does no checking that the arguments are the
#' proper length or type.
#'
#' @name survregVB.fit
#'
#' @inheritParams survregVB
#' @param Y A `Surv` object containing 2 columns: time and event.
#' @param X A design matrix including covariates with first column of ones
#'  to represent the intercept.
#' @returns A list containing results of the fit.
#'
#' @details
#' Implements the Variational Bayes algorithm proposed in the paper "Variational
#' Bayesian analysis of survival data using a log-logistic accelerated failure
#' time model."
#'
#' For right-censored survival time \eqn{T_i} of the \eqn{i_{th}} subject
#' in a sample, \eqn{i=1,...,n}, the log-logistic AFT model is specified
#' as follows:
#'
#' \eqn{\log(T_i)=X_i^T\beta+bz_i}, where
#' - \eqn{X_i} is a column vector of length \eqn{p, p\ge2} containing
#' \eqn{p-1} covariates and a constant one to incorporate the intercept
#' (i.e., \eqn{X_i=(1,x_{i1},...,x_{i(p-1)})^T}),
#' - \eqn{\beta} is the corresponding vector of coefficients for the fixed
#' effects,
#' - \eqn{z_i} is a random variable following a standard logistic
#' distribution, and
#' - \emph{b} is a scale parameter.
#'
#' @export
#' @references Xian, C., Souza, C. P. E. de, He, W., Rodrigues, F. F.,
#'   & Tian, R. (2024). "Variational Bayesian analysis of survival data
#'   using a log-logistic accelerated failure time model." Statistics and
#'   Computing, 34(2). https://doi.org/10.1007/s11222-023-10365-6
#' @examples
#' fit <- survregVB.fit(
#'   Y = survival::Surv(simulation_nofrailty$Time, simulation_nofrailty$delta),
#'   X = matrix(c(rep(1, 300), simulation_nofrailty$x1, simulation_nofrailty$x2), nrow = 300),
#'   alpha_0 = 11,
#'   omega_0 = 10,
#'   mu_0 = c(0, 0, 0),
#'   v_0 = 1
#' )
#'
#' @seealso \code{\link{survregVB}}
survregVB.fit <- function(Y, X, alpha_0, omega_0, mu_0, v_0,
                          max_iteration = 100, threshold = 0.0001) {
  y <- log(Y[, 1])
  delta <- Y[, 2]
  n <- nrow(X)

  alpha <- alpha_star(alpha_0, delta)
  omega <- omega_0
  mu <- mu_0

  converged <- FALSE
  iteration <- 0

  expectation_b <- omega / (alpha - 1)
  curr_mu <- mu_0
  curr_elbo <- 0

  while (!converged && iteration < max_iteration) {
    iteration <- iteration + 1
    Sigma <- Sigma_star(y, X, delta, v_0, alpha, omega, curr_mu, expectation_b)
    mu <- mu_star(
      y, X, delta, mu_0, v_0, alpha, omega, curr_mu, Sigma,
      expectation_b
    )
    omega <- omega_star(y, X, delta, omega_0, mu, expectation_b)

    elbo <- elbo(
      y, X, delta, alpha_0, omega_0, mu_0, v_0, alpha, omega,
      curr_mu, Sigma, expectation_b
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
  }

  mu <- c(mu)
  names(mu) <- colnames(X)
  dimnames(Sigma) <- list(colnames(X), colnames(X))

  return_list <- list(
    ELBO = unname(elbo),
    alpha = alpha,
    omega = unname(omega),
    mu = mu,
    Sigma = Sigma,
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
