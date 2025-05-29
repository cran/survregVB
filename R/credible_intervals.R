#' Calculates the credible interval for a posterior distribution,
#' \eqn{q^*(\beta)}, a \eqn{N_p(\mu,\Sigma)} density function.
#'
#' @inheritParams elbo
#' @param ci The significance level. (Default:0.95).
#' @returns Matrix containing the credible intervals for \eqn{q^*(\beta)}.
#'
#' @noRd
beta_ci <- function(mu, Sigma, ci = 0.95) {
  k <- length(mu)
  lower_bounds <- numeric(k)
  upper_bounds <- numeric(k)
  for (i in 1:k) {
    lower_bounds[i] <-
      qnorm((1 - ci) / 2, mu[i], sqrt(diag(Sigma)[i]))
    upper_bounds[i] <-
      qnorm(1 - (1 - ci) / 2, mu[i], sqrt(diag(Sigma)[i]))
  }
  cbind(CI.Lower = lower_bounds, CI.Upper = upper_bounds)
}

#' Calculates the credible interval for a posterior distribution,
#' \eqn{q^*(b)}, an \eqn{\text{Inverse-Gamma}(\alpha,\omega)} density
#' function.
#'
#' @inheritParams elbo
#' @param ci The significance level. (Default:0.95).
#' @returns The credible interval for \emph{b}.
#'
#' @importFrom bayestestR hdi
#' @importFrom invgamma rinvgamma
#' @noRd
b_ci <- function(alpha, omega, ci = 0.95, seed = 100) {
  set.seed(seed)
  posterior <- rinvgamma(100000, alpha, omega)
  lower <- hdi(posterior, ci)$CI_low
  upper <- hdi(posterior, ci)$CI_high
  (c(lower, upper))
}
