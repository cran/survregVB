## Parameter Calculations ==============================================

#' Calculates parameter \eqn{\alpha^*} of \eqn{q^*(b)} to optimize the
#' evidence based lower bound (ELBO) in \code{survregVB.fit} and
#' \code{survregVB.frailty.fit}.
#'
#' @inheritParams elbo
#' @returns Parameter \eqn{\alpha^*} of \eqn{q^*(b)}.
#'
#' @seealso \code{\link{survregVB.fit}}
#' @seealso \code{\link{survregVB.frailty.fit}}
alpha_star <- function(alpha_0, delta) {
  alpha_0 + sum(delta)
}

#' Calculates parameter \eqn{\omega^*} of \eqn{q^*(b)} to optimize the
#' evidence based lower bound (ELBO) in \code{survregVB.fit}.
#'
#' @inheritParams elbo
#' @returns Parameter \eqn{\omega^*} of \eqn{q^*(b)}.
#'
#' @seealso \code{\link{survregVB.fit}}
omega_star <- function(y, X, delta, omega_0, mu, expectation_b) {
  res <- 0
  for (i in seq_len(nrow(X))) {
    bz_i <- y[i] - sum(X[i, ] * mu)
    z_i <- bz_i / expectation_b
    phi <- ifelse(z_i <= -5, 0,
      ifelse(z_i <= -1.701, 0.0426,
        ifelse(z_i <= 0, 0.3052,
          ifelse(z_i <= 1.702, 0.6950,
            ifelse(z_i <= 5, 0.9574, 1)
          )
        )
      )
    )
    res <- res + (delta[i] - (1 + delta[i]) * phi) * bz_i
  }
  omega_0 - res
}

#' Calculates parameter \eqn{\mu^*} of \eqn{q^*(\beta)} to optimize the
#' evidence based lower bound (ELBO) in \code{survregVB.fit}.
#'
#' @inheritParams elbo
#' @returns Parameter \eqn{\mu^*} of \eqn{q^*(\beta)}.
#'
#' @seealso \code{\link{survregVB.fit}}
mu_star <- function(y, X, delta, mu_0, v_0, alpha, omega, mu, Sigma,
                    expectation_b) {
  p <- ncol(X)
  inv_b <- alpha / omega
  inv_b_2 <- (alpha + alpha^2) / omega^2

  yX_matrix <- matrix(0, nrow = 1, ncol = p)
  for (i in seq_len(nrow(X))) {
    z_i <- (y[i] - sum(X[i, ] * mu)) / expectation_b
    rho <- ifelse(z_i <= -5, 0,
      ifelse(z_i <= -1.7, 0.1696,
        ifelse(z_i <= 1.7, 0.5,
          ifelse(z_i <= 5, 0.8303, 1)
        )
      )
    )
    zeta <- ifelse(z_i <= -5, 0,
      ifelse(z_i <= -1.7, 0.0189,
        ifelse(z_i <= 1.7, 0.1138,
          ifelse(z_i <= 5, 0.0190, 0)
        )
      )
    )

    yX_matrix <- yX_matrix +
      (inv_b * (-delta[i] + (1 + delta[i]) * rho) +
        2 * inv_b_2 * (1 + delta[i]) * y[i] * zeta) * X[i, ]
  }

  (v_0 * mu_0 + yX_matrix) %*% Sigma
}

#' Calculates parameter \eqn{\Sigma^*} of \eqn{q^*(\beta)} to optimize the
#' evidence based lower bound (ELBO) in \code{survregVB.fit}.
#'
#' @inheritParams elbo
#' @returns Parameter \eqn{\Sigma^*} of \eqn{q^*(\beta)}.
#'
#' @seealso \code{\link{survregVB.fit}}
Sigma_star <- function(y, X, delta, v_0, alpha, omega, mu, expectation_b) {
  p <- ncol(X)
  X_matrix <- matrix(0, nrow = p, ncol = p)
  for (i in seq_len(nrow(X))) {
    z_i <- (y[i] - sum(X[i, ] * mu)) / expectation_b
    zeta <- ifelse(z_i <= -5, 0,
      ifelse(z_i <= -1.7, 0.0189,
        ifelse(z_i <= 1.7, 0.1138,
          ifelse(z_i <= 5, 0.0190, 0)
        )
      )
    )
    X_matrix <- X_matrix + (1 + delta[i]) * zeta * (X[i, ] %*% t(X[i, ]))
  }

  inv_b_2 <- (alpha + alpha^2) / omega^2
  Sigma_inv <- diag(v_0, p) + 2 * inv_b_2 * X_matrix

  if (nrow(Sigma_inv) == 1) {
    matrix(1 / Sigma_inv, nrow = 1)
  } else {
    solve(Sigma_inv)
  }
}

## With cluster index ==================================================

#' Calculates parameter \eqn{\omega^*} of \eqn{q^*(b)} to optimize the evidence
#' based lower bound (ELBO) in \code{survregVB.frailty.fit}.
#'
#' @inheritParams elbo_cluster
#' @returns Parameter \eqn{\omega^*} of \eqn{q^*(b)}.
#'
#' @seealso \code{\link{survregVB.frailty.fit}}
omega_star_cluster <- function(y, X, delta, omega_0, mu, tau, expectation_b,
                               cluster) {
  y_cluster <- y - tau[cluster]
  res <- 0
  for (i in seq_len(nrow(X))) {
    z_i <- (y_cluster[i] - sum(X[i, ] * mu)) / expectation_b
    phi <- ifelse(z_i <= -5, 0,
      ifelse(z_i <= -1.701, 0.0426,
        ifelse(z_i <= 0, 0.3052,
          ifelse(z_i <= 1.702, 0.6950,
            ifelse(z_i <= 5, 0.9574, 1)
          )
        )
      )
    )
    bz_i <- y[i] - sum(X[i, ] * mu)
    res <- res + (delta[i] - (1 + delta[i]) * phi) * bz_i
  }
  omega_0 - res
}

#' Calculates parameter \eqn{\mu^*} of \eqn{q^*(\beta)} to optimize the
#' evidence based lower bound (ELBO) in \code{survregVB.frailty.fit}.
#'
#' @inheritParams elbo_cluster
#' @returns Parameter \eqn{\mu^*} of \eqn{q^*(\beta)}
#'
#' @seealso \code{\link{survregVB.frailty.fit}}
mu_star_cluster <- function(y, X, delta, mu_0, v_0, alpha, omega, mu, Sigma,
                            tau, expectation_b, cluster) {
  y <- y - tau[cluster]
  mu_star(y, X, delta, mu_0, v_0, alpha, omega, mu, Sigma, expectation_b)
}

#' Calculates parameter \eqn{\Sigma^*} of \eqn{q^*(\beta)} to optimize the
#' evidence based lower bound (ELBO) in \code{survregVB.frailty.fit}.
#'
#' @inheritParams elbo_cluster
#' @returns Parameter \eqn{\Sigma^*} of \eqn{q^*(\beta)}.
#'
#' @seealso \code{\link{survregVB.frailty.fit}}
Sigma_star_cluster <- function(y, X, delta, v_0, alpha, omega, mu, tau,
                               expectation_b, cluster) {
  y <- y - tau[cluster]
  Sigma_star(y, X, delta, v_0, alpha, omega, mu, expectation_b)
}

#' Calculates parameter \eqn{\sigma^{2*}_i} of \eqn{q^*(\gamma_i)} for
#' \eqn{i=1,...,K} clusters to optimize the evidence based lower bound
#' (ELBO) in \code{survregVB.frailty.fit}.
#'
#' @inheritParams elbo_cluster
#' @returns Parameter vector \eqn{\sigma^{2*}_i} of \eqn{q^*(\gamma_i)}
#'  for all clusters.
#'
#' @seealso \code{\link{survregVB.frailty.fit}}
sigma_squared_star <- function(y, X, delta, alpha, omega, mu, tau, lambda,
                               eta, expectation_b, cluster) {
  y <- y - tau[cluster]
  zeta <- numeric(nrow(X))

  for (i in seq_len(nrow(X))) {
    z_i <- (y[i] - sum(X[i, ] * mu)) / expectation_b
    zeta[i] <- ifelse(z_i <= -5, 0,
      ifelse(z_i <= -1.7, 0.0189,
        ifelse(z_i <= 1.7, 0.1138,
          ifelse(z_i <= 5, 0.0190, 0)
        )
      )
    )
  }

  inv_b_2 <- (alpha + alpha^2) / omega^2
  inv_sigma <- lambda / eta
  K <- length(unique(cluster))

  sigma <- numeric(K)
  for (k in 1:K) {
    delta_k <- delta[cluster == k]
    zeta_k <- zeta[cluster == k]
    sigma[k] <- 1 / (inv_sigma + 2 * inv_b_2 * sum((1 + delta_k) * zeta_k))
  }
  sigma
}

#' Calculates parameter \eqn{\tau^*_i} of \eqn{q^*(\gamma_i)} for
#' \eqn{i=1,...,K} clusters to optimize the evidence based lower bound
#' (ELBO) in \code{survregVB.frailty.fit}.
#'
#' @inheritParams elbo_cluster
#' @returns Parameter vector \eqn{\tau^*_i} of \eqn{q^*(\gamma_i)} for
#'  \eqn{i=1,...,K} clusters.
#'
#' @seealso \code{\link{survregVB.frailty.fit}}
tau_star <- function(y, X, delta, alpha, omega, mu, tau, sigma,
                     expectation_b, cluster) {
  y_cluster <- y - tau[cluster]
  n <- nrow(X)
  zeta <- numeric(n)
  rho <- numeric(n)

  for (i in 1:n) {
    # zeta_i is the coefficient of quadratic approximation
    z_i <- (y_cluster[i] - sum(X[i, ] * mu)) / expectation_b
    zeta[i] <- ifelse(z_i <= -5, 0,
      ifelse(z_i <= -1.7, 0.0189,
        ifelse(z_i <= 1.7, 0.1138,
          ifelse(z_i <= 5, 0.0190, 0)
        )
      )
    )
    rho[i] <- ifelse(z_i <= -5, 0,
      ifelse(z_i <= -1.7, 0.1696,
        ifelse(z_i <= 1.7, 0.5,
          ifelse(z_i <= 5, 0.8303, 1)
        )
      )
    )
  }

  inv_b <- alpha / omega
  inv_b_2 <- (alpha + alpha^2) / omega^2

  K <- length(unique(cluster))
  tau <- numeric(K)
  for (k in 1:K) {
    delta_k <- delta[cluster == k]
    zeta_k <- zeta[cluster == k]
    rho_k <- rho[cluster == k]
    y_k <- y[cluster == k]
    X_k <- X[cluster == k, , drop = FALSE]
    tau_k <- 0

    for (i in seq_len(nrow(X_k))) {
      tau_sub <- inv_b * (rho_k[i] + delta_k[i] * (rho_k[i] - 1)) +
        2 * inv_b_2 * (1 + delta_k[i]) * zeta_k[i] *
          (y_k[i] - sum(X_k[i, ] * mu))

      tau_k <- tau_k + tau_sub
    }
    tau[k] <- tau_k * sigma[k]
  }
  tau
}

#' Calculates parameter \eqn{\lambda^*} of \eqn{q^*(\sigma^2_{\gamma})} to
#' optimize the evidence based lower bound (ELBO) in
#' \code{survregVB.frailty.fit}.
#'
#' @inheritParams elbo_cluster
#' @param K The number of clusters.
#' @return Parameter \eqn{\lambda^*} of \eqn{q^*(\sigma^2_{\gamma})}.
#'
#' @seealso \code{\link{survregVB.frailty.fit}}
lambda_star <- function(lambda_0, K) {
  lambda_0 + K / 2
}

#' Calculates parameter \eqn{\eta^*} of \eqn{q^*(\sigma^2_{\gamma})} to
#' optimize the evidence based lower bound (ELBO) in
#' \code{survregVB.frailty.fit}.
#'
#' @inheritParams elbo_cluster
#' @return Parameter \eqn{\eta^*} of \eqn{q^*(\sigma^2_{\gamma})}.
#'
#' @seealso \code{\link{survregVB.frailty.fit}}
eta_star <- function(eta_0, tau, sigma) {
  eta_0 + 0.5 * sum(sigma + tau^2)
}
