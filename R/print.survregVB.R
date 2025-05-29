#' @method print survregVB
#' @export
print.survregVB <- function(x, digits = max(options()$digits - 4, 3),
                            signif.stars = FALSE, ...) {
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  if (!is.null(x$not_converged)) {
    cat("\nThe VB algorithm did not converge.\n")
  }

  cat("\nPosterior distributions of the regression coefficients (Beta):")
  cat("\nmu=\n")
  print(x$mu, digits = digits)
  cat("\nSigma=\n")
  printCoefmat(x$Sigma,
    digits = digits, signif.stars = signif.stars,
    P.values = TRUE, has.Pvalue = TRUE
  )

  cat("\nPosterior distribution of the scale parameter (b):")
  cat("\nalpha= ", x$alpha, "  omega= ", x$omega, "\n")

  # for models with shared frailty
  if (!is.null(x$clustered)) {
    cat("\nPosterior distribution of the random intercept (sigma_gamma squared):")
    cat("\nlambda= ", x$lambda, "  eta= ", x$eta, "\n")

    cat("\nPosterior distributions of the random effects for each cluster (gamma):")
    cat("\ntau=\n")
    print(x$tau, digits = digits)
    cat("\nsigma=\n")
    print(x$sigma, digits = digits)
  }

  cat("\nELBO= ", round(x$ELBO, digits), "\n")
  cat("\nNumber of iterations= ", x$iterations, "\n")

  omit <- x$na.action
  if (length(omit)) {
    cat("\nn=", x$n, " (", naprint(omit), ")\n", sep = "")
  } else {
    cat("\nn=", x$n, "\n")
  }

  invisible(x)
}
