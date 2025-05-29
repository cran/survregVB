#' @method print summary.survregVB
#' @export
print.summary.survregVB <- function(x, digits =
                                      max(options()$digits - 4, 3),
                                    signif.stars = FALSE, ...) {
  cat("Call:\n")
  dput(x$call)

  printCoefmat(round(x$estimates, digits),
    signif.stars = signif.stars,
    P.values = TRUE, has.Pvalue = TRUE
  )

  cat("\nELBO= ", round(x$ELBO, digits), "\n")
  cat("\nNumber of iterations= ", x$iterations, "\n")

  omit <- x$na.action
  if (length(omit)) {
    cat("\nn=", x$n, " (", naprint(omit), ")\n", sep = "")
  } else {
    cat("\nn=", x$n, "\n")
  }

  invisible(NULL)
}
