#' Variational Bayesian Analysis of Survival Data Using a Log-Logistic
#' Accelerated Failure Time Model
#'
#' Applies a mean-field Variational Bayes (VB) algorithm to infer the
#' parameters of an accelerated failure time (AFT) survival model with
#' right-censored survival times following a log-logistic distribution.
#'
#' @name survregVB
#'
#' @param formula A formula object, with the response on the left of a `~`
#'  operator, and the covariates on the right. The response must be a survival
#'  object of type `right`, as returned by the \code{Surv} function.
#' @param data A `data.frame` in which to interpret the variables named in
#'  the `formula` and `cluster` arguments.
#' @param alpha_0 The shape hyperparameter \eqn{\alpha_0} of the prior
#'  distribution of the scale parameter, \emph{b}.
#' @param omega_0 The shape hyperparameter \eqn{\omega_0} of the prior
#'  distribution of the scale parameter, \emph{b}.
#' @param mu_0 Hyperparameter \eqn{\mu_0}, a vector of means, of the prior
#'  distribution of the vector of coefficients, \eqn{\beta}.
#' @param v_0 The precision (inverse variance) hyperparameter \eqn{v_0},
#'  of the prior distribution of the vector of coefficients, \eqn{\beta}.
#' @param lambda_0 The shape hyperparameter \eqn{\lambda_0} of the prior
#'  distribution of the frailty variance, \eqn{\sigma_\gamma^2}.
#' @param eta_0 The scale hyperparameter \eqn{\eta_0} of the prior distribution
#'  of the frailty variance, \eqn{\sigma_\gamma^2}.
#' @param na.action A missing-data filter function, applied to the
#'  \code{model.frame}, after any subset argument has been used.
#'  (Default:\code{options()$na.action}).
#' @param cluster An optional variable which clusters the observations to
#'  introduce shared frailty for correlated survival data.
#' @param max_iteration The maximum number of iterations for the variational
#'  inference optimization. If reached, iteration stops. (Default:100)
#' @param threshold The convergence threshold for the evidence based lower
#'  bound (ELBO) optimization. If the difference between the current and
#'  previous ELBO's is smaller than this threshold, iteration stops.
#'  (Default:0.0001)
#'
#' @returns An object of class \code{survregVB}.
#'
#' @details
#' The goal of \code{survregVB} is to maximize the evidence lower bound
#' (ELBO) to approximate posterior distributions of the AFT model parameters
#' using the VB algorithms with and without shared frailty proposed in Xian
#' et al. (2024) <doi:10.1007/s11222-023-10365-6> and
#' <doi:10.48550/ARXIV.2408.00177> respectively.
#'
#' @examples
#' # Data frame containing survival data
#' fit <- survregVB(
#'   formula = survival::Surv(time, infect) ~ trt + fev,
#'   data = dnase,
#'   alpha_0 = 501,
#'   omega_0 = 500,
#'   mu_0 = c(4.4, 0.25, 0.04),
#'   v_0 = 1,
#'   max_iteration = 100,
#'   threshold = 0.0005
#' )
#' summary(fit)
#'
#' # Call the survregVB function with shared frailty
#' fit2 <- survregVB(
#'   formula = survival::Surv(Time.15, delta.15) ~ x1 + x2,
#'   data = simulation_frailty,
#'   alpha_0 = 3,
#'   omega_0 = 2,
#'   mu_0 = c(0, 0, 0),
#'   v_0 = 0.1,
#'   lambda_0 = 3,
#'   eta_0 = 2,
#'   cluster = cluster,
#'   max_iteration = 100,
#'   threshold = 0.01
#' )
#' summary(fit2)
#' @import stats
#'
#' @export
#' @references
#'   Xian, C., Souza, C. P. E. de, He, W., Rodrigues, F. F.,
#'   & Tian, R. (2024). "Variational Bayesian analysis of survival data
#'   using a log-logistic accelerated failure time model." Statistics and
#'   Computing, 34(2). https://doi.org/10.1007/s11222-023-10365-6
#'
#'   Xian, C., Souza, C. P. E. de, He, W., Rodrigues, F. F.,
#'   & Tian, R. (2024). "Fast variational bayesian inference for correlated
#'   survival data: An application to invasive mechanical ventilation
#'   duration analysis." https://doi.org/10.48550/ARXIV.2408.00177
#'
#' @seealso \code{\link{survregVB.object}}
survregVB <- function(formula, data, alpha_0, omega_0, mu_0, v_0,
                      lambda_0, eta_0, na.action, cluster,
                      max_iteration = 100, threshold = 0.0001) {
  Call <- match.call() # save a copy of the call

  if (missing(formula)) stop("a formula argument is required")
  if (is.list(formula)) {
    stop("formula argument cannot be a list")
  }

  Terms <- if (missing(data)) {
    terms(formula)
  } else {
    terms(formula, data = data)
  }

  defined <- if (missing(cluster)) {
    c("alpha_0", "omega_0", "mu_0", "v_0")
  } else {
    c("alpha_0", "omega_0", "mu_0", "v_0", "lambda_0", "eta_0")
  }
  passed <- match(defined, names(Call), nomatch = 0)
  missing <- which(passed == 0)
  if (length(missing) > 0) {
    stop(paste("missing value(s) for",
               paste(defined[missing], collapse = ", ")))
  }

  indx <- match(c("formula", "data", "cluster", "na.action"),
    names(Call),
    nomatch = 0
  )

  temp <- Call[c(1, indx)] # only keep the arguments we wanted
  temp[[1L]] <- quote(stats::model.frame) # change the function called

  temp$formula <- if (missing(data)) {
    terms(formula)
  } else {
    terms(formula, data = data)
  }

  m <- eval(temp, parent.frame())

  Terms <- attr(m, "terms")
  Y <- model.extract(m, "response")

  if (!inherits(Y, "Surv")) {
    stop("response must be a survival object")
  }
  type <- attr(Y, "type")
  if (type != "right") {
    stop("only Survival objects of type right are supported")
  }

  X <- model.matrix(Terms, m)
  if (ncol(X) != length(mu_0)) {
    stop("the length of mu_0 must match the number of covariates")
  }

  cluster <- model.extract(m, "cluster")

  if (length(cluster)) {
    result <- survregVB.frailty.fit(
      Y, X, alpha_0, omega_0, mu_0, v_0,
      lambda_0, eta_0, cluster,
      max_iteration, threshold
    )
  } else {
    result <- survregVB.fit(
      Y, X, alpha_0, omega_0, mu_0, v_0, max_iteration,
      threshold
    )
  }

  na.action <- attr(m, "na.action")
  if (length(na.action)) result$na.action <- na.action
  result$call <- Call
  class(result) <- "survregVB"

  result
}

#' Variational Bayes Accelererated Failure Time Survival Model Object
#'
#' This class of objects is returned by the survregVB function to represent
#' a fitted parametric log-logistic accelerated failure time (AFT) survival
#' model. Objects of this class have methods for the functions \code{print}
#' and \code{summary}.
#'
#' @name survregVB.object
#' @aliases print.survregVB
#'
#' @details
#' For approximate posterior distributions:
#' - \eqn{q^*(\beta)}, a \eqn{N_p(\mu^*,\Sigma^*)} density function, and
#' - \eqn{q^*(b)}, an \eqn{\text{Inverse-Gamma}(\alpha^*,\omega^*)}
#'    density function,
#'
#' the components of this class are:
#' \itemize{
#'   \item \code{ELBO}: The final value of the Evidence Lower Bound (ELBO)
#'    at the last iteration.
#'   \item \code{alpha}: The shape parameter \eqn{\alpha^*} of \eqn{q^*(b)}.
#'   \item \code{omega}: The scale parameter \eqn{\omega^*} of \eqn{q^*(b)}.
#'   \item \code{mu}: Parameter \eqn{\mu^*} of \eqn{q^*(\beta)}, a vector
#'    of means.
#'   \item \code{Sigma}: Parameter \eqn{\Sigma^*} of \eqn{q^*(\beta)}, a
#'    covariance matrix.
#'   \item \code{na.action}: A missing-data filter function, applied to the
#'    \code{model.frame}, after any subset argument has been used.
#'   \item \code{iterations}: The number of iterations performed by the VB
#'    algorithm: before converging or reaching \code{max_iteration}.
#'   \item \code{n}: The number of observations.
#'   \item \code{call}: The function call used to invoke the \code{survregVB}
#'    method.
#'   \item \code{not_converged}: A boolean indicating if the algorithm
#'    converged.
#'   \itemize{
#'     \item \code{TRUE}: If the algorithm did not converge prior to
#'      achieving `max_iteration`.
#'     \item \code{NULL}: If the algorithm converged successfully.
#'   }
#' }
#'
#' If \code{survregVB} was called with shared frailty (with the `cluster`
#' argument), for approximate posterior distributions:
#' - \eqn{q^*(\sigma^2_\gamma)}, an \eqn{\text{Inverse-Gamma}(\lambda^*,\eta^*)}
#'   density function,
#' - \eqn{q^*(\gamma_i)}, a \eqn{N(\tau^*_i,\sigma^{2*}_i)} density function,
#'   for \eqn{i=1,...,K} clusters, and
#'
#' the additional components are present:
#'
#' \itemize{
#'   \item \code{lambda}: The shape parameter \eqn{\lambda^*} of
#'    \eqn{q^*(\sigma^2_\gamma)}.
#'   \item \code{eta}: The scale parameter \eqn{\eta^*} of
#'    \eqn{q^*(\sigma^2_\gamma)}.
#'   \item \code{tau}: Parameter \eqn{\tau^*_i} of \eqn{q^*(\gamma_i)}, a
#'    vector of means.
#'   \item \code{sigma}: Parameter \eqn{\sigma^{2*}_i} of \eqn{q^*(\gamma_i)},
#'    a vector of variance.
#'  }
#'
NULL
