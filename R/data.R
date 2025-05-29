## Datasets ============================================================

#' Subset of `rhDNase` from the `survival` package
#'
#' The `dnase` dataset is a subset of the `rhDNase` dataset from the
#' `survival` package.
#' It is included in this package under the LGPL (\eqn{\geq2}) license.
#'
#' @format A data frame with 767 observations on the following variables:
#' \describe{
#'   \item{trt}{treatment arm: 0=placebo, 1= rhDNase}
#'   \item{fev}{forced expriatory volume at enrollment, a measure of lung
#'    capacity}
#'   \item{infect}{an infection that required the use of intravenous
#'    antibiotics}
#'   \item{time}{difference between the date of entry into the study and
#'    the date of last follow-up capped at 169 days}
#' }
#'
#' @source `survival` package.
#'  \url{https://cran.r-project.org/package=survival}
"dnase"

#' Subset of GSE102287: African American (AA) Patients
#'
#' This dataset is a subset of the GSE102287 dataset that includes only
#' characteristics of patients who are identified as African American (AA).
#'
#' @format A data frame with 60 observations on selected patient
#'  characteristics:
#' \describe{
#'   \item{patient}{Patient identification number.}
#'   \item{age}{Patient age.}
#'   \item{Stage}{Lung cancer stage (I, II, III).}
#'   \item{time}{Survival time in days.}
#'   \item{gender}{Gender of the patient.}
#'   \item{smoking}{0 = Never smoked, 1 = Has smoked.}
#'   \item{status}{0 = Alive, 1 = Death due to lung cancer.}
#' }
#'
#' @source  Gene Expression Omnibus (GEO), Accession: GSE102287.
#'  \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102287}
#'
#' @references Mitchell, K. A., Zingone, A., Toulabi, L., Boeckelman, J.,
#' & Ryan, B. M. (2017). Comparative Transcriptome Profiling Reveals Coding
#' and Noncoding RNA Differences in NSCLC from African Americans and European
#' Americans. Clinical cancer research: an official journal of the American
#' Association for Cancer Research, 23(23), 7412–7425.
#' doi:10.1158/1078-0432.CCR-17-0527.
"lung_cancer"

#' Simulated data incorporating shared frailty effects to model clustered
#' time-to-event data.
#'
#' @format A dataframe with 75 observations grouped into 15 clusters, each
#'  with 5 individuals.
#'  \describe{
#'   \item{x1}{Continuous covariate from `N(1, 0.2^2)`}
#'   \item{x2}{Binary covariate from `Bernoulli(0.5)`}
#'   \item{Time}{True survival time}
#'   \item{Time.15}{Observed survival time accounting for uniformly distributed
#'    right censoring time from `uniform(0,u)`}
#'   \item{delta}{Event indicator for uncensored data (always 1 in this
#'    simulation.)}
#'   \item{delta.15}{Event indicator after censoring (1 = event, 0 =
#'    censored).}
#'   \item{cluster}{Cluster ID (1–15), indicating group-level frailty}
#'. }
#'   @references Xian, C., Souza, C. P. E. de, He, W., Rodrigues, F. F.,
#'   & Tian, R. (2024). Fast variational bayesian inference for correlated
#'   survival data: An application to invasive mechanical ventilation
#'   duration analysis. https://doi.org/10.48550/ARXIV.2408.00177
"simulation_frailty"

#' Simulated data without shared frailty effects to model unclustered
#' time-to-event data.
#'
#' @format A dataframe with 300 observations.
#'  \describe{
#'   \item{x1}{Continuous covariate from `N(1, 0.2^2)`}
#'   \item{x2}{Binary covariate from `Bernoulli(0.5)`}
#'   \item{Time}{True survival time}
#'   \item{Time.10}{Observed survival time accounting for uniformly distributed
#'    right censoring time from `uniform(0,48)`}
#'   \item{Time.30}{Observed survival time accounting for uniformly distributed
#'    right censoring time from `uniform(0,17)`}
#'   \item{delta}{Event indicator for uncensored data (always 1 in this
#'    simulation.)}
#'   \item{delta.10}{Event indicator for T.10 (1 = event, 0 =
#'    censored).}
#'   \item{delta.30}{Event indicator for T.30 (1 = event, 0 =
#'    censored).}
#'   }
#'   @references Xian, C., Souza, C. P. E. de, He, W., Rodrigues, F. F.,
#'   & Tian, R. (2024). Variational Bayesian analysis of survival data
#'   using a log-logistic accelerated failure time model. Statistics and
#'   Computing, 34(2). https://doi.org/10.1007/s11222-023-10365-6
"simulation_nofrailty"
