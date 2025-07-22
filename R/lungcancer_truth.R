#' Ground-Truth Potential Survival Outcomes for `lungcancer`
#'
#' A companion data set to \code{lungcancer}, containing the true potential
#' survival times under treatment and control for each individual.
#'
#' @format A data frame with n rows and 3 variables:
#' \describe{
#'   \item{patid}{Patient ID (numeric, matched to \code{lungcancer})}
#'   \item{T1}{Potential survival time (in years) under treatment}
#'   \item{T0}{Potential survival time (in years) under control}
#' }
#'
#' @details
#' This data set provides uncensored, ground-truth potential survival outcomes
#' for individuals in the \code{lungcancer} data set. It can be used to evaluate
#' the performance of survival models or causal inference methods in simulation
#' settings.
#'
#' For covariate definitions, cohort structure, and data generation details,
#' refer to the documentation for \code{lungcancer}.
#'
#' @examples
#' summary(lungcancer)
#' summary(lungcancer_truth)
#'
#' # library(summarytools)
#' # view(dfSummary(lungcancer))
#' # view(dfSummary(lungcancer_truth))
"lungcancer_truth"
