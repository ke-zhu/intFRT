#' Determine the Adaptive Threshold for Conformal Selective Borrowing
#'
#' This function adaptively determines the threshold (gamma) to minimize the
#' Mean Squared Error (MSE) when performing Conformal Selective Borrowing.
#'
#' @param Y Numeric vector of outcomes.
#' @param A Binary treatment indicator vector (1 for treatment, 0 for control).
#' @param S Binary indicator vector for source of data (1 for randomized trial
#'   data, 0 for external control data).
#' @param X Matrix of covariates.
#' @param family A description of the outcome type. Default is "gaussian" for
#'   continuous outcomes, but can also be "binomial" for binary outcomes.
#' @param gamma_grid Numeric vector representing the grid of potential gamma
#'   values to search over. Default is a sequence from 0 to 1 with 11 grids.
#' @param measure The performance metric to be minimized. Default is "mse_hat".
#'   Other options include "var_hat" and "bias2_hat", but these are only for
#'   comparison purposes and not recommended for actual selection.
#' @param n_rep_gamma Integer specifying the number of data-splitting iterations
#'   for estimating the measure given a certain gamma. Default is 100.
#' @param parallel Logical value indicating whether to use parallel computing.
#'   Default is `FALSE`.
#' @param n_cores Integer specifying the number of cores to use for parallel
#'   computing. Default is the number of available physical cores on the
#'   machine, as determined by `parallel::detectCores(logical = FALSE)`.
#' @param ... Additional arguments passed to the internal `ec_borrow` function.
#'
#' @return A numeric value representing the selected gamma that minimizes the
#'   specified performance measure (default is MSE).
#'
#' @examples
#' # This example illustrates the use of the Fisher Randomization Test (FRT)
#' # and different borrowing methods for hybrid controlled trials.
#'
#' library(intFRT)
#'
#' # Simulate data for a hybrid controlled trial
#' set.seed(1)
#' n_rct <- 50  # Number of observations in the randomized controlled trial
#' n_ec <- 100  # Number of external controls
#' n <- n_rct + n_ec  # Total number of observations
#'
#' # Covariates (2 covariates, uniformly distributed)
#' X <- matrix(runif(n * 2), n, 2)
#'
#' # Treatment assignment (1 = treatment, 0 = control)
#' A <- c(rep(0:1, each = n_rct / 2), rep(0, n_ec))
#'
#' # Data source indicator (1 = randomized trial, 0 = external control)
#' S <- c(rep(1, n_rct), rep(0, n_ec))
#'
#' # Generate potential outcomes (continuous)
#' Y1 <- 1 + 2 * X[,1] + 0.1 * X[,2] + rnorm(n)
#' Y0 <- 2 * X[,1] + 0.1 * X[,2] + rnorm(n)
#'
#' # Introduce bias in half of the external controls
#' id_biased_EC <- tail(which(S == 0), n = n_ec / 2)
#' Y0[id_biased_EC] <- Y0[id_biased_EC] - 10
#'
#' # Observed outcome
#' Y <- A * Y1 + (1 - A) * Y0
#'
#' # Compute adaptive gamma (with a small n_rep_gamma for illustration)
#' ada_g <- compute_ada_gamma(
#'   Y, A, S, X,
#'   # Use a small n_rep_gamma for fast illustration;
#'   # recommend n_rep_gamma = 100 for more stable results
#'   n_rep_gamma = 10
#' )
#' ada_g
#'
#' @import tibble
#' @import dplyr
#' @import purrr
#' @export
compute_ada_gamma <- function(Y, A, S, X,
                              family = "gaussian",
                              gamma_grid = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
                              measure = "mse_hat", n_rep_gamma = 100,
                              parallel = F,
                              n_cores = parallel::detectCores(logical = FALSE),
                              ...) {
  if (!parallel) {n_cores <- 1}
  res_grid <- parallel::mclapply(gamma_grid, function(g) {
    dat_rct <- tibble(Y, A, S, X) %>% filter(S == 1)
    dat_ec <- tibble(Y, A, S, X) %>% filter(S == 0)
    n_rt1 <- floor(0.5 * sum(dat_rct$A == 1))
    n_rc1 <- floor(0.5 * sum(dat_rct$A == 0))
    est_rep <- map(1:n_rep_gamma, ~ {
      fold1_id <- c(
        sample(which(dat_rct$A == 1), size = n_rt1),
        sample(which(dat_rct$A == 0), size = n_rc1)
      )
      dat1 <- bind_rows(dat_rct[fold1_id,], dat_ec)
      dat2 <- bind_rows(dat_rct[-fold1_id,], dat_ec)
      est_csb <- ec_borrow(
        dat1$Y, dat1$A, dat1$S, dat1$X,
        "Conformal Selective Borrow AIPW", family, n_fisher = NULL,
        gamma_sel = g, ...
      )$res$est[1]
      est_nb <- ec_borrow(
        dat2$Y, dat2$A, dat2$S, dat2$X,
        "No Borrow AIPW", family, n_fisher = NULL
      )$res$est[1]
      lst(est_csb, est_nb)
    })
    var_hat <- map_dbl(est_rep, ~ {.$est_csb}) %>% var
    bias2_hat <- mean(map_dbl(est_rep, ~ {.$est_csb - .$est_nb}), na.rm = T)^2
    mse_hat <- var_hat + bias2_hat
    # output
    cat(paste0("For gamma_sel = ", g,
               ", MSE = ", mse_hat, "\n\n"))
    lst(mse_hat, var_hat, bias2_hat)
  }, mc.cores = n_cores)
  gamma_grid[which.min(map_dbl(res_grid, measure))]
}
