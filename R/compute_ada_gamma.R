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
#'   continuous outcomes.
#' @param gamma_grid Numeric vector representing the grid of potential gamma
#'   values to search over. Default is a sequence from 0 to 1 with 11 grids.
#' @param n_rep_gamma Integer specifying the number of resampling iterations
#'   for estimating the MSE given a certain gamma. Default is `10`.
#' @param parallel Logical value indicating whether to use parallel computing.
#'   Default is `FALSE`.
#' @param n_cores Integer specifying the number of cores to use for parallel
#'   computing. Default is the number of available physical cores on the
#'   machine, as determined by `future::availableCores(logical = FALSE)`.
#' @param opt Character specifying the criterion to determine gamma. Options are
#' "power" or "mse". Default is "mse".
#' @param sig_level Numeric value indicating the significance level for
#' hypothesis tests when `opt = "power"`. Default is 0.05.
#' @param ... Additional arguments passed to the internal `ec_borrow` function.
#'
#' @return A numeric value representing the selected gamma that minimizes the
#'   empirical MSE.
#'
#' @examples
#' # This example illustrates the use of the Fisher Randomization Test (FRT)
#' # and different borrowing methods for hybrid controlled trials.
#' \dontrun{
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
#' # Compute adaptive gamma
#' ada_g <- compute_ada_gamma(
#'   Y, A, S, X,
#'   # Tuning with 20 replications for illustration purposes
#'   # Recommend setting `n_rep_gamma = 100` or higher with parallel computing
#'   n_rep_gamma = 20
#' )
#' ada_g
#' }
#'
#' @import tibble
#' @import dplyr
#' @import purrr
#' @export
compute_ada_gamma <- function(Y, A, S, X,
                              family = "gaussian",
                              gamma_grid = seq(0, 1, by = 0.1),
                              n_rep_gamma = 100,
                              parallel = FALSE,
                              n_cores = NULL,
                              opt = "mse",
                              sig_level = 0.05,
                              ...) {
  if (parallel) {
    if (is.null(n_cores)) {
      n_cores <- future::availableCores(logical = FALSE)
    }
    future::plan(multisession, workers = n_cores)
  }

  if (any(gamma_grid == 1)) {
    id_nb <- which(gamma_grid == 1)
  } else {
    gamma_grid <- c(gamma_grid, 1)
    id_nb <- which(gamma_grid == 1)
  }

  # data
  dat_rct <- tibble(Y, A, S, X) %>% filter(S == 1)
  dat_ec <- tibble(Y, A, S, X) %>% filter(S == 0)
  dat_full <- bind_rows(dat_rct, dat_ec)
  n_rct <- nrow(dat_rct)
  n_full <- nrow(dat_full)

  # main algorithm, the other two only for testing
  if (opt == "mse" & !is.null(n_rep_gamma)) {
    # bootstrap data sets
    dat_rep <- map(1:n_rep_gamma, ~ {
      dat_rct_boot <- dat_rct %>%
        group_by(A) %>%
        slice_sample(prop = 1, replace = T)
      bind_rows(dat_rct_boot, dat_ec)
    })
    # estimators
    if (parallel) {
      est_grid <- furrr::future_map(gamma_grid, function(g) {
        est_one <- tryCatch({
          ec_borrow(
            dat_full$Y, dat_full$A, dat_full$S, dat_full$X,
            "Conformal Selective Borrow AIPW", family, n_fisher = NULL,
            gamma_sel = g, ...
          )$res$est[1]
        }, error = function(e) NA)

        if (!is.na(est_one)) {
          est_rep <- map(dat_rep, ~ {
            est_csb <- tryCatch({
              ec_borrow(
                .$Y, .$A, .$S, .$X,
                "Conformal Selective Borrow AIPW", family, n_fisher = NULL,
                gamma_sel = g, ...
              )$res$est[1]
            }, error = function(e) est_one)

            est_csb
          })
        } else {
          est_rep <- NA
        }

        lst(est_one, est_rep)
      }, .options = furrr::furrr_options(seed = TRUE))
    } else {
      est_grid <- map(gamma_grid, function(g) {
        est_one <- tryCatch({
          ec_borrow(
            dat_full$Y, dat_full$A, dat_full$S, dat_full$X,
            "Conformal Selective Borrow AIPW", family, n_fisher = NULL,
            gamma_sel = g, ...
          )$res$est[1]
        }, error = function(e) NA)

        if (!is.na(est_one)) {
          est_rep <- map(dat_rep, ~ {
            est_csb <- tryCatch({
              ec_borrow(
                .$Y, .$A, .$S, .$X,
                "Conformal Selective Borrow AIPW", family, n_fisher = NULL,
                gamma_sel = g, ...
              )$res$est[1]
            }, error = function(e) est_one)

            est_csb
          })
        } else {
          est_rep <- NA
        }

        lst(est_one, est_rep)
      })
    }

    # MSE
    res_grid <- map2(est_grid, gamma_grid, function(est_g, g) {
      if (!is.na(est_g$est_one)) {
        var_hat <- map_dbl(est_g$est_rep, ~ .) %>% var
        if (g == 1) {
          bias2_hat <- 0
        } else {
          bias2_hat <- max(
            (est_g$est_one - est_grid[[id_nb]]$est_one)^2 -
              map2_dbl(est_g$est_rep, est_grid[[id_nb]]$est_rep,
                       function(x, y) {x - y}) %>% var,
            0
          )
        }
        mse_hat <- bias2_hat + var_hat

        # print
        cat(sprintf("Gamma = %.2f | MSE = %.4f\n", g, mse_hat))
        lst(mse_hat, bias2_hat, var_hat)
      } else {
        # print
        cat(sprintf("Gamma = %.2f | Error\n", g))
        lst(mse_hat = Inf, bias2_hat = Inf, var_hat = Inf)
      }
    })
    # output
    out <- gamma_grid[which.min(map_dbl(res_grid, "mse_hat"))]
  }

  # for testing
  if (opt == "mse" & is.null(n_rep_gamma)) {
    # sandwich variance estimator, not recommend for this purpose
    # est
    if (parallel) {
      est_grid <- furrr::future_map(gamma_grid, function(g) {
        fit <- ec_borrow(
          dat_full$Y, dat_full$A, dat_full$S, dat_full$X,
          "Conformal Selective Borrow AIPW", family, n_fisher = NULL,
          gamma_sel = g, ...
        )
        est_one <- fit$out$est
        est_se <- fit$out$se
        est_d <- fit$out$d[[1]]
        lst(est_one, est_se, est_d)
      }, .options = furrr::furrr_options(seed = TRUE))
    } else {
      est_grid <- map(gamma_grid, function(g) {
        fit <- ec_borrow(
          dat_full$Y, dat_full$A, dat_full$S, dat_full$X,
          "Conformal Selective Borrow AIPW", family, n_fisher = NULL,
          gamma_sel = g, ...
        )
        est_one <- fit$out$est
        est_se <- fit$out$se
        est_d <- fit$out$d[[1]]
        lst(est_one, est_se, est_d)
      })
    }

    # MSE
    res_grid <- map2(est_grid, gamma_grid, function(est_g, g) {
      d_g <- est_g$est_d
      d_0 <- est_grid[[id_nb]]$est_d
      d_dif <- d_g - c(d_0, rep(0, length(d_g) - length(d_0)))
      var_dif_hat <- sum((d_dif - mean(d_dif))^2) / n_rct / n_full
      #var_hat <- sum((d_g - mean(d_g))^2) / n_rct / n_full
      var_hat <- est_g$est_se^2
      if (g == 1) {
        bias2_hat <- 0
      } else {
        bias2_hat <- max(
          (est_g$est_one - est_grid[[id_nb]]$est_one)^2 - var_dif_hat,
          0
        )
      }
      mse_hat <- bias2_hat + var_hat

      # output
      cat(paste0("For gamma_sel = ", g, ", MSE = ", mse_hat, "\n\n"))
      lst(mse_hat, bias2_hat, var_hat)
    })
    # output
    out <- gamma_grid[which.min(map_dbl(res_grid, "mse_hat"))]
  }

  # for testing
  if (opt == "power") {
    if (parallel) {
      power_hat <- furrr::future_map(gamma_grid, function(g) {
        # distribution of T under H0
        res_H0 <- ec_borrow(
          dat_full$Y, dat_full$A, dat_full$S, dat_full$X,
          "Conformal Selective Borrow AIPW", family, n_fisher = n_rep_gamma,
          gamma_sel = g, output_frt = TRUE
        )
        T_H0 <- res_H0$out_frt$est
        # CATE estimation
        #taux <- mean(dat_rct$Y[dat_rct$A == 1]) - mean(dat_rct$Y[dat_rct$A == 0])
        cate_fit <- grf::causal_forest(X = dat_rct$X, Y = dat_rct$Y,
                                       W = dat_rct$A)
        cate <- predict(cate_fit)
        taux <- cate$predictions
        # impute potential outcome under H1
        dat_imp_H1 <- dat_full %>%
          mutate(
            Y1 = case_when(
              A == 1 ~ Y,
              A == 0 ~ Y + taux
            ),
            Y0 = case_when(
              A == 0 ~ Y,
              A == 1 ~ Y - taux
            )
          )
        # distribution of T under H0
        T_H1 <- map_dbl(1:n_rep_gamma, ~ {
          dat_H1 <- dat_imp_H1 %>%
            mutate(
              A = {A[S == 1] <- sample(A[S == 1]); A},
              Y = A * Y1 + (1 - A) * Y0
            )
          res_H1 <- ec_borrow(
            dat_H1$Y, dat_H1$A, dat_H1$S, dat_H1$X,
            "Conformal Selective Borrow AIPW", family, n_fisher = NULL,
            gamma_sel = g, output_frt = TRUE
          )
          res_H1$res$est[1]
        })
        # bind_rows(
        #   tibble(type = "H0", t = T_H0),
        #   tibble(type = "H1", t = T_H1)
        # ) %>%
        #   ggplot(aes(t, fill = type)) +
        #   geom_histogram()
        critical_value <- quantile(abs(T_H0), probs = 1 - sig_level)
        power_hat_g <- mean(T_H1 > critical_value) # estimated power
        # output
        power_hat_g
      }, .options = furrr::furrr_options(seed = TRUE)) %>% map_dbl(~.)
    } else {
      power_hat <- map(gamma_grid, function(g) {
        # distribution of T under H0
        res_H0 <- ec_borrow(
          dat_full$Y, dat_full$A, dat_full$S, dat_full$X,
          "Conformal Selective Borrow AIPW", family, n_fisher = n_rep_gamma,
          gamma_sel = g, output_frt = TRUE
        )
        T_H0 <- res_H0$out_frt$est
        # CATE estimation
        #taux <- mean(dat_rct$Y[dat_rct$A == 1]) - mean(dat_rct$Y[dat_rct$A == 0])
        cate_fit <- grf::causal_forest(X = dat_rct$X, Y = dat_rct$Y,
                                       W = dat_rct$A)
        cate <- predict(cate_fit)
        taux <- cate$predictions
        # impute potential outcome under H1
        dat_imp_H1 <- dat_full %>%
          mutate(
            Y1 = case_when(
              A == 1 ~ Y,
              A == 0 ~ Y + taux
            ),
            Y0 = case_when(
              A == 0 ~ Y,
              A == 1 ~ Y - taux
            )
          )
        # distribution of T under H0
        T_H1 <- map_dbl(1:n_rep_gamma, ~ {
          dat_H1 <- dat_imp_H1 %>%
            mutate(
              A = {A[S == 1] <- sample(A[S == 1]); A},
              Y = A * Y1 + (1 - A) * Y0
            )
          res_H1 <- ec_borrow(
            dat_H1$Y, dat_H1$A, dat_H1$S, dat_H1$X,
            "Conformal Selective Borrow AIPW", family, n_fisher = NULL,
            gamma_sel = g, output_frt = TRUE
          )
          res_H1$res$est[1]
        })
        # bind_rows(
        #   tibble(type = "H0", t = T_H0),
        #   tibble(type = "H1", t = T_H1)
        # ) %>%
        #   ggplot(aes(t, fill = type)) +
        #   geom_histogram()
        critical_value <- quantile(abs(T_H0), probs = 1 - sig_level)
        power_hat_g <- mean(T_H1 > critical_value) # estimated power
        # output
        power_hat_g
      }) %>% map_dbl(~.)
    }
    # print
    walk2(gamma_grid, power_hat, function(g, power_hat_g) {
      cat(paste0("For gamma_sel = ", g, ", power = ", power_hat_g, "\n\n"))
    })
    # output
    out <- gamma_grid[which.max(power_hat)]
  }

  if (parallel) {
    future::plan(sequential)
  }

  return(out)
}
