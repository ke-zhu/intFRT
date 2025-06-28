#' Fisher Randomization Test and Conformal Selective Borrowing
#'
#' This function implements Fisher randomization tests (FRT) and Conformal
#' Selective Borrowing (CSB) in hybrid controlled trials. This function also
#' includes a series of external control borrowing estimators and their
#' inference results based on asymptotic normality for comparison (see the
#' `method` argument).
#'
#' @param Y Numeric vector of outcomes.
#' @param A Binary treatment indicator vector (1 for treatment, 0 for control).
#' @param S Binary indicator vector for source of data (1 for randomized trial
#'   data, 0 for external control data).
#' @param X Matrix of covariates.
#' @param method Character string specifying the method for estimation. Options
#'   include:
#'   * "No Borrow DiM"
#'   * "No Borrow AIPW"
#'   * "Borrow OM"
#'   * "Borrow IPW"
#'   * "Borrow staIPW"
#'   * "Borrow CW"
#'   * "Borrow AIPW"
#'   * "Borrow ACW"
#'   * "Conformal Selective Borrow AIPW" (the proposed & defalut)
#'   * "Conformal Selective Borrow ACW" (the proposed)
#' @param family A description of the outcome type. Default is "gaussian" for
#'   continuous outcomes, but can also be "binomial" for binary outcomes.
#' @param n_fisher Integer specifying the number of randomizations for the
#'   Fisher randomization test (FRT). Default is `NULL`, meaning FRT is skipped
#'   to avoid long computation times without the user's permission. For reliable
#'   FRT, set `n_fisher > 5000`.
#' @param n_boot Integer specifying the number of bootstrap samples for
#'   computing standard errors (SE), only applicable for IPW/staIPW/CW/OM.
#'   Default is `NULL`, meaning no bootstrap sampling is performed. To obtain
#'   reliable SE, set `n_boot > 5000`.
#' @param gamma_sel Numeric value controlling the conformal selective borrowing
#'   threshold. Default is `0.6`. It can also be adaptively determined using
#'   [compute_ada_gamma()].
#' @param cf Character string specifying the type of conformal inference to use.
#'   Default is "cv+", with other options being "split", "full", and
#'   "jackknife+".
#' @param cf_score Character string specifying the scoring function for
#'   conformal inference. Default is "CQR" for conformal quantile regression.
#'   Another option is "AR" for absolute residuals.
#' @param cf_model Character string specifying the model to use for outcome
#'   prediction in conformal inference. Default is "glm" for generalized linear
#'   model. Other options are "rf" for random forest and "ral" for relaxed
#'   adaptive lasso.
#' @param split_train Numeric value (between 0 and 1) specifying the proportion
#'   of data to use for training when `cf = "split"`. Default is 0.75.
#' @param cv_fold Integer specifying the number of folds for cross-validation
#'   when `cf = "cv+"`. Default is 10.
#' @param outcome_model Character string specifying the model to use for outcome
#'   imputation in OM/AIPW/ACW methods. Default is "glm" for generalized linear
#'   model. Other options are "rf" for random forest and "ral" for relaxed
#'   adaptive lasso.
#' @param sampling_model Character string specifying the model to use for
#'   sampling model. Default is "glm" for generalized linear model.
#'   Other options are "rf" for random forest and "ral" for relaxed adaptive
#'   lasso.
#' @param X_cw_ind Integer vector of column indices in `X` to balance using
#' calibration weighting. Default is `NULL`: uses all columns, or applies
#' data-adaptive selection if `sampling_model = "ral"`.
#' @param max_r Numeric value limiting the ratio of residual variances between
#'   trial and external control groups. Default is `Inf` (no restriction).
#' @param sig_level Numeric significance level for hypothesis tests. Default is
#'   0.05.
#' @param small_n_adj Logical value indicating whether to adjust variance
#'   estimation for small sample sizes in asymptotic normality inference. If
#'   `TRUE`, the sandwich variance estimation will multiply by `n / (n - d)`,
#'   where `d` is the number of estimated parameters and `n` is the number of
#'   units used to compute the sandwich variance estimator. Default is `TRUE`.
#' @param parallel Logical value indicating whether to use parallel computing.
#'   Default is `FALSE`.
#' @param n_cores Integer specifying the number of cores to use for parallel
#'   computing. Default is the number of available physical cores on the
#'   machine, as determined by `future::availableCores(logical = FALSE)`.
#' @param output_frt Logical value indicating whether to output all test
#'   statistic values from the Fisher randomization test. Default is `FALSE`.
#'
#' @return A list containing the following elements:
#'    * `res`: A tibble containing estimation and inference results:
#'      - `est`: ATE estimate.
#'      - `se`: Standard error.
#'      - `ci_l`, `ci_u`: Confidence interval.
#'      - `p_value` (first row): P-value based on asymptotic normality.
#'      - `p_value` (second row): P-value based on FRT (If `n_fisher != NULL`).
#'      - `n_sel`: The number of borrowed external controls.
#'      - `ess_sel`: The effective sample size of borrowed external controls.
#'      - `runtime`: The computation time.
#'    * `id_sel`: The row numbers of selected external controls.
#'    * `dat_info`: A tibble containing main information about the original data:
#'      - `n_rt`: The number of units under randomized treatment.
#'      - `n_rc`: The number of units under randomized control.
#'      - `n_rct`: The total number of units in the randomized controlled trial.
#'      - `n_et`: The number of all external controls.
#'      - `id_et`: The row numbers of all external controls.
#'    * `gamma_sel`: Conformal selection threshold.
#'    * `out`: A tibble containing raw output.
#'    * `out_frt`: A tibble containing output of FRT, if `output_frt` is `TRUE`.
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
#' # Perform Fisher Randomization Test with No Borrowing
#' result_nb <- ec_borrow(
#'   Y = Y,
#'   A = A,
#'   S = S,
#'   X = X,
#'   method = "No Borrow AIPW",
#'   family = "gaussian",
#'   # FRT with 10 randomizations for illustration;
#'   # recommend n_fisher = 5000 for more powerful results
#'   n_fisher = 10
#' )
#'
#' # View results
#' print(result_nb$res)
#'
#' # View IDs of borrowed external controls (None for No Borrowing)
#' result_nb$id_sel
#'
#' # Compute adaptive gamma (with a small n_rep_gamma for illustration)
#' ada_g <- compute_ada_gamma(
#'   Y, A, S, X,
#'   # Tuning with 20 replications for illustration purposes
#'   # Recommend setting `n_rep_gamma = 100` or higher with parallel computing
#'   n_rep_gamma = 20
#' )
#'
#' # Perform Fisher Randomization Test with Conformal Selective Borrowing
#' result_csb <- ec_borrow(
#'   Y = Y,
#'   A = A,
#'   S = S,
#'   X = X,
#'   method = "Conformal Selective Borrow AIPW",
#'   family = "gaussian",
#'   gamma_sel = ada_g,
#'   # FRT with 10 randomizations for illustration
#'   # recommend n_fisher = 5000+ for more powerful results
#'   n_fisher = 10
#' )
#'
#' # View results
#' print(result_csb$res)
#'
#' # View IDs of borrowed external controls
#' result_csb$id_sel
#'
#' # Mark borrowed external controls for plotting
#' sel <- rep(0, length(Y))
#' sel[result_csb$id_sel] <- 1
#'
#' # Visualize borrowed external controls with ggplot
#' library(dplyr)
#' library(ggplot2)
#' library(forcats)
#' library(quantreg)
#'
#' # Fit sampling score for multiple covariates
#' `Sampling Score` <- glm(S ~ X) %>% predict(type = "response")
#'
#' # Prepare data for plotting
#' dat_plot <- tibble(Y, A, S, X, `Sampling Score`, sel) %>%
#'   mutate(
#'     Type = case_when(
#'       A == 1 & S == 1 ~ "RCT treated",
#'       A == 0 & S == 1 ~ "RCT control",
#'       A == 0 & S == 0 & sel == 1 ~ "External control (selected)",
#'       A == 0 & S == 0 & sel == 0 ~ "External control (unselected)"
#'     ) %>%
#'       as_factor %>%
#'       fct_relevel(
#'         "RCT treated", "RCT control",
#'         "External control (selected)", "External control (unselected)"
#'       )
#'   ) %>%
#'   filter(Type != "RCT treated")
#'
#' # Fit quantile regressions to visualize the main range of RCT controls
#' fit975 <- rq(Y ~ `Sampling Score`, tau = 0.975, data = dat_plot,
#'              subset = dat_plot$Type == "RCT control")
#' dat_plot$pred975 <- predict(fit975, newdata = dat_plot)
#'
#' fit025 <- rq(Y ~ `Sampling Score`, tau = 0.025, data = dat_plot,
#'              subset = dat_plot$Type == "RCT control")
#' dat_plot$pred025 <- predict(fit025, newdata = dat_plot)
#'
#' # Plot results
#' dat_plot %>%
#'   ggplot() +
#'   geom_ribbon(aes(`Sampling Score`, ymin = pred025, ymax = pred975),
#'               fill = "grey80", alpha = 0.5) +
#'   geom_point(aes(`Sampling Score`, Y, color = Type)) +
#'   scale_color_manual(values = c("#5A5A5A", "#00ADFA", "#F8766D"))
#' }
#'
#' @import tibble
#' @import dplyr
#' @import purrr
#' @importFrom stats as.formula fisher.test gaussian glm pnorm predict qnorm resid sd t.test var
#' @importFrom utils head tail
#' @export
ec_borrow <- function(
    Y, A, S, X,
    method = "Conformal Selective Borrow AIPW",
    family = "gaussian",
    n_fisher = NULL,
    # IPW/staIPW/OM/CW
    n_boot = NULL,
    # CSB
    gamma_sel = 0.6,
    cf = "cv+",
    cf_score = "CQR",
    cf_model = "glm",
    split_train = 0.75,
    cv_fold = 10,
    # AIPW
    outcome_model = "glm",
    sampling_model = "glm",
    X_cw_ind = NULL,
    max_r = Inf,
    # testing
    sig_level = 0.05,
    small_n_adj = TRUE,
    # computing & output
    parallel = FALSE,
    n_cores = NULL,
    output_frt = FALSE
) {
  if (parallel) {
    rlang::check_installed("future", reason = "to set execution plan for `furrr`")
    rlang::check_installed("furrr", reason = "to use `furrr::future_map()` for parallel execution")
    if (is.null(n_cores)) {
      n_cores <- future::availableCores(logical = FALSE)
    }
    future::plan(multisession, workers = n_cores)
  }

  if ("ral" %in% c(outcome_model, sampling_model))
    rlang::check_installed("glmnet", reason = "to use `glmnet::cv.glmnet()` for 'ral'")

  if ("rf" %in% c(outcome_model, sampling_model))
    rlang::check_installed("randomForest", reason = "to use `randomForest()` for 'rf'")

  if (cf_model == "ral")
    rlang::check_installed("rqPen", reason = "to use `rq.pen.cv()` for cf_model = 'ral'")

  if (cf_model == "rf")
    rlang::check_installed("grf", reason = "to use `quantile_forest()` for cf_model = 'rf'")

  if (cf_score == "NN")
    rlang::check_installed("RANN", reason = "to use `nn2()` for cf_score = 'NN'")

  dat_origin <- tibble(Y, A, S, X)
  rm(Y, A, S, X)

  # No Borrowing
  #   No Borrow DiM
  #   No Borrow AIPW

  # Borrowing
  #   Borrow OM
  #   Borrow IPW, Borrow staIPW, Borrow CW
  #   Borrow AIPW, Borrow ACW

  # Selective Borrowing
  #   AdaLasso Selective Borrow ACW (Chenyin)
  #   Conformal Selective Borrow AIPW (proposed)
  #   Conformal Selective Borrow ACW (proposed)

  if (identical(method, "No Borrow DiM")) {
    est_fun <- function(dat) {
      fit <- t.test(
        x = dat %>% filter(A == 1, S == 1) %>% pull(Y),
        y = dat %>% filter(A == 0, S == 1) %>% pull(Y)
      )
      tibble(
        est = -(fit$estimate %>% diff %>% unname),
        se = fit$stderr,
        ci_l = fit$conf.int[1],
        ci_u = fit$conf.int[2],
        p_value = fit$p.value,
        # borrow no EC
        ess_sel = 0,
        id_sel = list(NULL)
      )
    }
    gamma_sel <- 1
  }

  if (identical(method, "No Borrow AIPW")) {
    est_fun <- function(dat) {
      rct_aipw(dat, family, outcome_model, small_n_adj) %>%
        # borrow no EC
        mutate(
          ess_sel = 0,
          id_sel = list(NULL)
        )
    }
    gamma_sel <- 1
  }

  if (identical(method, "Borrow Naive")) {
    est_fun <- function(dat) {
      n_ec <- dat %>% filter(A == 0, S == 0) %>% nrow
      dat_naive <- dat %>% mutate(S = 1)
      rct_aipw(dat_naive, family, outcome_model, small_n_adj) %>%
        # borrow no EC
        mutate(
          ess_sel = n_ec,
          id_sel = list(which(dat$S == 0))
        )
    }
    gamma_sel <- 0
  }

  if (identical(method, "Borrow OM")) {
    est_fun <- function(dat) {
      n_ec <- dat %>% filter(A == 0, S == 0) %>% nrow

      m10 <- fit_outcome_model(dat, family, outcome_model)
      m1 <- m10$m1
      m0 <- m10$m0
      d <- dat %>%
        mutate(d_i = m1 - m0) %>%
        filter(S == 1) %>%
        pull(d_i)
      tibble(
        est = mean(d),
        # borrow all ECs
        ess_sel = n_ec,
        id_sel = list(which(dat$S == 0))
      )
    }
    gamma_sel <- 0
  }

  if (identical(method, "Borrow IPW")) {
    est_fun <- function(dat) {
      n_rt <- dat %>% filter(A == 1, S == 1) %>% nrow
      n_rc <- dat %>% filter(A == 0, S == 1) %>% nrow
      n_rct <- dat %>% filter(S == 1) %>% nrow
      n_all <- dat %>% nrow

      # treatment group
      # use true propensity score
      pA <- n_rt / n_rct
      w1 <- with(
        dat,
        S * A / pA
      )
      d1 <- with(
        dat,
        (n_all / n_rct) * w1 * Y
      )

      # control group
      # compute r
      if (family == "gaussian") {
        r1 <- glm(Y ~ X, family = family, dat %>% filter(A == 0, S == 1)) %>%
          resid(type = "response") %>% var
        r0 <- glm(Y ~ X, family = family, dat %>% filter(S == 0)) %>%
          resid(type = "response") %>% var
        r <- min(r1 / r0, max_r)
      } else if (family == "binomial") {
        # for binary outcome, under exchangeablity assumption, r=1 (Li et al., 2023)
        r <- 1
      }
      # compute qhat
      pS <- glm(S ~ X, family = "binomial", dat) %>% predict(dat, "response")
      qhat <- pS / (1 - pS)
      w0 <- with(
        dat,
        qhat * (S * (1 - A) + (1 - S) * r) / (qhat * (1 - pA) + r)
      )
      d0 <- with(
        dat,
        (n_all / n_rct) * w0 * Y
      )

      # compute est
      d <- d1 - d0
      # output
      tibble(
        est = mean(d),
        # borrow all ECs
        ess_sel = max(0, ESS(w0) - n_rc),
        id_sel = list(which(dat$S == 0))
      )
    }
    gamma_sel <- 0
  }

  if (identical(method, "Borrow staIPW")) {
    est_fun <- function(dat) {
      n_rt <- dat %>% filter(A == 1, S == 1) %>% nrow
      n_rc <- dat %>% filter(A == 0, S == 1) %>% nrow
      n_rct <- dat %>% filter(S == 1) %>% nrow
      n_all <- dat %>% nrow

      # treatment group
      # use true propensity score
      pA <- n_rt / n_rct
      w1 <- with(
        dat,
        S * A / pA
      )
      d1 <- with(
        dat,
        (n_all / n_rct) * w1 * Y
      )

      # control group
      # compute r
      if (family == "gaussian") {
        r1 <- glm(Y ~ X, family = family, dat %>% filter(A == 0, S == 1)) %>%
          resid(type = "response") %>% var
        r0 <- glm(Y ~ X, family = family, dat %>% filter(S == 0)) %>%
          resid(type = "response") %>% var
        r <- min(r1 / r0, max_r)
      } else if (family == "binomial") {
        # for binary outcome, under exchangeablity assumption, r=1 (Li et al., 2023)
        r <- 1
      }
      # compute qhat
      pS <- glm(S ~ X, family = "binomial", dat) %>% predict(dat, "response")
      qhat <- pS / (1 - pS)
      w0init <- with(
        dat,
        qhat * (S * (1 - A) + (1 - S) * r) / (qhat * (1 - pA) + r)
      )
      w0 <- w0init / sum(w0init) * n_rct
      d0 <- with(
        dat,
        (n_all / n_rct) * w0 * Y
      )

      # compute est
      d <- d1 - d0
      # output
      tibble(
        est = mean(d),
        # borrow all ECs
        ess_sel = max(0, ESS(w0) - n_rc),
        id_sel = list(which(dat$S == 0))
      )
    }
    gamma_sel <- 0
  }

  if (identical(method, "Borrow CW")) {
    est_fun <- function(dat) {
      n_rt <- dat %>% filter(A == 1, S == 1) %>% nrow
      n_rc <- dat %>% filter(A == 0, S == 1) %>% nrow
      n_rct <- dat %>% filter(S == 1) %>% nrow
      n_all <- dat %>% nrow

      # treatment group
      # use true propensity score
      pA <- n_rt / n_rct
      w1 <- with(
        dat,
        S * A / pA
      )
      d1 <- with(
        dat,
        (n_all / n_rct) * w1 * Y
      )

      # control group
      # compute r
      if (family == "gaussian") {
        r1 <- glm(Y ~ X, family = family, dat %>% filter(A == 0, S == 1)) %>%
          resid(type = "response") %>% var
        r0 <- glm(Y ~ X, family = family, dat %>% filter(S == 0)) %>%
          resid(type = "response") %>% var
        r <- min(r1 / r0, max_r)
      } else if (family == "binomial") {
        # for binary outcome, under exchangeablity assumption, r=1 (Li et al., 2023)
        r <- 1
      }
      # compute qhat
      if (is.null(X_cw_ind)) { # no user-specified X_cw
        if (sampling_model == "ral") { # data-adaptive X_cw by outcome-adaptive lasso
          if (is.null(X_cw_ind)) {
            X_cw_ind_rc <- pred_model(
              dat %>% filter(A == 0, S == 1),
              dat %>% filter(A == 0, S == 1),
              family, base_model = "ral", var_sel = TRUE
            )
            X_cw_ind_ec <- pred_model(
              dat %>% filter(A == 0, S == 0),
              dat %>% filter(A == 0, S == 0),
              family, base_model = "ral", var_sel = TRUE
            )
            X_cw_ind <- union(X_cw_ind_rc, X_cw_ind_ec)
            X_cw <- dat$X[, X_cw_ind, drop = FALSE]
          }
        } else { # use all X
          X_cw <- dat$X
        }
      } else { # user-specified X_cw
        X_cw <- dat$X[, X_cw_ind, drop = FALSE]
      }
      # calibration weighting
      qhat <- compute_cw(dat$S, X_cw)

      # compute weight
      w0init <- with(
        dat,
        qhat * (S * (1 - A) + (1 - S) * r) / (qhat * (1 - pA) + r)
      )
      w0 <- w0init / sum(w0init) * n_rct
      d0 <- with(
        dat,
        (n_all / n_rct) * w0 * Y
      )

      # compute est
      d <- d1 - d0
      # output
      tibble(
        est = mean(d),
        # borrow all ECs
        ess_sel = max(0, ESS(w0) - n_rc),
        id_sel = list(which(dat$S == 0))
      )
    }
    gamma_sel <- 0
  }

  if (identical(method, "Borrow AIPW")) {
    est_fun <- function(dat) {
      rct_ec_aipw_acw(dat, family, outcome_model, max_r, small_n_adj,
                      sampling_model) %>%
        # borrow all ECs
        mutate(id_sel = list(which(dat$S == 0)))
    }
    gamma_sel <- 0
  }

  if (identical(method, "Borrow ACW")) {
    est_fun <- function(dat) {
      rct_ec_aipw_acw(dat, family, outcome_model, max_r, small_n_adj,
                      sampling_model, cw = TRUE, X_cw_ind = X_cw_ind) %>%
        # borrow all ECs
        mutate(id_sel = list(which(dat$S == 0)))
    }
    gamma_sel <- 0
  }

  if (identical(method, "AdaLasso Selective Borrow ACW")) { # Chenyin's method
    rlang::check_installed("SelectiveIntegrative", reason = "to use `srEC()`")
    est_fun <- function(dat) {
      fit <- with(dat,
                  SelectiveIntegrative::srEC(
                    data_rt = list(X = X[S == 1,], Y = Y[S == 1], A = A[S == 1]),
                    data_ec = list(list(X = X[S == 0,], Y = Y[S == 0], A = A[S == 0])),
                    method = "glm"
                  )
      )
      tibble(
        est = as.vector(fit$est$ACW.final),
        se = as.vector(fit$sd$ACW.final) / sqrt(fit$n_c),
        ess_sel = NA,
        id_sel = list(which(dat$S == 0)[fit$subset.idx])
      )
    }
    gamma_sel <- NA
  }

  if (identical(method, c("Conformal Selective Borrow AIPW"))) { # proposed method
    # est_fun
    est_fun <- function(dat) {
      x <- dat %>% filter(A == 0) %>% pull(X)
      y <- dat %>% filter(A == 0) %>% pull(Y)
      s <- dat %>% filter(A == 0) %>% pull(S)
      p_cf <- conformal_p(
        x, y, s, family,
        cf, cf_score, cf_model, split_train, cv_fold, sig_level
      )
      # biased or unbiased
      bias_ec <- ifelse(p_cf > gamma_sel, 0, 1)
      # estimation
      if (sum(bias_ec == 0) < 5) {
        # if n_sel < 5, do not borrow anyone
        rct_aipw(dat, family, outcome_model, small_n_adj) %>%
          mutate(
            ess_sel = 0,
            id_sel = list(NULL)
          )
      } else {
        # if n_sel >= 5, borrow them
        bias <- rep(0, nrow(dat))
        bias[dat$S == 0] <- bias_ec
        dat_sel <- dat %>% filter(bias == 0)
        rct_ec_aipw_acw(dat_sel, family, outcome_model, max_r, small_n_adj,
                        sampling_model) %>%
          mutate(id_sel = list(which(dat$S == 0 & bias == 0)))
      }
    }
  }

  if (identical(method, c("Conformal Selective Borrow ACW"))) { # proposed method
    # est_fun
    est_fun <- function(dat) {
      x <- dat %>% filter(A == 0) %>% pull(X)
      y <- dat %>% filter(A == 0) %>% pull(Y)
      s <- dat %>% filter(A == 0) %>% pull(S)
      p_cf <- conformal_p(
        x, y, s, family,
        cf, cf_score, cf_model, split_train, cv_fold, sig_level
      )
      # biased or unbiased
      bias_ec <- ifelse(p_cf > gamma_sel, 0, 1)
      # estimation
      if (sum(bias_ec == 0) < 5) {
        # if n_sel < 5, do not borrow anyone
        rct_aipw(dat, family, outcome_model, small_n_adj) %>%
          mutate(
            ess_sel = 0,
            id_sel = list(NULL)
          )
      } else {
        # if n_sel >= 5, borrow them
        bias <- rep(0, nrow(dat))
        bias[dat$S == 0] <- bias_ec
        dat_sel <- dat %>% filter(bias == 0)
        rct_ec_aipw_acw(dat_sel, family, outcome_model, max_r, small_n_adj,
                        sampling_model, cw = TRUE, X_cw_ind = X_cw_ind) %>%
          mutate(id_sel = list(which(dat$S == 0 & bias == 0)))
      }
    }
  }

  # 1 Estimation
  # record run time
  runtime <- system.time(
    # raw output
    out <- est_fun(dat_origin)
  )[3] %>% unname

  # 2 Inference
  if (method %in% c("Borrow OM", "Borrow IPW", "Borrow staIPW", "Borrow CW")) {
    if (!is.null(n_boot)) {
      # compute bootstrap SE
      runtime_boot <- system.time(
        # bootstrap
        if (parallel) {
          cat(paste0("parallel computing enabled with ", n_cores,
                     " cores for bootstrap SE of ", method, "\n\n"))
          out_boot <- furrr::future_map(1:n_boot, function(i) {
            dat_boot <- dat_origin %>%
              group_by(A, S) %>%
              slice_sample(prop = 1, replace = TRUE) %>%
              ungroup()
            tryCatch({
              est_fun(dat_boot)$est
            }, error = function(e) {
              NA
            })
          }, .options = furrr::furrr_options(seed = TRUE)) %>%
            map_dbl(~.)
        } else {
          out_boot <- map_dbl(1:n_boot, ~ {
            dat_boot <- dat_origin %>%
              group_by(A, S) %>%
              slice_sample(prop = 1, replace = TRUE) %>%
              ungroup()
            tryCatch({
              est_fun(dat_boot)$est
            }, error = function(e) {
              NA
            })
          })
        }
      )[3] %>% unname
      out$se <- sd(out_boot, na.rm = TRUE)
      # warning for NA
      if (any(is.na(out_boot))) {
        warning(paste0("There are ", sum(is.na(out_boot)),
                       " NA in bootstrap for ", method))
      }
      # organize
      res <- tibble(
        method,
        est = out$est,
        se = out$se,
        ci_l = est - qnorm(1 - sig_level / 2) * se,
        ci_u = est + qnorm(1 - sig_level / 2) * se,
        p_value = (1 - pnorm(abs(est / se))) * 2,
        n_sel = map_dbl(out$id_sel, length),
        ess_sel = out$ess_sel,
        runtime = runtime + runtime_boot
      )
    } else {
      # organize
      res <- tibble(
        method,
        est = out$est,
        se = NA,
        ci_l = NA,
        ci_u = NA,
        p_value = NA,
        n_sel = map_dbl(out$id_sel, length),
        ess_sel = out$ess_sel,
        runtime = runtime
      )
    }
  } else {
    # for No Borrow AIPW; Borrow AIPW; Borrow ACW; AdaLasso Selective Borrow ACW; Conformal Selective Borrow AIPW
    # organize
    res <- tibble(
      method,
      est = out$est,
      se = out$se,
      ci_l = est - qnorm(1 - sig_level / 2) * se,
      ci_u = est + qnorm(1 - sig_level / 2) * se,
      p_value = (1 - pnorm(abs(est / se))) * 2,
      n_sel = map_dbl(out$id_sel, length),
      ess_sel = out$ess_sel,
      runtime
    )
  }

  # 3 Fisher randomization test
  if (!is.null(n_fisher)) {
    if (identical(method, "No Borrow DiM") & family == "binomial") {
      # special case of Fisher's exact test
      x1 <- dat_origin %>% filter(S == 1, A == 0, Y == 1) %>% nrow()
      n1 <- dat_origin %>% filter(S == 1, A == 0) %>% nrow()
      x2 <- dat_origin %>% filter(S == 1, A == 1, Y == 1) %>% nrow()
      n2 <- dat_origin %>% filter(S == 1, A == 1) %>% nrow()
      runtime_frt <- system.time(
        fit <- matrix(
          c(x2, x1, n2 - x2, n1 - x1), 2, 2,
          dimnames = list(c("t", "c"), c("Event", "No Event"))
        ) %>% fisher.test()
      )[3] %>% unname
      res_frt <- tibble(
        method = paste0(method, "+FRT"),
        est = NA,
        se = NA,
        ci_l = NA,
        ci_u = NA,
        p_value = fit$p.value,
        n_sel = NA,
        ess_sel = NA,
        runtime = runtime_frt
      )
    } else {
      runtime_frt <- system.time(
        if (parallel) {
          cat(paste0("parallel computing enabled with ", n_cores,
                     " cores for ", method, "+FRT\n\n"))
          # randomization
          out_frt <- furrr::future_map(1:n_fisher, function(i) {
            dat_rand <- dat_origin %>%
              mutate(A = {A[S == 1] <- sample(A[S == 1]); A})
            tryCatch({
              est_fun(dat_rand)
            }, error = function(e) {
              NULL
            })
          }, .options = furrr::furrr_options(seed = TRUE)) %>%
            map_dfr(~.)
        } else {
          # randomization
          out_frt <- map_dfr(1:n_fisher, ~{
            dat_rand <- dat_origin %>%
              mutate(A = {A[S == 1] <- sample(A[S == 1]); A})
            tryCatch({
              est_fun(dat_rand)
            }, error = function(e) {
              NULL
            })
          })
        }
      )[3] %>% unname
      # warning for NA
      if (nrow(out_frt) < n_fisher) {
        warning(paste0("nrow(out_frt) is ", nrow(out_frt), " for ",
                       method, "+FRT"))
      }
      if (any(is.na(out_frt$est))) {
        warning(paste0("There are ", sum(is.na(out_frt$est)), " NA in ",
                       method, "+FRT"))
      }
      # organize
      res_frt <- tibble(
        method = paste0(method, "+FRT"),
        est = NA,
        se = NA,
        ci_l = NA,
        ci_u = NA,
        p_value = mean(c(abs(out_frt$est) >= abs(out$est), 1), na.rm = T),
        n_sel = NA,
        ess_sel = NA,
        runtime = runtime_frt
      )
    }
    res <- rbind(res, res_frt)
  } else {
    out_frt <- NULL
  }
  dat_info <- tibble(
    n_rt = dat_origin %>% filter(A == 1, S == 1) %>% nrow,
    n_rc = dat_origin %>% filter(A == 0, S == 1) %>% nrow,
    n_rct = n_rt + n_rc,
    n_ec = dat_origin %>% filter(A == 0, S == 0) %>% nrow,
    id_ec = list(which(dat_origin$S == 0))
  )

  if (parallel) {
    future::plan(sequential)
  }

  if (output_frt) {
    lst(res, id_sel = out$id_sel[[1]], dat_info, gamma_sel, out, out_frt)
  } else {
    lst(res, id_sel = out$id_sel[[1]], dat_info, gamma_sel, out)
  }
}
