
# ATE estimator -----------------------------------------------------------

# fit_outcome_model <- function(dat, family, outcome_model) {
#   if (outcome_model == "glm") {
#     m1 <- glm(Y ~ X, family = family, dat %>% filter(A == 1)) %>%
#       predict(dat, "response")
#     m0 <- glm(Y ~ X, family = family, dat %>% filter(A == 0)) %>%
#       predict(dat, "response")
#   } else if (outcome_model == "rf") {
#     rlang::check_installed("randomForest", reason = "to use `randomForest()`")
#     if (family == "gaussian") {
#       m1 <- randomForest::randomForest(dat$X[dat$A == 1, , drop = FALSE], dat$Y[dat$A == 1]) %>%
#         predict(dat$X)
#       m0 <- randomForest::randomForest(dat$X[dat$A == 0, , drop = FALSE], dat$Y[dat$A == 0]) %>%
#         predict(dat$X)
#     } else if (family == "binomial") {
#       m1 <- randomForest::randomForest(
#         dat$X[dat$A == 1, , drop = FALSE], as.factor(dat$Y[dat$A == 1])
#       ) %>%
#         predict(dat$X, type = "prob") %>%
#         {.[, 2]}
#       m0 <- randomForest::randomForest(
#         dat$X[dat$A == 0, , drop = FALSE], as.factor(dat$Y[dat$A == 0])
#       ) %>%
#         predict(dat$X, type = "prob") %>%
#         {.[, 2]}
#     }
#   } else if (outcome_model == "ral") {
#     # Fit relaxed adaptive lasso for treated group
#     fit_ridge1 <- glmnet::cv.glmnet(
#       X[A == 1, ], Y[A == 1],
#       family = family, alpha = 0
#     )
#     fit_lasso1 <- glmnet::cv.glmnet(
#       X[A == 1, ], Y[A == 1],
#       family = family,
#       relax = TRUE, gamma = 0, # refitting
#       penalty.factor = abs(coef(fit_ridge1)[-1]) # adaptive weights
#     )
#     m1 <- predict(
#       fit_lasso1, newx = X,
#       s = "lambda.min", type = "response"
#     ) %>% as.vector()
#
#     # Fit relaxed adaptive lasso for control group
#     fit_ridge0 <- glmnet::cv.glmnet(
#       X[A == 0, ], Y[A == 0],
#       family = family, alpha = 0
#     )
#     fit_lasso0 <- glmnet::cv.glmnet(
#       X[A == 0, ], Y[A == 0],
#       family = family,
#       relax = TRUE, gamma = 0, # refitting
#       penalty.factor = abs(coef(fit_ridge0)[-1]) # adaptive weights
#     )
#     m0 <- predict(
#       fit_lasso0, newx = X,
#       s = "lambda.min", type = "response"
#     ) %>% as.vector()
#   }
#   lst(m1, m0)
# }

pred_model <- function(
    dat_train, dat_pred, family, base_model,
    var_sel = FALSE # if TRUE, output selected X rather than prediction
  ) {
  if (base_model != "ral") var_sel = FALSE # only work when base_model == "ral"

  if (base_model == "glm") {
    # GLM
    mu_pred <- glm(Y ~ X, family = family, dat_train) %>%
      predict(dat_pred, "response") %>%
      as.vector()
  } else if (base_model == "rf") {
    # random forest
    if (family == "gaussian") {
      mu_pred <- randomForest::randomForest(
        as.matrix(dat_train$X), dat_train$Y
      ) %>%
        predict(dat_pred$X) %>%
        as.vector()
    } else if (family == "binomial") {
      mu_pred <- randomForest::randomForest(
        as.matrix(dat_train$X), as.factor(dat_train$Y)
      ) %>%
        predict(as.matrix(dat_pred$X), type = "prob") %>%
        {.[, 2]} %>%
        as.vector()
    }
  } else if (base_model == "ral") {
    # relaxed adaptive lasso
    fit_ridge <- glmnet::cv.glmnet(
      as.matrix(dat_train$X), dat_train$Y,
      family = family,
      alpha = 0 # ridge
    )
    fit_ralasso <- glmnet::cv.glmnet(
      as.matrix(dat_train$X), dat_train$Y,
      family = family,
      relax = TRUE, gamma = 0, # refitting
      penalty.factor = abs(coef(fit_ridge)[-1]) # adaptive weights
    )
    mu_pred <- predict(
      fit_ralasso, newx = as.matrix(dat_pred$X),
      s = "lambda.min", type = "response"
    ) %>%
      as.vector()
    X_sel_ind <- which(coef(fit_ralasso, s = "lambda.min")[-1] != 0)
  }
  if (var_sel) {
    return(X_sel_ind)
  } else {
    return(mu_pred)
  }
}

compute_cw <- function(S, X) {
  n_rct <- sum(S == 1)
  n_ec <- sum(S == 0)

  X_try <- X
  success <- FALSE

  while (ncol(X_try) > 0) {
    # calibration weights (entropy balancing) for EC sample
    dat_df <- data.frame(X_try, S = S)
    tryCatch({
      W.out <- WeightIt::weightit(
        S ~ .,
        data = dat_df,
        estimand = "ATT",
        method = "ebal",
        include.obj = TRUE
      )
      success <- TRUE
    }, error = function(e) {
      # do nothing, just continue
    })
    if (success) {
      break
    } else {
      # remove last column and retry
      X_try <- X_try[, -ncol(X_try), drop = FALSE]
    }
  }

  if (success) {
    w <- W.out$weights / n_ec

    # reproduce calibration weights for all sample, since w[S == 1] is constant
    Xscale <- WeightIt:::.make_closer_to_1(X_try) %>% as.matrix()
    Xtarget <- sweep(Xscale, 2, colMeans(Xscale[S == 1, , drop = FALSE]))
    ww_init <- exp(-Xtarget %*% W.out$obj$`0`$par)
    ww <- ww_init / sum(ww_init[S == 0])
    # # check ww[S == 0] = w[S == 0]
    # max(abs(ww[S == 0] - w[S == 0]))
    # plot(ww[S == 0], w[S==0])
    # tibble(ww[S == 0], w[S==0])

    # calibration weights for RCT sample
    w[S == 1] <- ww[S == 1]
    qhat <- w * n_rct

    # if there is some strange weights, return constant weights
    if (any(!is.finite(qhat))) {
      pS <- n_rct / (n_rct + n_ec)
      qhat <- rep(pS / (1 - pS), n_rct + n_ec)
    }
  } else {
    # if all attempts fail, return constant weights
    pS <- n_rct / (n_rct + n_ec)
    qhat <- rep(pS / (1 - pS), n_rct + n_ec)
  }

  # check q_cw vs q_ipw
  # pS <- glm(S ~ X, family = "binomial", tibble(S,X)) %>%
  #   predict(tibble(S,X),"response")
  # tibble(q_cw = qhat, q_ipw = pS / (1 - pS)) %>%
  #   ggplot(aes(q_cw, q_ipw)) +
  #   geom_point()+
  #   geom_abline(slope = 1, intercept = 0, color = "red")
  # max(abs(qhat - pS / (1 - pS)))

  qhat
}

rct_aipw <- function(dat, family, outcome_model, small_n_adj) {
  n_rt <- dat %>% filter(A == 1, S == 1) %>% nrow
  n_rct <- dat %>% filter(S == 1) %>% nrow
  dat_rct <- dat %>% filter(S == 1)

  # treatment propensity score model (use true)
  pA <- n_rt / n_rct

  # outcome model
  # m10 <- fit_outcome_model(dat_rct, family, outcome_model)
  # m1 <- m10$m1
  # m0 <- m10$m0
  m1 <- pred_model(dat_rct %>% filter(A == 1), dat_rct, family, outcome_model)
  m0 <- pred_model(dat_rct %>% filter(A == 0), dat_rct, family, outcome_model)

  # EIF
  d <- with(
    dat_rct,
    m1 + A / pA * (Y - m1) - m0 - (1 - A) / (1 - pA) * (Y - m0)
  )

  # small sample adjustment for variance estimation
  if (small_n_adj) {
    dof <- max(n_rct - ncol(dat$X) * 2, 1)
    se <- sqrt(sum((d - mean(d))^2) / dof^2)
  } else {
    se <- sqrt(sum((d - mean(d))^2) / n_rct^2)
  }

  # output
  tibble(
    est = mean(d),
    se = se,
    d = list(d)
  )
}

rct_ec_aipw_acw <- function(dat, family, outcome_model, max_r, small_n_adj,
                            sampling_model, cw = FALSE, X_cw_ind = NULL) {
  n_rt <- dat %>% filter(A == 1, S == 1) %>% nrow
  n_rc <- dat %>% filter(A == 0, S == 1) %>% nrow
  n_rct <- dat %>% filter(S == 1) %>% nrow
  n_all <- dat %>% nrow

  # treatment propensity score model (use true)
  pA <- n_rt / n_rct

  # outcome model
  # m10 <- fit_outcome_model(dat, family, outcome_model)
  # m1 <- m10$m1
  # m0 <- m10$m0
  m1 <- pred_model(dat %>% filter(A == 1), dat, family, outcome_model)
  m0 <- pred_model(dat %>% filter(A == 0), dat, family, outcome_model)

  # conditional variance model (align with outcome model)
  if (family == "gaussian") {
    # r1 <- glm(Y ~ X, family = family, dat %>% filter(A == 0, S == 1)) %>%
    #   resid(type = "response") %>% var
    # r0 <- glm(Y ~ X, family = family, dat %>% filter(S == 0)) %>%
    #   resid(type = "response") %>% var

    m0_rc <- pred_model(
      dat %>% filter(A == 0, S == 1),
      dat %>% filter(A == 0, S == 1),
      family, outcome_model
    )
    y0_rc <- dat %>% filter(A == 0, S == 1) %>% pull(Y)
    r1 <- var(y0_rc - m0_rc)

    m0_ec <- pred_model(
      dat %>% filter(A == 0, S == 0),
      dat %>% filter(A == 0, S == 0),
      family, outcome_model
    )
    y0_ec <- dat %>% filter(A == 0, S == 0) %>% pull(Y)
    r0 <- var(y0_ec - m0_ec)

    r <- min(r1 / r0, max_r)
  } else if (family == "binomial") {
    # for binary outcome, under exchangeablity assumption, r=1 (Li et al., 2023)
    r <- 1
  }

  # weighting model
  if (cw) {
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
  } else {
    # sampling propensity score model
    # pS <- glm(S ~ X, family = "binomial", dat) %>% predict(dat, "response")
    pS <- pred_model(
      dat %>% mutate(Y = S),
      dat %>% mutate(Y = S),
      family = "binomial", base_model = sampling_model
    )
    qhat <- pS / (1 - pS)
  }
  w0init <- with(
    dat,
    qhat * (S * (1 - A) + (1 - S) * r) / (qhat * (1 - pA) + r)
  )
  w0 <- w0init / sum(w0init) * n_rct

  # treatment group EIF
  w1 <- with(
    dat,
    S * A / pA
  )
  d1 <- with(
    dat,
    (n_all / n_rct) * (S * m1 +  w1 * (Y - m1))
  )

  # control group EIF
  d0 <- with(
    dat,
    (n_all / n_rct) * (S * m0 +  w0 * (Y - m0))
  )

  # compute est
  d <- d1 - d0
  est <- mean(d)

  # compute se
  if (small_n_adj) {
    dof <- max(n_all - ncol(dat$X) * 3, 1)
    #se <- sqrt(sum((d - mean(d))^2) / dof^2)
  } else {
    #se <- sqrt(sum((d - mean(d))^2) / n_all^2)
    dof <- n_all
  }
  se_i <- d - dat$S / (n_rct / n_all) * mean(d)
  se <- sqrt(sum(se_i^2) / dof^2)

  # output
  tibble(est, se, ess_sel = max(0, ESS(w0) - n_rc), d = list(d))
}


# Conformal p-value -------------------------------------------------------


fit_cf_model_mean <- function(dat_train, dat_pred, family, cf_model) {
  if (cf_model == "glm") {
    fit <- glm(y ~ x, family = family, data = dat_train)
    predict(fit, newdata = dat_pred, type = "response")
  } else if (cf_model == "rf") {
    if (family == "gaussian") {
      fit <- randomForest::randomForest(x = dat_train$x, y = dat_train$y)
      predict(fit, newdata = dat_pred$x)
    } else if (family == "binomial") {
      fit <- randomForest::randomForest(x = dat_train$x, y = as.factor(dat_train$y))
      predict(fit, newdata = dat_pred$x, type = "prob")[,2]
    }
  }
}

fit_cf_model_mean_sd <- function(dat_train, dat_pred) { # family = "gaussian"
  # mean model
  fit1 <- glm(y ~ x, family = "gaussian", data = dat_train)
  ar_train <- abs(resid(fit1))
  # absolute residual model
  fit2 <- glm(ar_train ~ x, family = gaussian(link = "log"),
              data = dat_train %>% mutate(ar_train))
  yhat <- predict(fit1, newdata = dat_pred, type = "response")
  sighat <- predict(fit2, newdata = dat_pred, type = "response")
  lst(yhat, sighat)
}

fit_cf_model_quantile <- function(dat_train, dat_pred, a, cf_model) { # family = "gaussian"
  if (cf_model == "glm") {
    ci_mat <- tryCatch({
      fit_low <- quantreg::rq(y ~ x, tau = a/2, data = dat_train)
      fit_up  <- quantreg::rq(y ~ x, tau = 1 - a/2, data = dat_train)
      cbind(
        predict(fit_low, newdata = dat_pred),
        predict(fit_up,  newdata = dat_pred)
      )
    }, error = function(e) {
      warning(str_glue("quantreg::rq() failed. Using empirical quantiles instead."))
      y_q <- quantile(dat_train$y, probs = c(a/2, 1 - a/2), na.rm = TRUE)
      matrix(rep(y_q, each = nrow(dat_pred)), ncol = 2)
    })
    ci_mat
  } else if (cf_model == "rf") {
    ci_mat <- tryCatch({
      fit <- grf::quantile_forest(dat_train$x, dat_train$y, quantiles = c(a/2, 1 - a/2))
      predict(fit, dat_pred$x)$predictions
    }, error = function(e) {
      warning(str_glue("grf::quantile_forest() failed. Using empirical quantiles instead."))
      y_q <- quantile(dat_train$y, probs = c(a/2, 1 - a/2), na.rm = TRUE)
      matrix(rep(y_q, each = nrow(dat_pred$x)), ncol = 2)
    })
    ci_mat
  } else if (cf_model == "ral") {
    ci_mat <- tryCatch({
      x_train <- as.matrix(dat_train$x)
      y_train <- dat_train$y
      x_pred  <- as.matrix(dat_pred$x)
      taus <- c(a / 2, 1 - a / 2)

      # Step 1: Cross-validated adaptive lasso
      cv_fit <- rqPen::rq.pen.cv(
        x = x_train,
        y = y_train,
        tau = taus,
        penalty = "aLASSO",
        a = 1
      )

      # Step 2: Extract selected variables (non-zero coefficients)
      sel_low <- which(abs(coef(cv_fit)[-1, 1]) > 1e-6)
      sel_up  <- which(abs(coef(cv_fit)[-1, 2]) > 1e-6)

      if (length(sel_low) == 0 || length(sel_up) == 0) {
        stop("No variables selected by adaptive lasso rq.")
      }

      # Step 3: Refit unpenalized quantile regression
      relaxed_low <- quantreg::rq(
        y_train ~ .,
        data = data.frame(y_train, x_train[, sel_low, drop = FALSE]),
        tau = taus[1]
      )
      relaxed_up <- quantreg::rq(
        y_train ~ .,
        data = data.frame(y_train, x_train[, sel_up, drop = FALSE]),
        tau = taus[2]
      )

      X_low_pred <- x_pred[, sel_low, drop = FALSE]
      X_up_pred  <- x_pred[, sel_up,  drop = FALSE]
      colnames(X_low_pred) <- names(coef(relaxed_low))[-1]
      colnames(X_up_pred)  <- names(coef(relaxed_up))[-1]

      cbind(
        predict(relaxed_low, newdata = as.data.frame(X_low_pred)),
        predict(relaxed_up,  newdata = as.data.frame(X_up_pred))
      )
    }, error = function(e) {
      warning(str_glue("Relaxed adaptive lasso for rq failed. Using empirical quantiles instead."))
      y_q <- quantile(dat_train$y, probs = c(a / 2, 1 - a / 2), na.rm = TRUE)
      matrix(rep(y_q, each = nrow(dat_pred)), ncol = 2)
    })
    ci_mat
  }
}

conformal_p <- function(x, y, s, family, cf, cf_score, cf_model,
                        split_train, cv_fold, sig_level) {
  id_rc <- which(s == 1) # train & cal
  id_test <- which(s == 0) # test

  if (cf == "jackknife+") {
    cv_fold <- length(id_rc)
    cf <- "cv+"
  }

  if (cf == "full") {
    p_cf <- map_dbl(id_test, function(j) { # for every EC j
      # data augmentation
      id_train <- id_pred <- c(j, id_rc)
      dat_train <- dat_pred <- tibble(y = y[id_train], x = x[id_train, , drop = FALSE])
      # conformal score
      if (cf_score == "AR") {
        # prediction model
        yhat_pred <- fit_cf_model_mean(dat_train, dat_pred, family, cf_model)
        # compute score
        score_pred <- abs(yhat_pred - y[id_pred])
      } else if (cf_score == "local AR") {
        # prediction model
        hat_pred <- fit_cf_model_mean_sd(dat_train, dat_pred)
        # compute score
        score_pred <- abs(hat_pred$yhat - y[id_pred]) / hat_pred$sighat
      } else if (cf_score == "CQR") {
        # prediction model
        qhat_pred <- fit_cf_model_quantile(dat_train, dat_pred, sig_level, cf_model)
        # compute score
        score_pred <- pmax(qhat_pred[,1] - y[id_pred], y[id_pred] - qhat_pred[,2])
      } else if (cf_score == "NN") {
        # 1-nearest-neighbor
        d1 <- RANN::nn2(
          dat_train %>% filter(y == 1) %>% pull(x),
          k = 2
        )$nn.dists[,2]
        d0 <- RANN::nn2(
          dat_train %>% filter(y == 0) %>% pull(x),
          k = 2
        )$nn.dists[,2]
        # compute score
        score_pred <- rep(0, nrow(dat_pred))
        score_pred[dat_pred$y == 1] <- d1
        score_pred[dat_pred$y == 0] <- d0
      }
      score_test_j <- score_pred[1]
      score_cal <- score_pred[-1]
      # conformal p-values
      mean(c(score_cal >= score_test_j, 1))
    })
  } else if (cf == "split") {
    # data split
    n_train <- ceiling(split_train * length(id_rc))
    n_cal <- length(id_rc) - n_train
    split_id_rc <- split(
      id_rc,
      c(rep(1, n_cal), rep(2, n_train)) %>% sample
    )
    id_cal <- split_id_rc[[1]]
    id_train <- split_id_rc[[2]]
    dat_train <- tibble(y = y[id_train], x = x[id_train, , drop = FALSE])
    id_pred <- c(id_test, id_cal)
    dat_pred <- tibble(y = y[id_pred], x = x[id_pred, , drop = FALSE])
    # conformal score
    if (cf_score == "AR") {
      # prediction model
      yhat_pred <- fit_cf_model_mean(dat_train, dat_pred, family, cf_model)
      # compute score
      score_pred <- abs(yhat_pred - y[id_pred])
    } else if (cf_score == "local AR") {
      # prediction model
      hat_pred <- fit_cf_model_mean_sd(dat_train, dat_pred)
      # compute score
      score_pred <- abs(hat_pred$yhat - y[id_pred]) / hat_pred$sighat
    } else if (cf_score == "CQR") {
      # prediction model
      qhat_pred <- fit_cf_model_quantile(dat_train, dat_pred, sig_level, cf_model)
      # compute score
      score_pred <- pmax(qhat_pred[,1] - y[id_pred], y[id_pred] - qhat_pred[,2])
    } else if (cf_score == "NN") {
      # 1-nearest-neighbor
      d1 <- RANN::nn2(
        dat_train %>% filter(y == 1) %>% pull(x),
        dat_pred %>% filter(y == 1) %>% pull(x),
        k = 1
      )$nn.dists %>% as.vector()
      d0 <- RANN::nn2(
        dat_train %>% filter(y == 0) %>% pull(x),
        dat_pred %>% filter(y == 0) %>% pull(x),
        k = 1
      )$nn.dists %>% as.vector()
      # compute score
      score_pred <- rep(0, nrow(dat_pred))
      score_pred[dat_pred$y == 1] <- d1
      score_pred[dat_pred$y == 0] <- d0
    }
    score_test <- head(score_pred, length(id_test))
    score_cal <- tail(score_pred, length(id_cal))
    # conformal p-values
    p_cf <- map_dbl(score_test, function(score_test_j) {
      mean(c(score_cal >= score_test_j, 1))
    })
  } else if (cf == "cv+") {
    # data split
    split_id_rc <- split(
      id_rc,
      (rep(1:cv_fold, ceiling(length(id_rc) / cv_fold))[1:length(id_rc)]) %>% sample
    )
    # conformal score
    compare_all <- map(split_id_rc, function(id_cal) {
      id_train <- setdiff(id_rc, id_cal)
      dat_train <- tibble(y = y[id_train], x = x[id_train, , drop = FALSE])
      id_pred <- c(id_test, id_cal)
      dat_pred <- tibble(y = y[id_pred], x = x[id_pred, , drop = FALSE])
      if (cf_score == "AR") {
        # prediction model
        yhat_pred <- fit_cf_model_mean(dat_train, dat_pred, family, cf_model)
        # compute score
        score_pred <- abs(yhat_pred - y[id_pred])
      } else if (cf_score == "local AR") {
        # prediction model
        hat_pred <- fit_cf_model_mean_sd(dat_train, dat_pred)
        # compute score
        score_pred <- abs(hat_pred$yhat - y[id_pred]) / hat_pred$sighat
      } else if (cf_score == "CQR") {
        # prediction model
        qhat_pred <- fit_cf_model_quantile(dat_train, dat_pred, sig_level, cf_model)
        # compute score
        score_pred <- pmax(qhat_pred[,1] - y[id_pred], y[id_pred] - qhat_pred[,2])
      } else if (cf_score == "NN") {
        # 1-nearest-neighbor
        d1 <- RANN::nn2(
          dat_train %>% filter(y == 1) %>% pull(x),
          dat_pred %>% filter(y == 1) %>% pull(x),
          k = 1
        )$nn.dists %>% as.vector()
        d0 <- RANN::nn2(
          dat_train %>% filter(y == 0) %>% pull(x),
          dat_pred %>% filter(y == 0) %>% pull(x),
          k = 1
        )$nn.dists %>% as.vector()
        # compute score
        score_pred <- rep(0, nrow(dat_pred))
        score_pred[dat_pred$y == 1] <- d1
        score_pred[dat_pred$y == 0] <- d0
      }
      score_test_k <- head(score_pred, length(id_test))
      score_cal_k <- tail(score_pred, length(id_cal))
      # compare
      map_dbl(score_test_k, function(score_test_l) {
        sum(score_cal_k >= score_test_l)
      })
    }) %>% sapply(function(x) x)
    # conformal p-values
    p_cf <- map_dbl(1:length(id_test), function(l) {
      (sum(compare_all[l,]) + 1) / (length(id_rc) + 1)
    })
  }
  p_cf
}

ESS <- function (w) {
  sum(w)^2 / sum(w^2)
}

add_name <- function(l, new) {
  l$res <- l$res %>% mutate(method = paste0(method, " ", new))
  l
}
