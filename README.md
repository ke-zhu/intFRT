
<!-- README.md is generated from README.Rmd. Please edit that file -->

# intFRT

<!-- badges: start -->
<!-- badges: end -->

The goal of **intFRT** is to **int**egrate randomized controlled trials
(RCTs) with external controls (ECs) in hybrid controlled trials,
harnessing Fisher randomization tests (**FRT**) and Conformal Selective
Borrowing (CSB). It enhances the statistical efficiency of average
treatment effect (ATE) estimation and inference while ensuring valid
hypothesis testing. Key features include:

🎲 **Fisher Randomization Tests (FRT)**: Ensures strict control of the
Type I error rate in hybrid controlled trials.

🧲 **Conformal Selective Borrowing (CSB)**: Serves as both an ATE
estimator and a test statistic for FRT, selectively borrowing comparable
ECs to reduce hidden bias and enhance statistical power.

⚙️ **Adaptive Selection Threshold**: Dynamically determines the
selection threshold for CSB.

📊 **Flexible Estimator Options**: Offers No Borrowing and various EC
borrowing estimators (IPW, staIPW, CW, OM, AIPW, ACW) with inference
results based on asymptotic normality. These estimators can also be used
as test statistics in FRT.

## Installation

You can install the development version of intFRT from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ke-zhu/intFRT")
```

## Example

This example illustrates the use of the Fisher Randomization Test (FRT)
and Conformal Selective Borrowing methods in hybrid controlled trials.

### Simulate data for a hybrid controlled trial

``` r
set.seed(1)
n_rct <- 50  # Number of observations in the randomized controlled trial
n_ec <- 100  # Number of external controls
n <- n_rct + n_ec  # Total number of observations

# Covariates (2 covariates, uniformly distributed)
X <- matrix(runif(n * 2), n, 2)

# Treatment assignment (1 = treatment, 0 = control)
A <- c(rep(0:1, each = n_rct / 2), rep(0, n_ec))

# Data source indicator (1 = randomized trial, 0 = external control)
S <- c(rep(1, n_rct), rep(0, n_ec))

# Generate potential outcomes (continuous)
Y1 <- 1 + 2 * X[,1] + 0.1 * X[,2] + rnorm(n)
Y0 <- 2 * X[,1] + 0.1 * X[,2] + rnorm(n)

# Introduce bias in half of the external controls
id_biased_EC <- tail(which(S == 0), n = n_ec / 2)
Y0[id_biased_EC] <- Y0[id_biased_EC] - 10

# Observed outcome
Y <- A * Y1 + (1 - A) * Y0
```

### Adaptive selection threshold

``` r
library(intFRT)
ada_g <- compute_ada_gamma(Y, A, S, X)
#> For gamma_sel = 0, MSE = 3.02009406877217
#> Warning in est_g$est_one - est_grid[[id_nb]]$est_one: longer object length is
#> not a multiple of shorter object length
#> For gamma_sel = 0.1, MSE = 0.202123119515767
#> Warning in est_g$est_one - est_grid[[id_nb]]$est_one: longer object length is
#> not a multiple of shorter object length
#> For gamma_sel = 0.2, MSE = 0.188811976981897
#> Warning in est_g$est_one - est_grid[[id_nb]]$est_one: longer object length is
#> not a multiple of shorter object length
#> For gamma_sel = 0.3, MSE = 0.178675219443114
#> Warning in est_g$est_one - est_grid[[id_nb]]$est_one: longer object length is
#> not a multiple of shorter object length
#> For gamma_sel = 0.4, MSE = 0.144774075625022
#> Warning in est_g$est_one - est_grid[[id_nb]]$est_one: longer object length is
#> not a multiple of shorter object length
#> For gamma_sel = 0.5, MSE = 0.12768567184637
#> Warning in est_g$est_one - est_grid[[id_nb]]$est_one: longer object length is
#> not a multiple of shorter object length
#> For gamma_sel = 0.6, MSE = 0.110760514850767
#> Warning in est_g$est_one - est_grid[[id_nb]]$est_one: longer object length is
#> not a multiple of shorter object length
#> For gamma_sel = 0.7, MSE = 0.0999859347910325
#> Warning in est_g$est_one - est_grid[[id_nb]]$est_one: longer object length is
#> not a multiple of shorter object length
#> For gamma_sel = 0.8, MSE = 0.0659669326984849
#> Warning in est_g$est_one - est_grid[[id_nb]]$est_one: longer object length is
#> not a multiple of shorter object length
#> For gamma_sel = 0.9, MSE = 0.0465028607316254
#> 
#> For gamma_sel = 1, MSE = 0.076003315378369
ada_g
#> [1] 0.9
```

### Fisher Randomization Test with Conformal Selective Borrowing

``` r
result_csb <- ec_borrow(
  Y = Y,
  A = A,
  S = S,
  X = X,
  method = "Conformal Selective Borrow AIPW",
  family = "gaussian",
  gamma_sel = ada_g,
  # FRT with 10 randomizations for illustration
  n_fisher = 10
  # Recommend `n_fisher = 5000` or more
  # To perform only Conformal Selective Borrowing, set `n_fisher = NULL`
)
#> Warning in abs(out_frt$est) >= abs(out$est): longer object length is not a
#> multiple of shorter object length

# View results
print(result_csb$res, width = Inf)
#> # A tibble: 51 × 9
#>    method                            est    se  ci_l  ci_u p_value n_sel ess_sel
#>    <chr>                           <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl>   <dbl>
#>  1 Conformal Selective Borrow AIPW 0.793 0.300 0.206  1.38 0.00812     0       0
#>  2 Conformal Selective Borrow AIPW 0.793 0.300 0.206  1.38 0.00812     0       0
#>  3 Conformal Selective Borrow AIPW 0.793 0.300 0.206  1.38 0.00812     0       0
#>  4 Conformal Selective Borrow AIPW 0.793 0.300 0.206  1.38 0.00812     0       0
#>  5 Conformal Selective Borrow AIPW 0.793 0.300 0.206  1.38 0.00812     0       0
#>  6 Conformal Selective Borrow AIPW 0.793 0.300 0.206  1.38 0.00812     0       0
#>  7 Conformal Selective Borrow AIPW 0.793 0.300 0.206  1.38 0.00812     0       0
#>  8 Conformal Selective Borrow AIPW 0.793 0.300 0.206  1.38 0.00812     0       0
#>  9 Conformal Selective Borrow AIPW 0.793 0.300 0.206  1.38 0.00812     0       0
#> 10 Conformal Selective Borrow AIPW 0.793 0.300 0.206  1.38 0.00812     0       0
#>    runtime
#>      <dbl>
#>  1  0.0460
#>  2  0.0460
#>  3  0.0460
#>  4  0.0460
#>  5  0.0460
#>  6  0.0460
#>  7  0.0460
#>  8  0.0460
#>  9  0.0460
#> 10  0.0460
#> # ℹ 41 more rows
```

- `est`: ATE estimate.
- `se`, `ci_l`, `ci_u`, `p_value` (first row): Standard error,
  confidence interval, p-value based on asymptotic normality.
- `p_value` (second row): P-value based on FRT.
- `n_sel`: The number of borrowed external controls.
- `ess_sel`: The effective sample size of borrowed external controls.
- `runtime`: The computation time.

``` r
# View IDs of borrowed external controls
result_csb$id_sel
#> NULL
```

### Visualize borrowed external controls

``` r
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(ggplot2)
library(forcats)
library(quantreg)
#> Loading required package: SparseM

# Mark borrowed external controls for plotting
sel <- rep(0, length(Y))
sel[result_csb$id_sel] <- 1

# Fit sampling score for multiple covariates 
`Sampling Score` <- glm(S ~ X) %>% predict(type = "response")

# Prepare data for plotting
dat_plot <- tibble(Y, A, S, X, `Sampling Score`, sel) %>% 
  mutate(
    Type = case_when(
      A == 1 & S == 1 ~ "RCT treated",
      A == 0 & S == 1 ~ "RCT control",
      A == 0 & S == 0 & sel == 1 ~ "External control (selected)",
      A == 0 & S == 0 & sel == 0 ~ "External control (unselected)"
    ) %>% 
      as_factor %>% 
      fct_relevel(
        "RCT treated", "RCT control", 
        "External control (selected)", "External control (unselected)"
      )
  ) %>% 
  filter(Type != "RCT treated")
#> Warning: There was 1 warning in `mutate()`.
#> ℹ In argument: `Type = `%>%`(...)`.
#> Caused by warning:
#> ! 1 unknown level in `f`: External control (selected)

# Fit quantile regressions to visualize the main range of RCT controls
fit975 <- rq(Y ~ `Sampling Score`, tau = 0.975, data = dat_plot, 
             subset = dat_plot$Type == "RCT control")
dat_plot$pred975 <- predict(fit975, newdata = dat_plot)

fit025 <- rq(Y ~ `Sampling Score`, tau = 0.025, data = dat_plot, 
             subset = dat_plot$Type == "RCT control")
dat_plot$pred025 <- predict(fit025, newdata = dat_plot)

# Plot results
dat_plot %>% 
  ggplot() +
  geom_ribbon(aes(`Sampling Score`, ymin = pred025, ymax = pred975), 
              fill = "grey80", alpha = 0.5) +
  geom_point(aes(`Sampling Score`, Y, color = Type)) +
  scale_color_manual(values = c("#5A5A5A", "#00ADFA", "#F8766D"))
```

<img src="man/figures/README-example-1.svg" width="100%" />

The above figure shows that Conformal Selective Borrowing discards
uncomparable ECs and borrows most of the comparable ECs, except for a
few suspicious ones.

## Reference

Zhu, Ke, Shu Yang, and Xiaofei Wang. (2024) “Enhancing Statistical
Validity and Power in Hybrid Controlled Trials: A Randomization
Inference Approach with Conformal Selective Borrowing.” arXiv.
<http://arxiv.org/abs/2410.11713>.
