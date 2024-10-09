
<!-- README.md is generated from README.Rmd. Please edit that file -->

# intFRT

<!-- badges: start -->
<!-- badges: end -->

The goal of **intFRT** is to **int**egrate randomized controlled trials
(RCTs) with external controls (ECs) in hybrid controlled trials,
harnessing Fisher randomization tests (**FRT**) and Conformal Selective
Borrowing (CSB).

It improves the statistical efficiency of average treatment effect (ATE)
estimation and inference while ensuring valid hypothesis testing.

The package offers functions to perform FRTs in hybrid trials, ensuring
strict control of the Type I error rate.

It implements CSB as both an estimator of ATE and a test statistic for
FRT, enabling selective borrowing of comparable ECs to mitigate hidden
bias and enhance statistical power.

It also includes a function to adaptively determine the selection
threshold for CSB.

Additionally, the package provides No Borrowing and various EC borrowing
estimators (IPW, staIPW, CW, OM, AIPW, ACW), along with their inference
results based on asymptotic normality, for comparison.

## Installation

You can install the development version of intFRT from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("ke-zhu/intFRT")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(intFRT)
# This example illustrates the use of Fisher Randomization Test (FRT) and
# different borrowing methods for hybrid controlled trials.

# Simulate data for a hybrid controlled trial
set.seed(123)
n_rct <- 50  # Number of observations in randomized controlled trial
n_ec <- 100  # Number of external controls
n <- n_rct + n_ec  # Total number of observations

# Covariates (2 covariates, normally distributed)
X <- matrix(rnorm(n * 2), n, 2)

# Treatment assignment (1 = treatment, 0 = control)
A <- c(rep(0:1, each = n_rct / 2), rep(0, n_ec))

# Data source indicator (1 = randomized trial, 0 = external control)
S <- c(rep(1, n_rct), rep(0, n_ec))

# Generate outcome (continuous)
Y <- 1 + 0.5 * A + X[,1] + X[,2] + rnorm(n)

# Perform Fisher Randomization Test and No Borrowing
result_nb <- ec_borrow(
  Y = Y,
  A = A,
  S = S,
  X = X,
  method = "No Borrow AIPW",
  family = "gaussian",
  n_fisher = 10  # FRT with 10 randomizations for illustration
)

# View results
print(result_nb$res)
#> # A tibble: 2 × 9
#>   method                est     se   ci_l  ci_u p_value n_sel ess_sel runtime
#>   <chr>               <dbl>  <dbl>  <dbl> <dbl>   <dbl> <dbl>   <dbl>   <dbl>
#> 1 No Borrow AIPW      0.503  0.316 -0.115  1.12   0.111     0       0  0.0260
#> 2 No Borrow AIPW+FRT NA     NA     NA     NA      0.273    NA      NA  0.114

# View ID of borrowed external controls (For No Borrowing, it is NULL)
result_nb$out$id_sel[[1]]
#> NULL

# Perform Fisher Randomization Test and Conformal Selective Borrowing
result_csb <- ec_borrow(
  Y = Y,
  A = A,
  S = S,
  X = X,
  method = "Conformal Selective Borrow AIPW",
  family = "gaussian",
  n_fisher = 10  # FRT with 10 randomizations for illustration
)

# View results
print(result_csb$res)
#> # A tibble: 2 × 9
#>   method                 est     se    ci_l   ci_u p_value n_sel ess_sel runtime
#>   <chr>                <dbl>  <dbl>   <dbl>  <dbl>   <dbl> <dbl>   <dbl>   <dbl>
#> 1 Conformal Selectiv…  0.467  0.223  0.0295  0.905  0.0364    63    44.7  0.0860
#> 2 Conformal Selectiv… NA     NA     NA      NA      0.0909    NA    NA    0.461

# View ID of borrowed external controls (For No Borrowing, it is NULL)
result_csb$out$id_sel[[1]]
#>  [1]  51  53  54  55  56  57  58  61  62  63  66  67  68  70  73  74  75  76  77
#> [20]  78  79  81  82  83  85  87  89  90  91  93  94  95 100 101 103 104 105 107
#> [39] 109 110 111 112 114 115 117 118 119 120 122 123 126 128 131 132 133 135 138
#> [58] 142 143 144 145 146 147
```
