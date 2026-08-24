# LFVCR

LFVCR (**L**atent **F**actor **V**arying-**C**oefficient **R**egression)
is an R package for regression with high-dimensional predictors whose
effects may vary with a low-dimensional latent structure. It estimates
latent factors, constructs factor-by-feature interactions, and fits the
expanded design with sparse group lasso.

Two factor-extraction methods are available:

- [`LF_VCR()`](https://ruofan-wang.github.io/LFVCR/reference/LF_VCR.md)
  estimates factors with Principal Orthogonal complEment Thresholding
  (POET).
- [`LF_VCR_ae()`](https://ruofan-wang.github.io/LFVCR/reference/LF_VCR_ae.md)
  estimates factors with a single-hidden-layer autoencoder.

## Installation

``` r

# If needed: install.packages("remotes")
remotes::install_github("ruofan-wang/LFVCR")
```

## Quick Start

``` r

library(LFVCR)

set.seed(1)
n <- 100
p <- 20
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
y <- 0.8 * X[, 1] - 0.5 * X[, 2] + rnorm(n)

fit <- LF_VCR(
  X = X,
  y = y,
  number.K = 2,
  nfold = 5,
  matrix = "vad"
)

fit$beta
```

For a binary outcome, set `categorical = TRUE` and provide an outcome
with exactly two distinct values. Additional adjustment variables can be
supplied through `covariate`.

## Main Functions

| Function | Latent-factor method | Model |
|----|----|----|
| [`LF_VCR()`](https://ruofan-wang.github.io/LFVCR/reference/LF_VCR.md) | POET factor model | Sparse group lasso |
| [`LF_VCR_ae()`](https://ruofan-wang.github.io/LFVCR/reference/LF_VCR_ae.md) | H2O autoencoder | Sparse group lasso |

Both functions return the cross-validated sparse group lasso model and
its coefficients at `lambda.min`.

## Documentation

The complete function reference is available on the [LFVCR package
website](https://ruofan-wang.github.io/LFVCR/).
