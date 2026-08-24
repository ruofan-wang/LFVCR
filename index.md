# LFVCR

LFVCR (**L**atent **F**actor **V**arying-**C**oefficient **R**egression)
is an R package for continuous and binary outcomes with high-dimensional
predictors. It estimates a low-dimensional latent structure, constructs
factor-by-predictor interactions, and fits the expanded design using
cross-validated sparse group lasso.

## Method Overview

Both LFVCR implementations follow the same workflow:

1.  Center `X` without scaling its columns.
2.  Estimate the number of latent factors, `p_hat`, with
    [`GrFA::est_num()`](https://rdrr.io/pkg/GrFA/man/est_num.html). The
    default criterion is `"BIC3"`.
3.  Extract exactly `p_hat` factors:
    - [`LF_VCR()`](https://ruofan-wang.github.io/LFVCR/reference/LF_VCR.md)
      uses PCA with `scale. = FALSE`.
    - [`LF_VCR_ae()`](https://ruofan-wang.github.io/LFVCR/reference/LF_VCR_ae.md)
      uses a single-hidden-layer H2O autoencoder with `hidden = p_hat`.
4.  Construct every predictor-by-factor interaction.
5.  Combine the interactions, original predictors, and optional
    covariates.
6.  Fit Gaussian or binomial sparse group lasso with cross-validation.

`p_hat` is estimated automatically in both functions. There is no
`number.K` argument.

| Function | Factor extraction | Supported outcomes |
|----|----|----|
| [`LF_VCR()`](https://ruofan-wang.github.io/LFVCR/reference/LF_VCR.md) | PCA without scaling | Continuous and binary |
| [`LF_VCR_ae()`](https://ruofan-wang.github.io/LFVCR/reference/LF_VCR_ae.md) | H2O autoencoder | Continuous and binary |

## Installation

``` r

# If needed: install.packages("remotes")
remotes::install_github("ruofan-wang/LFVCR")
```

The package dependencies are installed automatically.
[`LF_VCR_ae()`](https://ruofan-wang.github.io/LFVCR/reference/LF_VCR_ae.md)
also requires Java because it starts or reuses a local H2O cluster.

## Example Data

The examples below use a simulated data set with two latent factors and
two adjustment covariates.

``` r

library(LFVCR)

set.seed(42)
n <- 120
p <- 24

true_factors <- matrix(rnorm(n * 2), nrow = n)
loadings <- matrix(rnorm(p * 2), nrow = p)
X <- true_factors %*% t(loadings) +
  matrix(rnorm(n * p, sd = 0.35), nrow = n)

covariate <- cbind(
  age = rnorm(n),
  sex = rbinom(n, 1, 0.5)
)
```

`X` must be a numeric matrix with observations in rows and predictors in
columns. `covariate`, when supplied, must also be numeric and have the
same number of rows as `X`.

## Continuous Outcome

``` r

y_continuous <-
  1.2 * true_factors[, 1] -
  0.7 * true_factors[, 2] +
  0.3 * covariate[, "age"] +
  rnorm(n)

fit_pca <- LF_VCR(
  X = X,
  y = y_continuous,
  covariate = covariate,
  nfold = 5,
  kmax = 8,
  factor_criterion = "BIC3"
)

fit_pca$p_hat
fit_pca$beta
```

Use the autoencoder version with the same arguments. `p_hat` is still
estimated by
[`GrFA::est_num()`](https://rdrr.io/pkg/GrFA/man/est_num.html); it
becomes the autoencoder hidden-layer width.

``` r

fit_ae <- LF_VCR_ae(
  X = X,
  y = y_continuous,
  covariate = covariate,
  nfold = 5,
  kmax = 8,
  factor_criterion = "BIC3"
)

fit_ae$p_hat
fit_ae$beta
```

## Binary Outcome

Set `categorical = TRUE` and supply an outcome with exactly two observed
levels. LFVCR maps the first factor level to 0 and the second factor
level to 1. Set the factor levels explicitly when the positive class
matters.

``` r

probability <- plogis(
  true_factors[, 1] -
  0.8 * true_factors[, 2] +
  0.2 * covariate[, "sex"]
)

y_binary <- factor(
  ifelse(rbinom(n, 1, probability) == 1, "case", "control"),
  levels = c("control", "case")
)

fit_binary <- LF_VCR(
  X = X,
  y = y_binary,
  covariate = covariate,
  nfold = 5,
  kmax = 8,
  factor_criterion = "BIC3",
  categorical = TRUE
)

fit_binary$outcome_levels  # "control" is 0; "case" is 1
fit_binary$beta
```

The same binary-outcome call works with
[`LF_VCR_ae()`](https://ruofan-wang.github.io/LFVCR/reference/LF_VCR_ae.md):

``` r

fit_binary_ae <- LF_VCR_ae(
  X = X,
  y = y_binary,
  covariate = covariate,
  nfold = 5,
  kmax = 8,
  factor_criterion = "BIC3",
  categorical = TRUE
)
```

Binary cross-validation folds are stratified by outcome level.

## Arguments

| Argument | Description |
|----|----|
| `X` | Numeric `n` by `p` predictor matrix. Missing or infinite values are not allowed. |
| `y` | Numeric response for a continuous model, or a two-level response for a binary model. |
| `covariate` | Optional numeric adjustment matrix with `n` rows. |
| `nfold` | Number of cross-validation folds; must be between 2 and `n`. |
| `kmax` | Maximum factor count considered by [`GrFA::est_num()`](https://rdrr.io/pkg/GrFA/man/est_num.html); default is 8. |
| `factor_criterion` | Factor-number criterion: `PC1`, `PC2`, `PC3`, `IC1`, `IC2`, `IC3`, `AIC3`, `BIC3`, `ER`, or `GR`. |
| `categorical` | Use `FALSE` for a Gaussian model and `TRUE` for a binomial model. |

## Returned Values

Both functions return:

| Value | Description |
|----|----|
| `p_hat` | Factor count estimated by [`GrFA::est_num()`](https://rdrr.io/pkg/GrFA/man/est_num.html). |
| `factors` | Extracted `n` by `p_hat` factor matrix. |
| `model` | Cross-validated `sparsegl` model. |
| `beta` | Coefficients evaluated at `lambda.min`. |
| `outcome_levels` | Binary level order; `NULL` for a continuous model. |

[`LF_VCR()`](https://ruofan-wang.github.io/LFVCR/reference/LF_VCR.md)
additionally returns the fitted `pca` object.
[`LF_VCR_ae()`](https://ruofan-wang.github.io/LFVCR/reference/LF_VCR_ae.md)
additionally returns the fitted H2O `autoencoder` object. If `p_hat` is
zero, factor extraction is skipped and the model is fitted with the
original predictors and optional covariates.

## Documentation

The complete function reference is available on the [LFVCR package
website](https://ruofan-wang.github.io/LFVCR/).
