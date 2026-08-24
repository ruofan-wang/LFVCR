# Fit LF-VCR with PCA latent factors

Uses a manually supplied factor count or estimates it with
[`GrFA::est_num()`](https://rdrr.io/pkg/GrFA/man/est_num.html), extracts
that many principal components without scaling the predictors,
constructs factor-by-feature interactions, and fits the expanded design
with cross-validated sparse group lasso.

## Usage

``` r
LF_VCR(
  X,
  y,
  covariate = NULL,
  nfold = 10,
  p_hat = NULL,
  kmax = 15,
  factor_criterion = "IC2",
  categorical = FALSE
)
```

## Arguments

- X:

  A numeric matrix with observations in rows and predictors in columns.

- y:

  A numeric outcome vector with one value per row of `X`.

- covariate:

  An optional matrix of adjustment covariates.

- nfold:

  The number of cross-validation folds used by
  [`sparsegl::cv.sparsegl()`](https://dajmcdon.github.io/sparsegl/reference/cv.sparsegl.html).

- p_hat:

  Optional nonnegative integer specifying the number of factors. The
  default, `NULL`, estimates the number with
  [`GrFA::est_num()`](https://rdrr.io/pkg/GrFA/man/est_num.html). When a
  value is supplied, automatic factor-number estimation is skipped.

- kmax:

  The maximum number of factors considered by
  [`GrFA::est_num()`](https://rdrr.io/pkg/GrFA/man/est_num.html) when
  `p_hat = NULL`. The default is 15. If the data have fewer dimensions,
  the largest allowable value is used.

- factor_criterion:

  The criterion passed to
  [`GrFA::est_num()`](https://rdrr.io/pkg/GrFA/man/est_num.html) when
  estimating `p_hat`. The default is `"IC2"`.

- categorical:

  Logical; use a binomial model when `TRUE` and a Gaussian model when
  `FALSE`. For a binary outcome, the second factor level is modeled as 1
  and both levels are returned in `outcome_levels`.

## Value

A list containing the selected factor count `p_hat`, extracted PCA
`factors`, fitted `pca` object, cross-validated sparse group lasso
`model`, its `beta` coefficients at `lambda.min`, and binary
`outcome_levels` when applicable. `p_hat_source` reports whether the
factor count was `"estimated"` or `"manual"`; `kmax_used` reports the
estimation limit.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
F <- matrix(rnorm(200), nrow = 100)
L <- matrix(rnorm(40), nrow = 20)
X <- F %*% t(L) + matrix(rnorm(2000, sd = 0.3), nrow = 100)
y <- F[, 1] + rnorm(100)
fit <- LF_VCR(X, y, nfold = 5)
fit$p_hat
fit$beta
} # }
```
