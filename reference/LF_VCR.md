# Fit LF-VCR with PCA latent factors

Estimates the number of factors with
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
  kmax = 8,
  factor_criterion = "BIC3",
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

- kmax:

  The maximum number of factors considered by
  [`GrFA::est_num()`](https://rdrr.io/pkg/GrFA/man/est_num.html).

- factor_criterion:

  The criterion passed to
  [`GrFA::est_num()`](https://rdrr.io/pkg/GrFA/man/est_num.html) when
  estimating `p_hat`. The default is `"BIC3"`.

- categorical:

  Logical; use a binomial model when `TRUE` and a Gaussian model when
  `FALSE`.

## Value

A list containing the estimated factor count `p_hat`, extracted PCA
`factors`, fitted `pca` object, cross-validated sparse group lasso
`model`, and its `beta` coefficients at `lambda.min`.

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
