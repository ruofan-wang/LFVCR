# Fit LF-VCR with POET latent factors

Extracts latent factors with Principal Orthogonal complEment
Thresholding (POET), constructs factor-by-feature interactions, and fits
the expanded design with cross-validated sparse group lasso.

## Usage

``` r
LF_VCR(
  X,
  y,
  number.K,
  covariate = NULL,
  nfold = 10,
  matrix = "vad",
  categorical = FALSE
)
```

## Arguments

- X:

  A numeric matrix with observations in rows and predictors in columns.

- y:

  A numeric outcome vector with one value per row of `X`.

- number.K:

  The number of latent factors to extract.

- covariate:

  An optional matrix of adjustment covariates.

- nfold:

  The number of cross-validation folds used by
  [`sparsegl::cv.sparsegl()`](https://dajmcdon.github.io/sparsegl/reference/cv.sparsegl.html).

- matrix:

  The POET thresholding scale: `"cor"` for the correlation matrix or
  `"vad"` for the covariance matrix.

- categorical:

  Logical; use a binomial model when `TRUE` and a Gaussian model when
  `FALSE`.

## Value

A list containing the cross-validated sparse group lasso `model` and its
`beta` coefficients at `lambda.min`.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(rnorm(1000), nrow = 100)
y <- X[, 1] + rnorm(100)
fit <- LF_VCR(X, y, number.K = 2, nfold = 5)
fit$beta
} # }
```
