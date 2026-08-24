# Fit LF-VCR with autoencoder latent factors

Extracts latent factors with an H2O autoencoder, constructs
factor-by-feature interactions, and fits the expanded design with
cross-validated sparse group lasso.

## Usage

``` r
LF_VCR_ae(X, y, number.K, covariate = NULL, nfold = 10, categorical = FALSE)

LF_VCR.ae(X, y, number.K, covariate = NULL, nfold = 10, categorical = FALSE)
```

## Arguments

- X:

  A numeric matrix with observations in rows and predictors in columns.

- y:

  A numeric outcome vector with one value per row of `X`.

- number.K:

  The number of autoencoder latent factors to extract.

- covariate:

  An optional matrix of adjustment covariates.

- nfold:

  The number of cross-validation folds used by
  [`sparsegl::cv.sparsegl()`](https://dajmcdon.github.io/sparsegl/reference/cv.sparsegl.html).

- categorical:

  Logical; use a binomial model when `TRUE` and a Gaussian model when
  `FALSE`.

## Value

A list containing the fitted autoencoder `autoencoder`, extracted
`factors`, cross-validated sparse group lasso `model`, and its `beta`
coefficients at `lambda.min`.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(rnorm(1000), nrow = 100)
y <- X[, 1] + rnorm(100)
fit <- LF_VCR_ae(X, y, number.K = 2, nfold = 5)
fit$beta
} # }
```
