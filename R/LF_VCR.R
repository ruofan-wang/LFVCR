
#' Fit LF-VCR with POET latent factors
#'
#' Extracts latent factors with Principal Orthogonal complEment Thresholding
#' (POET), constructs factor-by-feature interactions, and fits the expanded
#' design with cross-validated sparse group lasso.
#'
#' @param X A numeric matrix with observations in rows and predictors in
#'   columns.
#' @param y A numeric outcome vector with one value per row of `X`.
#' @param number.K The number of latent factors to extract.
#' @param covariate An optional matrix of adjustment covariates.
#' @param nfold The number of cross-validation folds used by
#'   [sparsegl::cv.sparsegl()].
#' @param matrix The POET thresholding scale: `"cor"` for the correlation
#'   matrix or `"vad"` for the covariance matrix.
#' @param categorical Logical; use a binomial model when `TRUE` and a Gaussian
#'   model when `FALSE`.
#'
#' @return A list containing the cross-validated sparse group lasso `model` and
#'   its `beta` coefficients at `lambda.min`.
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rnorm(1000), nrow = 100)
#' y <- X[, 1] + rnorm(100)
#' fit <- LF_VCR(X, y, number.K = 2, nfold = 5)
#' fit$beta
#' }
LF_VCR <- function(X, y, number.K, covariate = NULL, nfold = 10,
                   matrix = "vad", categorical = FALSE) {
  X <- as.matrix(X)
  p <- ncol(X)
  n <- nrow(X)
  if (is.null(covariate)) {
    p1 <- 0
  } else {
    covariate <- as.matrix(covariate)
    p1 <- ncol(covariate)
  }

  Z.hat <- POET::POET(t(X), number.K, 0.5, "soft", matrix)$factors
  Z.t <- t(Z.hat) 
  result <- matrix(nrow = n, ncol = p * number.K)
  for (i in 1:n) {
    temp <- numeric()
    for (j in 1:p) {
      temp <- c(temp, X[i, j] * Z.t[i, ])
    }
    result[i, ] <- temp
  }
  XZ <- result
  if (is.null(covariate)) {
    cbind.X.total <- cbind(XZ, X)
    groups <- c(rep(1:p, each = number.K), (p + 1):(p + p))
  } else {
    cbind.X.total <- cbind(XZ, X, covariate)
    groups <- c(rep(1:p, each = number.K), (p + 1):(p + p + p1))
  }
  
  if (!categorical) {
    # If y is not categorical, run without family argument
    cv_fit <- sparsegl::cv.sparsegl(
      cbind.X.total, y, group = groups, nfolds = nfold
    )
  } else {
    # Check that y has exactly two levels for a binomial model
    if (length(unique(y)) != 2) {
      stop("Error: For a binomial model, y must have exactly two levels.")
    }
    cv_fit <- sparsegl::cv.sparsegl(
      cbind.X.total, y, group = groups, nfolds = nfold,
      family = "binomial"
    )
  }
  beta <- stats::coef(cv_fit, s = "lambda.min")

  return(list(
    model = cv_fit,
    beta = beta
  ))
}
