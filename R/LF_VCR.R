
#' Fit LF-VCR with PCA latent factors
#'
#' Estimates the number of factors with [GrFA::est_num()], extracts that many
#' principal components without scaling the predictors, constructs
#' factor-by-feature interactions, and fits the expanded design with
#' cross-validated sparse group lasso.
#'
#' @param X A numeric matrix with observations in rows and predictors in
#'   columns.
#' @param y A numeric outcome vector with one value per row of `X`.
#' @param covariate An optional matrix of adjustment covariates.
#' @param nfold The number of cross-validation folds used by
#'   [sparsegl::cv.sparsegl()].
#' @param kmax The maximum number of factors considered by
#'   [GrFA::est_num()].
#' @param factor_criterion The criterion passed to [GrFA::est_num()] when
#'   estimating `p_hat`. The default is `"BIC3"`.
#' @param categorical Logical; use a binomial model when `TRUE` and a Gaussian
#'   model when `FALSE`.
#'
#' @return A list containing the estimated factor count `p_hat`, extracted PCA
#'   `factors`, fitted `pca` object, cross-validated sparse group lasso `model`,
#'   and its `beta` coefficients at `lambda.min`.
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' F <- matrix(rnorm(200), nrow = 100)
#' L <- matrix(rnorm(40), nrow = 20)
#' X <- F %*% t(L) + matrix(rnorm(2000, sd = 0.3), nrow = 100)
#' y <- F[, 1] + rnorm(100)
#' fit <- LF_VCR(X, y, nfold = 5)
#' fit$p_hat
#' fit$beta
#' }
LF_VCR <- function(X, y, covariate = NULL, nfold = 10, kmax = 8,
                   factor_criterion = "BIC3", categorical = FALSE) {
  X <- as.matrix(X)
  p <- ncol(X)
  n <- nrow(X)
  if (is.null(covariate)) {
    p1 <- 0
  } else {
    covariate <- as.matrix(covariate)
    p1 <- ncol(covariate)
  }

  if (kmax < 1 || kmax >= min(n, p)) {
    stop("`kmax` must be at least 1 and smaller than both nrow(X) and ncol(X).")
  }

  X.centered <- scale(X, center = TRUE, scale = FALSE)
  p_hat <- GrFA::est_num(
    X.centered, kmax = kmax, type = factor_criterion
  )

  if (p_hat > 0) {
    pca_fit <- stats::prcomp(
      X, center = TRUE, scale. = FALSE, rank. = p_hat
    )
    factors <- pca_fit$x[, seq_len(p_hat), drop = FALSE]
    XZ <- matrix(nrow = n, ncol = p * p_hat)
    for (i in seq_len(n)) {
      temp <- numeric()
      for (j in seq_len(p)) {
        temp <- c(temp, X[i, j] * factors[i, ])
      }
      XZ[i, ] <- temp
    }
  } else {
    pca_fit <- NULL
    factors <- matrix(numeric(), nrow = n, ncol = 0)
    XZ <- matrix(numeric(), nrow = n, ncol = 0)
  }

  if (is.null(covariate)) {
    cbind.X.total <- cbind(XZ, X)
    groups <- c(rep(seq_len(p), each = p_hat), p + seq_len(p))
  } else {
    cbind.X.total <- cbind(XZ, X, covariate)
    groups <- c(
      rep(seq_len(p), each = p_hat),
      p + seq_len(p + p1)
    )
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
    p_hat = p_hat,
    factors = factors,
    pca = pca_fit,
    model = cv_fit,
    beta = beta
  ))
}
