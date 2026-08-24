
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
#'   model when `FALSE`. For a binary outcome, the second factor level is
#'   modeled as 1 and both levels are returned in `outcome_levels`.
#'
#' @return A list containing the estimated factor count `p_hat`, extracted PCA
#'   `factors`, fitted `pca` object, cross-validated sparse group lasso `model`,
#'   its `beta` coefficients at `lambda.min`, and binary `outcome_levels` when
#'   applicable.
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
  inputs <- .lfvcr_prepare_inputs(
    X, y, covariate, nfold, kmax, factor_criterion, categorical
  )

  if (inputs$p_hat > 0L) {
    pca_fit <- stats::prcomp(
      inputs$X, center = TRUE, scale. = FALSE, rank. = inputs$p_hat
    )
    factors <- pca_fit$x[, seq_len(inputs$p_hat), drop = FALSE]
  } else {
    pca_fit <- NULL
    factors <- matrix(numeric(), nrow = nrow(inputs$X), ncol = 0L)
  }
  regression <- .lfvcr_fit_model(inputs, factors)

  return(list(
    p_hat = inputs$p_hat,
    factors = factors,
    pca = pca_fit,
    model = regression$model,
    beta = regression$beta,
    outcome_levels = regression$outcome_levels
  ))
}
