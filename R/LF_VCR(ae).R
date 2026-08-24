
#' Fit LF-VCR with autoencoder latent factors
#'
#' Uses a manually supplied factor count or estimates it with
#' [GrFA::est_num()], uses the resulting `p_hat` as the hidden-layer width of
#' an H2O autoencoder, constructs the same factor-by-feature interactions as
#' [LF_VCR()], and fits the same cross-validated sparse group lasso model.
#'
#' @inheritParams LF_VCR
#'
#' @return A list containing the selected factor count `p_hat`, fitted
#'   `autoencoder`, extracted `factors`, cross-validated sparse group lasso
#'   `model`, its `beta` coefficients at `lambda.min`, and binary
#'   `outcome_levels` when applicable. `p_hat_source` reports whether the
#'   factor count was `"estimated"` or `"manual"`; `kmax_used` reports the
#'   estimation limit.
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' F <- matrix(rnorm(200), nrow = 100)
#' L <- matrix(rnorm(40), nrow = 20)
#' X <- F %*% t(L) + matrix(rnorm(2000, sd = 0.3), nrow = 100)
#' y <- F[, 1] + rnorm(100)
#' fit <- LF_VCR_ae(X, y, nfold = 5)
#' fit$p_hat
#' fit$beta
#' }
LF_VCR_ae <- function(X, y, covariate = NULL, nfold = 10, p_hat = NULL,
                      kmax = 15,
                      factor_criterion = "IC2", categorical = FALSE) {
  inputs <- .lfvcr_prepare_inputs(
    X, y, covariate, nfold, p_hat, kmax, factor_criterion, categorical
  )

  if (inputs$p_hat > 0L) {
    h2o::h2o.init()
    training_frame <- h2o::as.h2o(inputs$X)
    ae_model <- h2o::h2o.deeplearning(
      x = seq_len(ncol(inputs$X)),
      training_frame = training_frame,
      ignore_const_cols = FALSE,
      activation = "Tanh",
      hidden = inputs$p_hat,
      reproducible = TRUE,
      seed = 1,
      autoencoder = TRUE
    )
    factors <- as.matrix(h2o::h2o.deepfeatures(
      ae_model, training_frame, layer = 1
    ))
    if (!identical(dim(factors), c(nrow(inputs$X), inputs$p_hat))) {
      stop("The autoencoder returned an unexpected factor-matrix dimension.")
    }
  } else {
    ae_model <- NULL
    factors <- matrix(numeric(), nrow = nrow(inputs$X), ncol = 0L)
  }
  regression <- .lfvcr_fit_model(inputs, factors)

  return(list(
    p_hat = inputs$p_hat,
    p_hat_source = inputs$p_hat_source,
    kmax_used = inputs$kmax_used,
    autoencoder = ae_model,
    factors = factors,
    model = regression$model,
    beta = regression$beta,
    outcome_levels = regression$outcome_levels
  ))
}

#' @rdname LF_VCR_ae
#' @export
LF_VCR.ae <- LF_VCR_ae
