
#' Fit LF-VCR with autoencoder latent factors
#'
#' Extracts latent factors with an H2O autoencoder, constructs
#' factor-by-feature interactions, and fits the expanded design with
#' cross-validated sparse group lasso.
#'
#' @inheritParams LF_VCR
#' @param number.K The number of autoencoder latent factors to extract.
#'
#' @return A list containing the fitted autoencoder `autoencoder`, extracted
#'   `factors`, cross-validated sparse group lasso `model`, and its `beta`
#'   coefficients at `lambda.min`.
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rnorm(1000), nrow = 100)
#' y <- X[, 1] + rnorm(100)
#' fit <- LF_VCR_ae(X, y, number.K = 2, nfold = 5)
#' fit$beta
#' }
LF_VCR_ae <- function(X, y, number.K, covariate = NULL, nfold = 10,
                      categorical = FALSE) {
  X <- as.matrix(X)
  p <- ncol(X)
  n <- nrow(X)
  if (is.null(covariate)) {
    p1 <- 0
  } else {
    covariate <- as.matrix(covariate)
    p1 <- ncol(covariate)
  }
  h2o::h2o.init()
  training_frame <- h2o::as.h2o(X)
  ae_model <- h2o::h2o.deeplearning(
    x = 1:p,
    training_frame = training_frame,
    ignore_const_cols = FALSE,
    activation = "Tanh",
    hidden = c(number.K),
    reproducible = TRUE,
    seed = 1,
    autoencoder = TRUE
  )
  Z.hat <- t(as.matrix(h2o::h2o.deepfeatures(
    ae_model, training_frame, layer = 1
  )))
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
    autoencoder = ae_model,
    factors = Z.t,
    model = cv_fit,
    beta = beta
  ))
}

#' @rdname LF_VCR_ae
#' @export
LF_VCR.ae <- LF_VCR_ae
