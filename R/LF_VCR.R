
#' LF-VCR
#'
#' extracts the latet factors and fit them with sparse group lasso
#'
#' @param X the data matrix
#' @param number.K the number of factors
#' @param covariate the covariate to be adjusted
#' @param matrix the option of thresholding either correlation or covairance matrix
#' @param categorical denoting whether y is continuous or categorical
#' @param y the outcome
#'
#'
LF_VCR <- function(X, y, number.K, covariate = NULL, nfold=10, matrix='vad', categorical=FALSE) {
  
  X=as.matrix(X)
  p=ncol(X)
  n=nrow(X)
  if (is.null(covariate)) {
    p1 <- 0
  } else {
    covariate <- as.matrix(covariate)
    p1 <- ncol(covariate)
  }
  
  Z.hat<-POET(t(X),number.K,0.5,'soft',matrix)$factors
  Z.t <- t(Z.hat) 
  result <- matrix(nrow = n, ncol = p*number.K)
  for (i in 1:n) {
    temp <- numeric()
    for (j in 1:p) {
      temp <- c(temp, X[i, j] * Z.t[i, ])
    }
    result[i, ] <- temp
  }
  XZ=result
  if (is.null(covariate)) {
    cbind.X.total <- cbind(XZ, X)
    groups <- c(rep(1:p, each = number.K), (p + 1):(p + p))
  } else {
    cbind.X.total <- cbind(XZ, X, covariate)
    groups <- c(rep(1:p, each = number.K), (p + 1):(p + p + p1))
  }
  
  if (!categorical) {
    # If y is not categorical, run without family argument
    cv_fit <- cv.sparsegl(cbind.X.total, y, group = groups, nfolds = nfold)
  } else {
    # Check that y has exactly two levels for a binomial model
    if (length(unique(y)) != 2) {
      stop("Error: For a binomial model, y must have exactly two levels.")
    }
    cv_fit <- cv.sparsegl(cbind.X.total, y, group = groups, nfolds = nfold, family = 'binomial')
  }
  
  return(list(
    model = cv_fit,
    beta = beta
  ))
}
