.lfvcr_prepare_inputs <- function(X, y, covariate, nfold, kmax,
                                  factor_criterion, categorical) {
  X <- as.matrix(X)
  if (!is.numeric(X) || length(dim(X)) != 2L) {
    stop("`X` must be a numeric matrix.")
  }
  if (anyNA(X) || any(!is.finite(X))) {
    stop("`X` must not contain missing or non-finite values.")
  }

  n <- nrow(X)
  p <- ncol(X)
  if (n < 2L || p < 2L) {
    stop("`X` must have at least two rows and two columns.")
  }
  if (length(y) != n) {
    stop("`y` must have one value for each row of `X`.")
  }
  if (length(categorical) != 1L || is.na(categorical) ||
      !is.logical(categorical)) {
    stop("`categorical` must be either TRUE or FALSE.")
  }

  if (length(nfold) != 1L || is.na(nfold) || nfold != as.integer(nfold) ||
      nfold < 2L || nfold > n) {
    stop("`nfold` must be an integer between 2 and nrow(X).")
  }
  nfold <- as.integer(nfold)

  if (length(kmax) != 1L || is.na(kmax) || kmax != as.integer(kmax) ||
      kmax < 1L || kmax >= min(n, p)) {
    stop("`kmax` must be an integer of at least 1 and smaller than both dimensions of `X`.")
  }
  kmax <- as.integer(kmax)

  valid_criteria <- c(
    "PC1", "PC2", "PC3", "IC1", "IC2", "IC3",
    "AIC3", "BIC3", "ER", "GR"
  )
  if (length(factor_criterion) != 1L ||
      !factor_criterion %in% valid_criteria) {
    stop(
      "`factor_criterion` must be one of: ",
      paste(valid_criteria, collapse = ", "), "."
    )
  }

  if (is.null(covariate)) {
    covariate_matrix <- NULL
  } else {
    covariate_matrix <- as.matrix(covariate)
    if (!is.numeric(covariate_matrix)) {
      stop("`covariate` must be a numeric matrix.")
    }
    if (nrow(covariate_matrix) != n) {
      stop("`covariate` must have the same number of rows as `X`.")
    }
    if (anyNA(covariate_matrix) || any(!is.finite(covariate_matrix))) {
      stop("`covariate` must not contain missing or non-finite values.")
    }
  }

  if (categorical) {
    if (anyNA(y)) {
      stop("Binary `y` must not contain missing values.")
    }
    y_factor <- droplevels(factor(y))
    if (nlevels(y_factor) != 2L) {
      stop("For a binomial model, `y` must have exactly two observed levels.")
    }
    y_model <- as.numeric(y_factor) - 1
    if (any(table(y_model) < 2L)) {
      stop("Each binary outcome level must contain at least two observations.")
    }
    outcome_levels <- levels(y_factor)
  } else {
    if (!is.numeric(y) || anyNA(y) || any(!is.finite(y))) {
      stop("For a Gaussian model, `y` must be numeric and finite.")
    }
    y_model <- as.numeric(y)
    outcome_levels <- NULL
  }

  X_centered <- scale(X, center = TRUE, scale = FALSE)
  p_hat <- as.integer(GrFA::est_num(
    X_centered, kmax = kmax, type = factor_criterion
  ))
  if (length(p_hat) != 1L || is.na(p_hat) || p_hat < 0L || p_hat > kmax) {
    stop("`GrFA::est_num()` returned an invalid factor count.")
  }

  list(
    X = X,
    y = y_model,
    covariate = covariate_matrix,
    nfold = nfold,
    categorical = categorical,
    outcome_levels = outcome_levels,
    p_hat = p_hat
  )
}

.lfvcr_build_design <- function(X, factors, covariate = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  p_hat <- ncol(factors)

  if (nrow(factors) != n) {
    stop("The factor matrix must have one row for each row of `X`.")
  }

  if (p_hat == 0L) {
    interactions <- matrix(numeric(), nrow = n, ncol = 0L)
  } else {
    interactions <- do.call(
      cbind,
      lapply(seq_len(p), function(j) factors * X[, j])
    )
  }

  design <- cbind(interactions, X, covariate)
  groups <- c(
    rep(seq_len(p), each = p_hat),
    p + seq_len(p + if (is.null(covariate)) 0L else ncol(covariate))
  )
  groups <- match(groups, unique(groups))

  if (ncol(design) != length(groups)) {
    stop("Internal error: sparse-group assignments do not match the design matrix.")
  }

  list(x = design, groups = groups)
}

.lfvcr_fit_model <- function(inputs, factors) {
  design <- .lfvcr_build_design(
    inputs$X, factors, covariate = inputs$covariate
  )

  foldid <- NULL
  if (inputs$categorical) {
    foldid <- integer(length(inputs$y))
    for (class_value in 0:1) {
      class_index <- which(inputs$y == class_value)
      foldid[class_index] <- sample(
        rep(seq_len(inputs$nfold), length.out = length(class_index))
      )
    }
  }

  cv_fit <- sparsegl::cv.sparsegl(
    design$x,
    inputs$y,
    group = design$groups,
    family = if (inputs$categorical) "binomial" else "gaussian",
    nfolds = inputs$nfold,
    foldid = foldid
  )

  list(
    model = cv_fit,
    beta = stats::coef(cv_fit, s = "lambda.min"),
    outcome_levels = inputs$outcome_levels
  )
}
