## usethis namespace: start
#' @useDynLib fqr, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL


## usethis namespace: start
#' @importFrom RcppParallel RcppParallelLibs
## usethis namespace: end
NULL

#' Quantile Loss Function
#' @param x residuals to be evaluated
#' @param tau quantile to use
check <- function(x,tau=.5){
  x * (tau - (x < 0))
}


#' Finds which column of X has the intercept term
#' @param X design matrix for regression
get_intercept <- function(X) {
  which_cols <- 1:ncol(X)
  i = 1
  while(length(which_cols) > 1 & i <= nrow(X)) {
    which_cols <- which(X[i, which_cols] == 1)
    i = i + 1
  }

  if (length(which_cols) > 1) {
    stop("Two intercept columns specified in design matrix. Please remove one.")
  } else if (length(which_cols) == 0) {
    0
  } else {
    if(all(X[, which_cols] == 1)) {
      which_cols
    } else {
      0
    }
  }
}

#' @rdname fqr
#' @param y outcome variable
#' @param X design matrix
#' @param tau target quantile
#' @param se whether to calculate standard errors
#' @param init_beta initial betas for gradient descent
#' (optional, default is random normal initial values)
#' @param smoothing_window neighborhood around 0 to smooth with
#' tilted least-squares loss function
#' @param beta_tol stopping criterion based on largest value of the
#' gradient
#' @param check_tol stopping criterion based on the change in the
#' value of the check function between iterations
#' @param intercept what column the intercept is, defaults to 1, 0 to indicate
#' no intercept
#' @param nwarmup_samples number of samples to use for warmup regression
#' @param nsubsamples number of subsamples to use when calculating standard errors
#' @param warm_start whether to run initial warmup regression or just
#' default inits for full data gradient descent
#' @param maxiter maximum number of allowed iterations for gradient descent
#' @param labmda lasso penalty parameter for penalized regression
#' @export
#' @importFrom stats model.matrix
#' @importFrom stats rnorm
#' @importFrom stats cov
fit_fqr <- function(X, y, tau,
                    se = T,
                    init_beta = rep(0, ncol(X)),
                    intercept = 1,
                    nsubsamples = 100,
                    nwarmup_samples = max(0.1 * nrow(X), 100),
                    lambda = 0,
                    warm_start = 1,
                    lp_size = 1000) {

  n <- nrow(X)
  if(nwarmup_samples > n) {
    nwarmup_samples <- min(n, nwarmup_samples/10)
  }
  if(lp_size > n) {
    lp_size <- min(n, lp_size/10)
  }

  model_fits <- list()
  boot_list <- list()

  # if it's too small we just bootstrap
  subsample_size <- min(min(0.2 * n, 100), n)

  ss_list <- list()
  for(t in 1:length(tau)) {

    model_fits[[t]] <- post_processed_grad_descent(X, y,
                                tau, lambda = lambda,
                                nwarmup_samples = nwarmup_samples,
                                lp_size = lp_size)

    init_beta = model_fits[[t]]

    if(se) {
      ss_mat <- matrix(NA, nrow = nsubsamples,
                       ncol = length(model_fits[[t]]))
      for(i in 1:nsubsamples) {
        ss_rows <- sample(1:n, subsample_size)

        ss_mat[i,] <- post_processed_grad_descent(X[ss_rows,], y[ss_rows],
                                                  tau, lambda = lambda,
                                                  nwarmup_samples = nwarmup_samples,
                                                  lp_size = lp_size)

        ss_list[[t]] <- ss_mat
      }

    }


  }

    est_betas <- do.call("cbind", model_fits)
    res <- do.call("cbind", args = lapply(model_fits, function(b) y - X %*% b))

    if(se) {

      ss_percent = subsample_size/n

      vcov_list <- lapply(ss_list, stats::cov)
      vcov_list <- lapply(vcov_list, function(x) x * ss_percent)

      standard_errors <- do.call("cbind",
                                 lapply(vcov_list,
                                        function(x) sqrt(diag(x))))
      structure(list(
        coefficients = est_betas,
        standard_errors = standard_errors,
        residuals = res,
        ss_coefs = ss_list,
        ss_percent = ss_percent,
        args = list(
          tau = tau,
          coefnames = colnames(X),
          X = X
        )
      ), class = "fqr")
    } else {
      structure(list(
        coefficients = est_betas,
        standard_errors = matrix(NA, ncol = ncol(est_betas),
                                 nrow = nrow(est_betas)),
        residuals = res,
        ss_coefs = NA,
        ss_percent = NA,
        args = list(
          tau = tau,
          coefnames = colnames(X),
          X = X
        )
      ), class = "fqr")
    }
}

#' Fast Quantile Regression
#' @param formula Regression formula
#' @param data data to use when fitting regression
#' @param tau vector of target quantile(s)
#' @param se whether to calculate standard errors
#' @param init_beta initial betas for gradient descent
#' (optional, default is random normal initial values)
#' @param smoothing_window neighborhood around 0 to smooth with
#' tilted least-squares loss function
#' @param intercept what column the intercept is, defaults to 1, 0 to indicate
#' no intercept
#' @param nwarmup_samples number of samples to use for warmup regression
#' @param nsubsamples number of subsamples to use when calculating standard errors
#' @param warm_start whether to run initial warmup regression or just
#' default inits for full data gradient descent
#' @param maxiter maximum number of allowed iterations for gradient descent
#' @export
#' @importFrom stats model.matrix
#' @importFrom stats terms
#' @details
#' This package performs quantile regression by approximating the check loss
#' function with a least-squares loss function in a small neighbor around 0.
#' Since the only point where the check function is not differentiable is at 0,
#' this allows for first-order gradient descent methods to work.
#'
#' This package uses "accelerated" gradient descent, which moves the future guess
#' at the coefficients not only by a step size * the gradient, but also based
#' on the prior momentum of the changes in the coefficient which leads to
#' faster convergence. Gradient-based methods work at scale (both in terms of
#' observations and dimension), and are much faster than the interior point
#' algorithms in the quantreg package for large problems. Still, they are
#' sometimes less exact for small problems.
#'
#' The algorithm employs two early stopping rules: the `check_tol` argument
#' stops based on the scaled change in the check function loss. The `beta_tol`
#' argument stops based on the largest value of the gradient vector.
#'
#' Before using the full dataset, the optimizer "warms up" on a random subset
#' of the dataset. nwarmup_samples controls the size of that. `warm_start` is
#' an integer which controls whether that happens at all
#' (it is _strongly_ recommended).
#' @examples
#'
#' fit <- fqr(area ~ peri, data = rock, tau = c(0.25, 0.5, 0.75))
#'
#' # print coefficients & SEs
#' print(fit)
#'
#' # grab coefficient vector
#' coef(fit)
#'
#' # predict values
#' predict(fit)
#'
#' # predict values with new data
#' predict(fit, newdata = head(rock))
fqr <- function(formula, data, tau = 0.5,
                se = T,
                nwarmup_samples = pmin(pmax(100, 0.1 * nrow(data)), nrow(data)),
                warm_start = 1,
                nsubsamples = 100) {
  m = stats::model.matrix(formula, data = data)
  intercept = get_intercept(m)
  y = stats::model.response(stats::model.frame(formula, data),
                                type = "numeric")
  fit = fit_fqr(X = m, y = y,
                tau = tau, se = se,
                init_beta = rep(0, ncol(m)),
                intercept = intercept,
                nwarmup_samples = nwarmup_samples,
                warm_start = warm_start,
                nsubsamples = nsubsamples)
  fit$args$formula = formula
  fit
}

#' Get coefficients from fast quantile regression
#' @param object object to get coefficients from
#' @param ... additional arguments, ignored for now
#' @export
#' @importFrom stats coef
coef.fqr <- function(object, ...) {
  object$coefficients
}

#' @rdname coef.fqr
#' @param object object to get coefficients from
#' @param ... additional arguments, ignored for now
#' @importFrom stats coefficients
#' @export
coefficients.fqr <- function(object, ...) {
  coef(object)
}

#' Get model residuals
#' @rdname resid.fqr
#' @param object object to get coefficients from
#' @param ... additional arguments, ignored for now
#' @importFrom stats resid
#' @export
resid.fqr <- function(object, ...) {
  object$residuals
}

#' @rdname resid.fqr
#' @param object object to get coefficients from
#' @param ... additional arguments, ignored for now
#' @importFrom stats residuals
#' @export
residuals.fqr <- function(object, ...) {
  object$residuals
}

#' Print fast quantile regression
#' @param x fit to print
#' @param digits digits to use when printing
#' @param ... other arguments to print
#' @export
print.fqr <- function(x, digits = getOption("digits"), ...) {
  coef_mat <- x$coefficients
  se_mat <- x$standard_errors
  for(i in 1:ncol(coef_mat)) {
    out <- cbind(coef_mat[,i], se_mat[,i])
    colnames(out) <- c("Coefficient", "SE")
    rownames(out) <- x$args$coefnames
    round(out, digits)
    cat("Tau: ", x$args$tau[i], "\n")
    print(out, ...)
  }
}

#' Predict from fast quantile regression fit
#' @param object a model object to use for prediction
#' @param newdata a data.frame with the same columns as data used to fit model
#' @param ... other parameters, ignored for now
#' @export
predict.fqr <- function(object, newdata = NULL, ...) {
  if(is.null(newdata)) {
    object$args$X %*% object$coefficients
  } else if(is.null(object$args$formula)) {
    stop("Must specify formula in fqr call to use predict function with new data")
  } else {
    ff <- stats::as.formula(object$args$formula)
    tt <- stats::terms(ff, data = newdata)
    tt <- stats::delete.response(tt)
    m <- stats::model.matrix(tt, data = newdata)
    m %*% object$coefficients
  }
}
