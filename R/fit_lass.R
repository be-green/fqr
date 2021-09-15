
#' @param x design matrix
#' @param y target variable
#' @param tau quantile to target
#' @param lambda penalty parameter
#' @param eps determines upper bound of possible target quantiles
#' @export
#' @importFrom quantreg LassoLambdaHat
fit_lasso = function (x, y, tau = 0.5, lambda = NULL,
                      eps = 1e-06)
{
  n <- length(y)
  p <- ncol(x)
  lambda = lambda * n
  if (n != nrow(x))
    stop("x and y don't match n")
  if (length(lambda) == 1)
    lambda <- c(0, rep(lambda, p - 1))
  else if (length(lambda) != p)
    stop(paste("lambda must be either of length ", p, " or length one"))

  if (any(lambda < 0))
    stop("negative lambdas disallowed")
  R <- diag(lambda, nrow = length(lambda))
  R <- R[which(lambda != 0), , drop = FALSE]
  r <- rep(0, nrow(R))
  if (tau < eps || tau > 1 - eps)
    stop("Set tau in (0,1)")
  X <- rbind(x, R)
  Y <- c(y, r)
}
