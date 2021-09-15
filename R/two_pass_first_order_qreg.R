
post_processed_grad_descent = function(X, y,
                                       tau, lambda,
                                       nwarmup_samples = 0.1 * nrow(X),
                                       lp_size = 10000) {
  n <- nrow(X)
  p <- ncol(X)

  if(nwarmup_samples > n) {
    nwarmup_samples <- min(n, nwarmup_samples/10)
  }

  if(nwarmup_samples < p) {
    nwarmup_samples <- p * 2
  }

  samples <- sample(1:n, nwarmup_samples, replace = F)

  X_sub <- X[samples,]
  y_sub <- y[samples]
  init_beta = rnorm(ncol(X))

  tol = 1e-4
  smooth = 0.01
  init_fit = fit_approx_quantile_model(X,
                            y,
                            X_sub,
                            y_sub,
                            tau,
                            init_beta,
                            mu = smooth,
                            maxiter = 1000,
                            beta_tol = tol,
                            check_tol = 0,
                            1,
                            nwarmup_samples,
                            1, lambda)

  pred = X %*% init_fit
  res = y - pred
  checked_res = check(res, tau)

  thresh = quantile(checked_res, lp_size / n)

  globbed_x = glob_obs_mat(X, res, thresh)
  globbed_y = glob_obs_vec(y, res, thresh)
  which_globbed = which(res > thresh | res < -thresh)

  new_fit = quantreg::rq.fit.br(globbed_x, globbed_y, tau)
  new_res = y - X %*% coef(new_fit)

  # if any residual now has the wrong sign, repeat the process w/ a
  # stricter tolerance for beta
  # if they don't, then the constraint is dominated
  i = 0
  max_times = 10
  while(any(sign(new_res[which_globbed]) != sign(res[which_globbed])) & i < max_times) {
    i = i + 1
    tol = tol / 5
    smooth = smooth / 2
    init_fit = fit_approx_quantile_model(X,
                                         y,
                                         X_sub,
                                         y_sub,
                                         tau,
                                         init_fit,
                                         mu = smooth,
                                         maxiter = 1000,
                                         beta_tol = tol,
                                         check_tol = 0,
                                         1,
                                         nwarmup_samples,
                                         0, lambda)

    if(lambda > 0) {

    }

    pred = X %*% init_fit
    res = y - pred
    checked_res = check(res)

    globbed_x = glob_obs_mat(X, res, thresh)
    globbed_y = glob_obs_vec(y, res, thresh)
    which_globbed = which(res > thresh | res < -thresh)

    new_fit = quantreg::rq.fit.br(globbed_x, globbed_y, tau)
    new_res = y - X %*% coef(new_fit)
  }

  coef(new_fit)
}
