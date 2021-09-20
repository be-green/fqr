
post_processed_grad_descent = function(X, y,
                                       tau, lambda,
                                       nwarmup_samples = 0.1 * nrow(X),
                                       lp_size = 50000) {
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
  init_beta = rep(0, ncol(X))

  init_fit = fit_approx_quantile_model(X,
                            y,
                            X_sub,
                            y_sub,
                            tau,
                            init_beta,
                            maxiter = 1000,
                            mu = 1e-15,
                            beta_tol = 1e-5,
                            check_tol = 1e-5,
                            intercept = 1,
                            num_samples = nwarmup_samples,
                            warm_start = 1, scale = 1, lambda,
                            min_delta = 1e-10)

  res = as.vector(init_fit$residuals)
  checked_res = check(res, tau)

  thresh = quantile(checked_res, lp_size / n)

  globbed_x = glob_obs_mat(X, res, thresh)
  globbed_y = glob_obs_vec(y, res, thresh)
  which_globbed = which(res > thresh | res < -thresh)

  if(lambda > 0) {
    # based on the quantreg implementation rq.lasso.fit
    lambdan = lambda * nrow(X)
    lambdan <- c(0, rep(lambdan, ncol(globbed_x) - 1))
    R <- diag(lambdan, nrow = length(lambdan))
    R <- R[which(lambdan != 0), , drop = FALSE]
    r <- rep(0, nrow(R))
    neg_R  <- -R
    globbed_x = rbind(globbed_x, R, neg_R)
    globbed_y = c(globbed_y, r, r)
  }

  suppressWarnings({
     new_fit = quantreg::rq.fit.br(globbed_x, globbed_y, tau)
  })
  if(lambda > 0) {
    added_rows = (nrow(globbed_x) - nrow(R) * 2 + 1):nrow(globbed_x)
    globbed_x = globbed_x[-added_rows,]
    globbed_y = globbed_y[-added_rows]
  }

  new_res = y - X %*% coef(new_fit)


  # if any residual now has the wrong sign, repeat the process w/ a
  # stricter tolerance for beta
  # if they don't, then the constraint is dominated
  i = 0
  max_times = 10
  while(any(sign(new_res[which_globbed]) != sign(res[which_globbed])) & i < max_times) {
    i = i + 1
    init_fit = fit_approx_quantile_model(X,
                                         y,
                                         X_sub,
                                         y_sub,
                                         tau,
                                         new_fit$coefficients,
                                         mu = 1e-15,
                                         maxiter = 1000,
                                         beta_tol = 1e-10,
                                         check_tol = 1e-10,
                                         intercept = 1,
                                         num_samples = nwarmup_samples,
                                         warm_start = 0, scale = 1, lambda,
                                         min_delta = 1e-10)

    res = as.vector(init_fit$residuals)
    checked_res = check(res)
    thresh = quantile(checked_res, lp_size / n)
    globbed_x = glob_obs_mat(X, res, thresh)
    globbed_y = glob_obs_vec(y, res, thresh)

    which_globbed = which(res > thresh | res < -thresh)

    if(lambda > 0) {
      # based on the quantreg implementation rq.lasso.fit
      lambdan = lambda * nrow(X)
      lambdan <- c(0, rep(lambdan, ncol(globbed_x) - 1))
      R <- diag(lambdan, nrow = length(lambdan))
      R <- R[which(lambdan != 0), , drop = FALSE]
      r <- rep(0, nrow(R))
      neg_R  <- -R
      globbed_x = rbind(globbed_x, R, neg_R)
      globbed_y = c(globbed_y, r, r)
    }

    suppressWarnings({
      new_fit = quantreg::rq.fit.br(globbed_x, globbed_y, tau)
    })

    if(lambda > 0) {
      added_rows = (nrow(globbed_x) - nrow(R) * 2 + 1):nrow(globbed_x)
      globbed_x = globbed_x[-added_rows,]
      globbed_y = globbed_y[-added_rows]
    }

    new_res = y - X %*% coef(new_fit)
  }
  coef(new_fit)
}
