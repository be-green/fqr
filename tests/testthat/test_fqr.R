n <- 1e4
b <- cbind(matrix(rnorm(n), ncol = 2))
betas <- rnorm(ncol(b))
y <- as.vector(b %*% betas) + rt(nrow(b), df = 1)
y[y > quantile(y, 0.9)] <- quantile(y, 0.9)
data <- data.frame(y = y, b)

fit <- fqr(y ~ X1 + X2, data = data, tau = 0.5, se = F)
fit_with_se <- fqr(y ~ X1 + X2, data = data, tau = 0.5, se = T)
fit_mult_taus <- fqr(y ~ X1 + X2, data = data, tau = c(0.5, 0.7, 0.9), se = F)
fit_custom_num_ss <- fqr(y ~ X1 + X2, data = data,
                         tau = 0.5, se = T, nsubsamples = 10)
fit_no_int <- fqr(y ~ 0 + X1 + X2, data = data,
                         tau = 0.5, se = F)

fit_no_warmup <- fqr(y ~ X1 + X2, data = data,
                  tau = 0.5, se = F, warm_start = F)

testthat::test_that("Regression interface works properly", {
  testthat::expect_equal(fit$standard_errors, matrix(nrow = 3))
  testthat::expect_type(fit_with_se$standard_errors, "double")
  testthat::expect_true(inherits(fit_with_se$standard_errors, "matrix"))
  testthat::expect_s3_class(fit, "fqr")
  testthat::expect_equal(nrow(fit_custom_num_ss$ss_coefs[[1]]), 10)
  testthat::expect_equal(ncol(fit_mult_taus$coefficients), 3)
  testthat::expect_equal(nrow(fit_no_int$coefficients), 2)
  testthat::expect_s3_class(fit_no_warmup, "fqr") # just making sure that this still runs
})

fit_manual <- fit_fqr(b, y, 0.2, se = F)
new <- cbind(matrix(rnorm(n/10), ncol = 2))
testthat::test_that("Predict interface works properly", {
  testthat::expect_equal(nrow(predict(fit_manual)), length(y))
  testthat::expect_error(predict(fit_manual, newdata = new))
  testthat::expect_equal(nrow(predict(fit, newdata = data.frame(new))), nrow(new))
})

testthat::test_that("Coefficient interface works properly", {
  testthat::expect_equal(coef(fit), fit$coefficients)
  testthat::expect_equal(coefficients(fit), fit$coefficients)
})

testthat::test_that("Fast exponential weights work", {
  testthat::expect_lt(var(fast_rexp(100000)), 2)
  testthat::expect_lt(mean(fast_rexp(100000)), 2)
})

testthat::test_that("Print functions generate output", {
  testthat::expect_output(print(fit), "Tau")
})

x = matrix(rnorm(10000, mean = 1, sd = 2), ncol = 1)
z_score(x, 1, 2, p = 1)
testthat::test_that("Z scoring code works", {
  testthat::expect_lt(abs(sd(x) - 1), 0.05)
  testthat::expect_lt(abs(mean(x)), 0.01)

})
