
# fqr

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/be-green/fqr/branch/main/graph/badge.svg)](https://codecov.io/gh/be-green/fqr?branch=main)
[![R-CMD-check](https://github.com/be-green/fqr/workflows/R-CMD-check/badge.svg)](https://github.com/be-green/fqr/actions)
<!-- badges: end -->

The fqr package makes quantile regression fast and scaleable using
accelerated gradient descent. While the quantile loss function isn’t
differentiable, you can get an arbitrarily close smooth approximation by
replacing the “check” function with an appropriately tilted least
squares approximation for a small neighborhood around the origin. As the
size of that window goes to zero, you have your check function back!

The package uses 2 stopping rules to assess convergence: the maximum
value of the gradient vector (for the coefficients of the quantile
regression) and the relative change in the loss function (scaled by the
step size).

`fqr` is substantially faster than the `quantreg` package’s interior
point methods (e.g. “br” or “sfn”), especially for large problems. The
algorithm implemented via the Armadillo library for linear algebra in
C++. It also has no dependencies other than base R and (if building from
source) a C++ compiler.

## Installation

You can install the fqr package from github by running

``` r
# get remotes if needed:
# install.packages("remotes")

remotes::install_github("be-green/fqr")
```

## Basic Use

The `fqr` package uses the same basic formula interface that `lm` does,
with standard errors calculated based on subsampling.

``` r
library(fqr)
data(rocks)
#> Warning in data(rocks): data set 'rocks' not found

fqr(area ~ peri, data = rock, tau = c(0.25, 0.5, 0.75))
#> Tau:  0.25 
#>              Coefficient           SE
#> (Intercept) 5.214371e+03 498.53198357
#> peri        4.107717e-02   0.01090661
#> Tau:  0.5 
#>             Coefficient           SE
#> (Intercept) 7350.836410 4.168024e+02
#> peri           0.048653 6.823924e-03
#> Tau:  0.75 
#>              Coefficient           SE
#> (Intercept) 8.740056e+03 3.981218e+02
#> peri        4.033071e-02 5.394272e-03
```

To turn off standard errors (and just get point predictions), you can
set `se = F`.

``` r
fqr(area ~ peri, data = rock, se = F, tau = c(0.25, 0.5, 0.75))
#> Tau:  0.25 
#>              Coefficient SE
#> (Intercept) 5.213284e+03 NA
#> peri        4.190239e-02 NA
#> Tau:  0.5 
#>              Coefficient SE
#> (Intercept) 7.348424e+03 NA
#> peri        4.947822e-02 NA
#> Tau:  0.75 
#>              Coefficient SE
#> (Intercept) 8.737112e+03 NA
#> peri        4.114631e-02 NA
```

## Benchmarks

Ok, but *how* fast is this approach? Let’s just take some point
estimates and see how it goes.

### Medium N, Medium P

But with all of this done, let’s compare to some benchmarks from the
`sfn` and `pfn` algorithms, which are currently the fastest in the
`quantreg` package.

``` r
# simulate some data, 101 x 100,000
p <- 20
n <- 1e6
beta <- rnorm(p + 1)

x <- cbind(1, matrix(rnorm(p * n), ncol = p, nrow = n))
y <- x %*% beta + exp(rnorm(n, sd = 2))

# let's take a look at what this looks like
hist(y)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

Ok so we have some *very* skewed data! Perfect for median regression.

``` r
# remove the intercept since it's implied in the formula
start = proc.time()
# lower level version that just takes design matrix
fit <- fit_fqr(x, y, tau = 0.5, se = F)
end = proc.time()
end - start
#>    user  system elapsed 
#>  17.628   0.357   3.663
```

I attempted to run the same thing with the `quantreg` package, with the
method advised for large datasets, like so:

``` r
# newton interior point method w/ pre-processing
start <- proc.time()
fit_pfn <- quantreg::rq.fit.pfn(x, y, tau = 0.5)
end <- proc.time()
end - start
```

but I killed it after 20 minutes (feel free to try this yourself!).

# Big N, Big P

Let’s benchmark with a bigger set of columns.

``` r
p <- 100
n <- 1e6
beta <- rnorm(p + 1)

x <- cbind(1, matrix(rnorm(p * n), ncol = p, nrow = n))
y <- 10 + x %*% beta + exp(rnorm(n, sd = 2))
```

``` r
start = proc.time()
fit <- fit_fqr(x, y, tau = 0.5, se = F)
end = proc.time()
end - start
#>    user  system elapsed 
#> 164.423   3.336  24.863
```

I’m not going to run the quantreg `pfn` algorithm since it was so slow
for the last problem.

### Big N, Small P

Let’s try a more manageable set of dimensions, with *lots* of
observations.

``` r
p <- 10
n <- 1e7
beta <- rnorm(p + 1)

x <- cbind(1, matrix(rnorm(p * n), ncol = p, nrow = n))
y <- x %*% beta + exp(rnorm(n, sd= 2))

start = proc.time()
fit <- fit_fqr(X = x, y = y, tau = 0.5, se = F)
end = proc.time()
end - start
#>    user  system elapsed 
#> 112.843   1.916  22.045
```

I attempted to do the comparable thing for the `pfn` algorithm:

``` r
start = proc.time()
fit <- quantreg::rq.fit.pfn(x = x, y = y, tau = 0.5)
end = proc.time()
end - start
```

…but I killed the process after 15 minutes or so.

# Medium-scale Problem

Ok, so we haven’t been able to run quantreg on these datasets, let’s see
how it does with a sort of medium-scale problem. Let’s use the same DGP.

``` r
p <- 10
n <- 1e5
beta <- rnorm(p + 1)

x <- cbind(1, matrix(rnorm(p * n), ncol = p, nrow = n))
y <- x %*% beta + exp(rnorm(n, sd= 2))

start = proc.time()
fit <- fit_fqr(X = x, y = y, tau = 0.5, se = F)
end = proc.time()
end - start
#>    user  system elapsed 
#>   0.832   0.069   0.231
```

``` r
start = proc.time()
fit_pfn <- quantreg::rq.fit.pfn(x = x, y = y, tau = 0.5)
end = proc.time()
end - start
#>    user  system elapsed 
#>   1.217   0.263   1.426
```

The coefficients match out to the 4th or 5th decimal place:

``` r
max(abs(fit$coefficients - fit_pfn$coefficients))
#> [1] 8.977448e-05
```

``` r
min(abs(fit$coefficients - fit_pfn$coefficients))
#> [1] 2.726844e-06
```

### Small Problems

It can also be faster for small problems, and with conservative
tolerance parameters will come extremely close to the default `quantreg`
outputs. Here’s an example:

``` r
# simulate some data, 101 x 10,000,000
p <- 3
n <- 10000
beta <- rnorm(p)
x <- cbind(matrix(rnorm(p * n), ncol = p, nrow = n))
y <- 10 + x %*% beta + exp(rnorm(n, sd = 2))

microbenchmark::microbenchmark(
  fqr_fit <- fqr(y ~ ., se = F, beta_tol = 0, check_tol = 0,
                 data = data.frame(y = y, x)),
  br_fit <- quantreg::rq(y ~ ., tau = 0.5, 
                    data = data.frame(y = y, x), method = "br"),
                    times = 100
)
#> Unit: milliseconds
#>                                                                                          expr
#>  fqr_fit <- fqr(y ~ ., se = F, beta_tol = 0, check_tol = 0, data = data.frame(y = y,      x))
#>     br_fit <- quantreg::rq(y ~ ., tau = 0.5, data = data.frame(y = y,      x), method = "br")
#>       min       lq     mean   median       uq      max neval
#>  19.05055 19.98623 20.92135 20.65108 21.24894 26.97846   100
#>  67.57361 68.66483 70.66433 69.19325 70.23377 83.04198   100
```

The coefficients match out to 4 decimal places:

``` r
fqr_fit$coefficients - br_fit$coefficients
#>               [,1]
#> [1,]  1.053418e-05
#> [2,] -4.118319e-05
#> [3,] -3.619511e-05
#> [4,]  1.184613e-05
```

And the check loss is nearly identical:

``` r
check <- function (x, tau = 0.5) {
  sum(x * (tau - (x < 0)))
}

check(fqr_fit$residuals) -  check(br_fit$residuals)
#> [1] 4.954552e-05
```

Still, though, the speed gains are most noticeable once N and P are
“medium” or larger (e.g. for N &lt; 300, probably just use quantreg).
