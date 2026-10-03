# Core functions for the variance conformity test experiments.

variance_test_statistic <- function(x, sigma2_0 = 1) {
  n <- length(x)
  (n - 1) * var(x) / sigma2_0
}

variance_test_reject <- function(x, sigma2_0 = 1, alpha = 0.05) {
  n <- length(x)
  statistic <- variance_test_statistic(x, sigma2_0)
  lower <- qchisq(alpha / 2, df = n - 1)
  upper <- qchisq(1 - alpha / 2, df = n - 1)
  statistic < lower || statistic > upper
}

simulate_gaussian_rejection <- function(n, sigma2_actual, sigma2_0 = 1,
                                        alpha = 0.05, n_sim = 400, mu = 0) {
  rejected <- replicate(n_sim, {
    x <- rnorm(n, mean = mu, sd = sqrt(sigma2_actual))
    variance_test_reject(x, sigma2_0, alpha)
  })
  100 * mean(rejected)
}

simulate_exponential_rejection <- function(n, sigma2_0 = 1,
                                           alpha = 0.05, n_sim = 400) {
  # For Exp(lambda), Var(X) = 1 / lambda^2.
  lambda <- 1 / sqrt(sigma2_0)
  rejected <- replicate(n_sim, {
    x <- rexp(n, rate = lambda)
    variance_test_reject(x, sigma2_0, alpha)
  })
  100 * mean(rejected)
}

simulate_ar1_rejection <- function(n, rho, sigma2_0 = 1,
                                   alpha = 0.05, n_sim = 400) {
  sigma2_eps <- sigma2_0 * (1 - rho^2)
  rejected <- replicate(n_sim, {
    x <- numeric(n)
    x[1] <- rnorm(1, mean = 0, sd = sqrt(sigma2_0))
    for (i in 2:n) {
      eps <- rnorm(1, mean = 0, sd = sqrt(sigma2_eps))
      x[i] <- rho * x[i - 1] + eps
    }
    variance_test_reject(x, sigma2_0, alpha)
  })
  100 * mean(rejected)
}
