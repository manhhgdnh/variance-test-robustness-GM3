source("R/variance_test.R")

set.seed(2026)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

sigma2_0 <- 1
alpha <- 0.05
n_sim <- 400
n_values <- c(50, 100, 200, 500, 1000)
epsilons <- c(0.05, 0.2, 0.5)
rhos <- c(0.1, 0.4, 0.8)

# H0 calibration -----------------------------------------------------------
null_rejection <- sapply(
  n_values,
  function(n) simulate_gaussian_rejection(n, sigma2_0, sigma2_0, alpha, n_sim)
)
null_results <- data.frame(n = n_values, rejection_percent = null_rejection)
write.csv(null_results, "results/tables/null_calibration_simulated.csv", row.names = FALSE)

png("results/figures/null_calibration_simulated.png", width = 1000, height = 650)
plot(n_values, null_rejection, type = "b", log = "x", ylim = c(0, max(10, null_rejection)),
     xlab = "Sample size n", ylab = "Rejection rate (%)",
     main = "Empirical type-I error under H0")
abline(h = 100 * alpha, lty = 2)
dev.off()

# Chi-square fit -----------------------------------------------------------
n_graph <- 1000
t_values <- replicate(n_sim, {
  x <- rnorm(n_graph, mean = 0, sd = sqrt(sigma2_0))
  variance_test_statistic(x, sigma2_0)
})

png("results/figures/chi_square_fit_simulated.png", width = 1000, height = 650)
hist(t_values, breaks = 30, probability = TRUE,
     main = "Test statistic under H0",
     xlab = "T = (n - 1) S^2 / sigma0^2")
curve(dchisq(x, df = n_graph - 1), add = TRUE, lwd = 2)
legend("topleft", legend = c("Empirical histogram", "Chi-square density"),
       lty = c(NA, 1), pch = c(15, NA), bty = "n")
dev.off()

# Power --------------------------------------------------------------------
power_matrix <- sapply(epsilons, function(epsilon) {
  sapply(n_values, function(n) {
    simulate_gaussian_rejection(n, sigma2_0 + epsilon, sigma2_0, alpha, n_sim)
  })
})
colnames(power_matrix) <- paste0("epsilon_", epsilons)
power_results <- data.frame(n = n_values, power_matrix)
write.csv(power_results, "results/tables/power_simulated.csv", row.names = FALSE)

png("results/figures/power_curves_simulated.png", width = 1000, height = 650)
matplot(n_values, power_matrix, type = "b", log = "x", ylim = c(0, 100),
        xlab = "Sample size n", ylab = "Rejection rate / power (%)",
        main = "Empirical power under H1")
legend("bottomright", legend = paste("epsilon =", epsilons), lty = 1:3, pch = 1:3, bty = "n")
dev.off()

# Robustness: exponential --------------------------------------------------
exp_rejection <- sapply(
  n_values,
  function(n) simulate_exponential_rejection(n, sigma2_0, alpha, n_sim)
)
exp_results <- data.frame(n = n_values, rejection_percent = exp_rejection)
write.csv(exp_results, "results/tables/exponential_simulated.csv", row.names = FALSE)

png("results/figures/exponential_robustness_simulated.png", width = 1000, height = 650)
plot(n_values, exp_rejection, type = "b", log = "x", ylim = c(0, max(40, exp_rejection)),
     xlab = "Sample size n", ylab = "Rejection rate (%)",
     main = "Robustness to non-Gaussian data: Exponential distribution")
abline(h = 100 * alpha, lty = 2)
dev.off()

# Robustness: AR(1) --------------------------------------------------------
ar1_matrix <- sapply(rhos, function(rho) {
  sapply(n_values, function(n) {
    simulate_ar1_rejection(n, rho, sigma2_0, alpha, n_sim)
  })
})
colnames(ar1_matrix) <- paste0("rho_", rhos)
ar1_results <- data.frame(n = n_values, ar1_matrix)
write.csv(ar1_results, "results/tables/ar1_simulated.csv", row.names = FALSE)

png("results/figures/ar1_robustness_simulated.png", width = 1000, height = 650)
matplot(n_values, ar1_matrix, type = "b", log = "x", ylim = c(0, max(45, ar1_matrix)),
        xlab = "Sample size n", ylab = "Rejection rate (%)",
        main = "Robustness to dependence: AR(1)")
abline(h = 100 * alpha, lty = 2)
legend("topleft", legend = paste("rho =", rhos), lty = 1:3, pch = 1:3, bty = "n")
dev.off()

cat("Analysis complete. Tables and figures are in results/.\n")
