mu <- 0
sigma2_0 <- 1
alpha <- 0.05
n_sim <- 400
n_values <- c(50, 100, 200, 500, 1000)

simulate_exponentielle <- function(n, sigma2_0, alpha, n_sim) {
  rejets <- 0
  
  # Pour avoir Var = sigma2_0, il faut lambda = 1 / sqrt(sigma2_0)
  lambda <- 1 / sqrt(sigma2_0) 
  
  for(i in 1:n_sim) {
    echantillon <- rexp(n, rate = lambda)
    S2 <- var(echantillon)
    T_stat <- (n - 1) * S2 / sigma2_0
    q_inf <- qchisq(alpha / 2, df = n - 1)
    q_sup <- qchisq(1 - alpha / 2, df = n - 1)
    
    if(T_stat < q_inf | T_stat > q_sup) {
      rejets <- rejets + 1
    }
  }
  return((rejets / n_sim) * 100)
}

empirique <- sapply(n_values, function(n) simulate_exponentielle(n, sigma2_0, alpha, n_sim))
resultats_exponentielle <- data.frame(Taille_n = n_values, Rejet_Pourcentage = empirique)
print(resultats_exponentielle)

simulate_dependance <- function(n, rho, sigma2_0, alpha, n_sim) {
  rejets <- 0
  sigma2_eps <- sigma2_0 * (1 - rho^2)
  
  for(i in 1:n_sim) {
    X <- numeric(n)
    X[1] <- rnorm(1, mean = 0, sd = sqrt(sigma2_0))
    
    for(j in 2:n) {
      eps <- rnorm(1, mean = 0, sd = sqrt(sigma2_eps))
      X[j] <- rho * X[j-1] + eps
    }
    
    S2 <- var(X)
    T_stat <- (n - 1) * S2 / sigma2_0
    
    q_inf <- qchisq(alpha / 2, df = n - 1)
    q_sup <- qchisq(1 - alpha / 2, df = n - 1)
    
    if(T_stat < q_inf | T_stat > q_sup) {
      rejets <- rejets + 1
    }
  } 
  return((rejets / n_sim) * 100)
}

valeurs_rho <- c(0.1, 0.4, 0.8)

matrice <- matrix(0, nrow = length(valeurs_rho), ncol = length(n_values))
rownames(matrice) <- paste("Rho =", valeurs_rho) 
colnames(matrice) <- paste("n =", n_values)

for(i in 1:length(valeurs_rho)) {
  for(j in 1:length(n_values)) {
    matrice[i, j] <- simulate_dependance(n_values[j], valeurs_rho[i], sigma2_0, alpha, n_sim)
  }
}

print(matrice)