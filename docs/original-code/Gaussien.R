#Paramètres initiaux
mu <- 0
sigma2_0 <- 1
alpha <- 0.05
n_sim  <- 400
n_values <- c(50, 100, 200, 500, 1000)


simulate_test <- function(n, sigma2_actual, sigma2_0, alpha, n_sim) {
  rejets <- 0
  for(i in 1:n_sim) {
    echantillon <- rnorm(n, mean = mu, sd = sqrt(sigma2_actual))
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

empirique <- sapply(n_values, function(n) simulate_test(n, sigma2_0, sigma2_0, alpha, n_sim))
resultats_H0 <- data.frame(Taille_n = n_values, Rejet_Pourcentage = empirique)
print(resultats_H0)


n_graph <- 1000
T_valeurs <- numeric(n_sim)

for(i in 1:n_sim) {
  echantillon <- rnorm(n_graph, mean = mu, sd = sqrt(sigma2_0))
  T_valeurs[i] <- (n_graph - 1) * var(echantillon) / sigma2_0
}

hist(T_valeurs, breaks = 30, probability = TRUE,
     main = "Distribution de la statistique de test sous H0",
     xlab = "Valeurs de la statistique T", col = "lightblue"
     )

curve(dchisq(x, df = n_graph - 1), add = TRUE, col = "red", lwd=2)
legend("topleft", legend = c("Histogramme empirique", "Densité du Chi-deux"), fill = c("lightblue", "red"))


epsilons <- c(0.05, 0.2, 0.5)
matrice <- matrix(0, nrow = length(epsilons), ncol = length(n_values))
rownames(matrice) <- paste("Epsilon = ", epsilons)
colnames(matrice) <- paste("n =", n_values)

for(i in 1:length(epsilons)) {
  for (j in 1: length(n_values)) {
    sigma2_actual <- sigma2_0 + epsilons[i]
    matrice[i, j] <- simulate_test(n_values[j], sigma2_actual, sigma2_0, alpha, n_sim)
  }
}

print(matrice)

