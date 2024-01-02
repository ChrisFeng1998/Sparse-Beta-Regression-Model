
# Test out buster ---------------------------------------------------------

library(glmnet)
library(matrixStats)
library(tidyverse)


# n = 300 -----------------------------------------------------------------

# Number nodes
n <- 300
p <- 20

# Number Monte Carlo iterations (uneven so that median exists)
M <- 20

xi <- 1.5
mu_dagger <- 1
gamma_star <- c(1.5, 1.2, 0.8,rep(1,p-3))
mu <- -xi*log(n) + mu_dagger

# Helper functions --------------------------------------------------------

p_ij <- function(mu, i, j){
  # must use j,i here, because we use lower.tri later on
  Z_ij <- sapply(Z_matrix, function(M) M[j,i])%*%gamma_star
  1/(1 + exp(-(mu + Z_ij)))
}

# Environment -------------------------------------------------------------

setwd("/storage")
SGE_ID <- as.numeric(Sys.getenv("SLURM_JOB_ID"))
set.seed(SGE_ID)

# Monte Carlo Start -------------------------------------------------------

start_time <- Sys.time()
for (m in 1:M) {
  print(m)
  
  
  Z_matrix <- lapply(1:p, function(i) matrix(nrow = n, ncol = n, data = rbeta(n*n, 2, 2) - 1/2))
  Z <- do.call(cbind, lapply(Z_matrix, function(M) M[lower.tri(M)]))

  probabilities <- matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in i:n) {
      probabilities[j,i] <- p_ij(mu, i, j)
    }
  }
  
  adj_matrix <- matrix(nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      #so far this is only the upper half of the adj matrix.
      adj_matrix[j,i] <- sample(0:1, 1, prob=c(1-probabilities[j,i],probabilities[j,i]))
      adj_matrix[i,j] <- adj_matrix[j,i]
    }
  }
  # Record Network Summaries ------------------------------------------------
  degrees <- rowSums(adj_matrix, na.rm = TRUE)
  network_summaries <- c(length(unique(degrees)),
                         sum(degrees)/2,
                         sum(degrees)/2/choose(n,2),
                         summary(probabilities[lower.tri(probabilities)]))
  rm(probabilities)
  # Create fit --------------------------------------------------------------
  
  # Values Y; I have to take lower.tri to get the right order of entries
  Y <- adj_matrix[lower.tri(adj_matrix)]
  # Intercept mu fitted automatically
  fit <- glm(Y~., family = binomial(link='logit'), data = data.frame(Y,Z))
  # summary(fit)
  # abs(fit$coefficients - c(mu, gamma_star))%>%sum()

  # Update Network summaries ------------------------------------------------
  
  gamma_l1_error <- sum(abs(fit$coefficients[2:(p+1)] - gamma_star))
  gamma_l2_error <- sqrt(sum((fit$coefficients[2:(p+1)] - gamma_star)^2))
  # notice that this is also the error for mu_dagger when xi is known
  mu_error <- fit$coefficients[1] - mu
  
  
  # log likelihood
  fitted_probabilities <- predict(fit, newx = Z, type = "response")
  # This is the POSITIVE log likelihood for each fitted model.

  
  # Inference
  hat_gamma <- fit$coefficients[2:(p+1)]
  D <- cbind(1, Z)
  sample_gram_matrix <- 1/choose(n,2)*t(D)%*%(fitted_probabilities*(1 - fitted_probabilities)*D)
  root_Theta_hat_diag <- sample_gram_matrix%>%solve()%>%diag()%>%sqrt()
  root_Theta_hat_diag <- root_Theta_hat_diag[2:(p+1)]
  CI_gamma <- rbind(hat_gamma + qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)),
                    hat_gamma - qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)))
  CI_lengths <- CI_gamma[1,] - CI_gamma[2,]
  coverage <- (gamma_star <= CI_gamma[1,]) & (gamma_star >= CI_gamma[2,])
  
  standardized_gamma_matrix <- sqrt(choose(n,2))*(hat_gamma - gamma_star)/root_Theta_hat_diag
  
  
  readr::write_csv(as.data.frame(t(c(standardized_gamma_matrix, CI_lengths, coverage))), 
                   path = paste0(SGE_ID, "_", n, "_gamma_inference.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(mu_error)), path = paste0(SGE_ID, "_", n, "mu_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(gamma_l1_error)), path = paste0(SGE_ID, "_", n, "_gamma_l1_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(gamma_l2_error)), path = paste0(SGE_ID, "_", n, "_gamma_l2_error.csv"), append = TRUE)
  
  

  
  # Other summaries
  
  readr::write_csv(as.data.frame(t(network_summaries)), 
                   path = paste0(SGE_ID, "_", n, "_network_summaries.csv"), append = TRUE)

}
end_time <- Sys.time()
print(end_time - start_time)


# Clear workspace ---------------------------------------------------------

rm(list = ls())


# n = 500 -----------------------------------------------------------------

# Number nodes
n <- 500
p <- 20

# Number Monte Carlo iterations (uneven so that median exists)
M <- 20

xi <- 1.5
mu_dagger <- 1
gamma_star <- c(1.5, 1.2,0.8,rep(1,p-3))
mu <- -xi*log(n) + mu_dagger

# Helper functions --------------------------------------------------------

p_ij <- function(mu, i, j){
  # must use j,i here, because we use lower.tri later on
  Z_ij <- sapply(Z_matrix, function(M) M[j,i])%*%gamma_star
  1/(1 + exp(-(mu + Z_ij)))
}

# Environment -------------------------------------------------------------

SGE_ID <- as.numeric(Sys.getenv("SLURM_JOB_ID"))
set.seed(SGE_ID+100)

# Monte Carlo Start -------------------------------------------------------

start_time <- Sys.time()
for (m in 1:M) {
  print(m)
  
  
  Z_matrix <- lapply(1:p, function(i) matrix(nrow = n, ncol = n, data = rbeta(n*n, 2, 2) - 1/2))
  Z <- do.call(cbind, lapply(Z_matrix, function(M) M[lower.tri(M)]))
  
  probabilities <- matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in i:n) {
      probabilities[j,i] <- p_ij(mu, i, j)
    }
  }
  
  adj_matrix <- matrix(nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      #so far this is only the upper half of the adj matrix.
      adj_matrix[j,i] <- sample(0:1, 1, prob=c(1-probabilities[j,i],probabilities[j,i]))
      adj_matrix[i,j] <- adj_matrix[j,i]
    }
  }
  # Record Network Summaries ------------------------------------------------
  degrees <- rowSums(adj_matrix, na.rm = TRUE)
  network_summaries <- c(length(unique(degrees)),
                         sum(degrees)/2,
                         sum(degrees)/2/choose(n,2),
                         summary(probabilities[lower.tri(probabilities)]))
  rm(probabilities)
  
  # Create fit --------------------------------------------------------------
  
  # Values Y; I have to take lower.tri to get the right order of entries
  Y <- adj_matrix[lower.tri(adj_matrix)]
  # Intercept mu fitted automatically
  fit <- glm(Y~., family = binomial(link='logit'), data = data.frame(Y,Z))
  # summary(fit)
  # abs(fit$coefficients - c(mu, gamma_star))%>%sum()
  
  # Update Network summaries ------------------------------------------------
  
  gamma_l1_error <- sum(abs(fit$coefficients[2:(p+1)] - gamma_star))
  gamma_l2_error <- sqrt(sum((fit$coefficients[2:(p+1)] - gamma_star)^2))
  # notice that this is also the error for mu_dagger when xi is known
  mu_error <- fit$coefficients[1] - mu
  
  
  # log likelihood
  fitted_probabilities <- predict(fit, newx = Z, type = "response")
  # This is the POSITIVE log likelihood for each fitted model.
  
  
  # Inference
  hat_gamma <- fit$coefficients[2:(p+1)]
  D <- cbind(1, Z)
  sample_gram_matrix <- 1/choose(n,2)*t(D)%*%(fitted_probabilities*(1 - fitted_probabilities)*D)
  root_Theta_hat_diag <- sample_gram_matrix%>%solve()%>%diag()%>%sqrt()
  root_Theta_hat_diag <- root_Theta_hat_diag[2:(p+1)]
  CI_gamma <- rbind(hat_gamma + qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)),
                    hat_gamma - qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)))
  CI_lengths <- CI_gamma[1,] - CI_gamma[2,]
  coverage <- (gamma_star <= CI_gamma[1,]) & (gamma_star >= CI_gamma[2,])
  
  standardized_gamma_matrix <- sqrt(choose(n,2))*(hat_gamma - gamma_star)/root_Theta_hat_diag
  
  
  readr::write_csv(as.data.frame(t(c(standardized_gamma_matrix, CI_lengths, coverage))), 
                   path = paste0(SGE_ID, "_", n, "_gamma_inference.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(mu_error)), path = paste0(SGE_ID, "_", n, "mu_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(gamma_l1_error)), path = paste0(SGE_ID, "_", n, "_gamma_l1_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(gamma_l2_error)), path = paste0(SGE_ID, "_", n, "_gamma_l2_error.csv"), append = TRUE)
  
  
  
  
  # Other summaries
  
  readr::write_csv(as.data.frame(t(network_summaries)), 
                   path = paste0(SGE_ID, "_", n, "_network_summaries.csv"), append = TRUE)
  
}
end_time <- Sys.time()
print(end_time - start_time)



# n=800 -------------------------------------------------------------------

# Number nodes
n <- 800
p <- 20

# Number Monte Carlo iterations (uneven so that median exists)
M <- 20

xi <- 1.5
mu_dagger <- 1
gamma_star <- c(1.5, 1.2,0.8,rep(1,p-3))
mu <- -xi*log(n) + mu_dagger

# Helper functions --------------------------------------------------------

p_ij <- function(mu, i, j){
  # must use j,i here, because we use lower.tri later on
  Z_ij <- sapply(Z_matrix, function(M) M[j,i])%*%gamma_star
  1/(1 + exp(-(mu + Z_ij)))
}

# Environment -------------------------------------------------------------

SGE_ID <- as.numeric(Sys.getenv("SLURM_JOB_ID"))
set.seed(SGE_ID + 200)

# Monte Carlo Start -------------------------------------------------------

start_time <- Sys.time()
for (m in 1:M) {
  print(m)
  
  
  Z_matrix <- lapply(1:p, function(i) matrix(nrow = n, ncol = n, data = rbeta(n*n, 2, 2) - 1/2))
  Z <- do.call(cbind, lapply(Z_matrix, function(M) M[lower.tri(M)]))
  
  probabilities <- matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in i:n) {
      probabilities[j,i] <- p_ij(mu, i, j)
    }
  }
  
  adj_matrix <- matrix(nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      #so far this is only the upper half of the adj matrix.
      adj_matrix[j,i] <- sample(0:1, 1, prob=c(1-probabilities[j,i],probabilities[j,i]))
      adj_matrix[i,j] <- adj_matrix[j,i]
    }
  }
  # Record Network Summaries ------------------------------------------------
  degrees <- rowSums(adj_matrix, na.rm = TRUE)
  network_summaries <- c(length(unique(degrees)),
                         sum(degrees)/2,
                         sum(degrees)/2/choose(n,2),
                         summary(probabilities[lower.tri(probabilities)]))
  rm(probabilities)
  
  # Create fit --------------------------------------------------------------
  
  # Values Y; I have to take lower.tri to get the right order of entries
  Y <- adj_matrix[lower.tri(adj_matrix)]
  # Intercept mu fitted automatically
  fit <- glm(Y~., family = binomial(link='logit'), data = data.frame(Y,Z))
  # summary(fit)
  # abs(fit$coefficients - c(mu, gamma_star))%>%sum()
  
  # Update Network summaries ------------------------------------------------
  
  gamma_l1_error <- sum(abs(fit$coefficients[2:(p+1)] - gamma_star))
  gamma_l2_error <- sqrt(sum((fit$coefficients[2:(p+1)] - gamma_star)^2))
  # notice that this is also the error for mu_dagger when xi is known
  mu_error <- fit$coefficients[1] - mu
  
  
  # log likelihood
  fitted_probabilities <- predict(fit, newx = Z, type = "response")
  # This is the POSITIVE log likelihood for each fitted model.
  
  
  # Inference
  hat_gamma <- fit$coefficients[2:(p+1)]
  D <- cbind(1, Z)
  sample_gram_matrix <- 1/choose(n,2)*t(D)%*%(fitted_probabilities*(1 - fitted_probabilities)*D)
  root_Theta_hat_diag <- sample_gram_matrix%>%solve()%>%diag()%>%sqrt()
  root_Theta_hat_diag <- root_Theta_hat_diag[2:(p+1)]
  CI_gamma <- rbind(hat_gamma + qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)),
                    hat_gamma - qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)))
  CI_lengths <- CI_gamma[1,] - CI_gamma[2,]
  coverage <- (gamma_star <= CI_gamma[1,]) & (gamma_star >= CI_gamma[2,])
  
  standardized_gamma_matrix <- sqrt(choose(n,2))*(hat_gamma - gamma_star)/root_Theta_hat_diag
  
  
  readr::write_csv(as.data.frame(t(c(standardized_gamma_matrix, CI_lengths, coverage))), 
                   path = paste0(SGE_ID, "_", n, "_gamma_inference.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(mu_error)), path = paste0(SGE_ID, "_", n, "mu_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(gamma_l1_error)), path = paste0(SGE_ID, "_", n, "_gamma_l1_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(gamma_l2_error)), path = paste0(SGE_ID, "_", n, "_gamma_l2_error.csv"), append = TRUE)
  
  
  
  
  # Other summaries
  
  readr::write_csv(as.data.frame(t(network_summaries)), 
                   path = paste0(SGE_ID, "_", n, "_network_summaries.csv"), append = TRUE)
  
}
end_time <- Sys.time()
print(end_time - start_time)


# n = 1000 ----------------------------------------------------------------

# Number nodes
n <- 1000
p <- 20

# Number Monte Carlo iterations (uneven so that median exists)
M <- 20

xi <- 1.5
mu_dagger <- 1
gamma_star <- c(1.5, 1.2, 0.8,rep(1,p-3))
mu <- -xi*log(n) + mu_dagger

# Helper functions --------------------------------------------------------

p_ij <- function(mu, i, j){
  # must use j,i here, because we use lower.tri later on
  Z_ij <- sapply(Z_matrix, function(M) M[j,i])%*%gamma_star
  1/(1 + exp(-(mu + Z_ij)))
}

# Environment -------------------------------------------------------------

SGE_ID <- as.numeric(Sys.getenv("SLURM_JOB_ID"))
set.seed(SGE_ID + 300)

# Monte Carlo Start -------------------------------------------------------

start_time <- Sys.time()
for (m in 1:M) {
  print(m)
  
  
  Z_matrix <- lapply(1:p, function(i) matrix(nrow = n, ncol = n, data = rbeta(n*n, 2, 2) - 1/2))
  Z <- do.call(cbind, lapply(Z_matrix, function(M) M[lower.tri(M)]))
  
  probabilities <- matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in i:n) {
      probabilities[j,i] <- p_ij(mu, i, j)
    }
  }
  
  adj_matrix <- matrix(nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      #so far this is only the upper half of the adj matrix.
      adj_matrix[j,i] <- sample(0:1, 1, prob=c(1-probabilities[j,i],probabilities[j,i]))
      adj_matrix[i,j] <- adj_matrix[j,i]
    }
  }
  # Record Network Summaries ------------------------------------------------
  degrees <- rowSums(adj_matrix, na.rm = TRUE)
  network_summaries <- c(length(unique(degrees)),
                         sum(degrees)/2,
                         sum(degrees)/2/choose(n,2),
                         summary(probabilities[lower.tri(probabilities)]))
  rm(probabilities)
  
  # Create fit --------------------------------------------------------------
  
  # Values Y; I have to take lower.tri to get the right order of entries
  Y <- adj_matrix[lower.tri(adj_matrix)]
  # Intercept mu fitted automatically
  fit <- glm(Y~., family = binomial(link='logit'), data = data.frame(Y,Z))
  # summary(fit)
  # abs(fit$coefficients - c(mu, gamma_star))%>%sum()
  
  # Update Network summaries ------------------------------------------------
  
  gamma_l1_error <- sum(abs(fit$coefficients[2:(p+1)] - gamma_star))
  gamma_l2_error <- sqrt(sum((fit$coefficients[2:(p+1)] - gamma_star)^2))
  # notice that this is also the error for mu_dagger when xi is known
  mu_error <- fit$coefficients[1] - mu
  
  
  # log likelihood
  fitted_probabilities <- predict(fit, newx = Z, type = "response")

  
  # Inference
  hat_gamma <- fit$coefficients[2:(p+1)]
  D <- cbind(1, Z)
  sample_gram_matrix <- 1/choose(n,2)*t(D)%*%(fitted_probabilities*(1 - fitted_probabilities)*D)
  root_Theta_hat_diag <- sample_gram_matrix%>%solve()%>%diag()%>%sqrt()
  root_Theta_hat_diag <- root_Theta_hat_diag[2:(p+1)]
  CI_gamma <- rbind(hat_gamma + qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)),
                    hat_gamma - qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)))
  CI_lengths <- CI_gamma[1,] - CI_gamma[2,]
  coverage <- (gamma_star <= CI_gamma[1,]) & (gamma_star >= CI_gamma[2,])
  
  standardized_gamma_matrix <- sqrt(choose(n,2))*(hat_gamma - gamma_star)/root_Theta_hat_diag
  
  
  readr::write_csv(as.data.frame(t(c(standardized_gamma_matrix, CI_lengths, coverage))), 
                   path = paste0(SGE_ID, "_", n, "_gamma_inference.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(mu_error)), path = paste0(SGE_ID, "_", n, "mu_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(gamma_l1_error)), path = paste0(SGE_ID, "_", n, "_gamma_l1_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(gamma_l2_error)), path = paste0(SGE_ID, "_", n, "_gamma_l2_error.csv"), append = TRUE)
  
  
  
  
  # Other summaries
  
  readr::write_csv(as.data.frame(t(network_summaries)), 
                   path = paste0(SGE_ID, "_", n, "_network_summaries.csv"), append = TRUE)
  
}
end_time <- Sys.time()
print(end_time - start_time)

