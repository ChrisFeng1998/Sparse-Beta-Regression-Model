

# Run model simulation ---------------------------------------------------------

library(glmnet)
library(matrixStats)
library(tidyverse)

# Number nodes
n <- 300
p <- 2
s <- 7 #300
#s <- 9 # 500
#s <- 11 #800
#s <- 12 # 1000

M <- 100

betas <- c(1.2,0.8,rep(1,s-2),rep(0,n-s))
mu <- -0.5*log(log(n))
gamma_star <- c(1,0.8)

# For fixed lambda
# Penaltiy
a <- sqrt(2*log(2*(n+p+1))/choose(n,2))
t <- 2
lambda_0 <- 4*a + 2*sqrt(2*t/choose(n,2) * (1+ sqrt(2*n) * a)) + t*sqrt(2*n)/(3*choose(n,2))
lambda_bar <- 8*lambda_0
lambda <- sqrt(2/n)*lambda_0

rm(a,t,lambda_0, lambda_bar)

# GIC parameter
a_n <- log(log(choose(n,2)))*log(n+1)

# Helper functions --------------------------------------------------------

p_ij <- function(mu, beta_i, beta_j, i, j){
  Z_ij <- c(Z_1_matrix[j,i], Z_2_matrix[j,i])%*%gamma_star
  1/(1 + exp(-(beta_i + beta_j + mu + Z_ij)))
}

log_like <- function(adj, probs){
  return(sum(adj%*%log(probs) + (1-adj)%*%log(1 - probs)))
}

actives_as_inactive <- function(vec){
  return(sum(which(vec == 0) <= s))
}
inactives_as_active <- function(vec){
  return(sum(which(vec > 0) > s))
}

# Create sparse design matrix ---------------------------------------------

# Sparse design matrix (same for each iteration, so created outside of loop)

# row indices: we take all rows (no empty rows) and two entries per row
rows <- rep(1:choose(n,2), each = 2)
# column indices
a <- 2:n
b <- rep(n, n-1)
c <- NULL
for (i in 1:(n-1)) {
  c <- c(c, seq(a[i],b[i]))
}
cols <- c(rbind(rep(1:(n-1), (n-1):1),c))

X_deterministic <- sparseMatrix(i = rows, j = cols, x = 1)
rm(a,b,c, rows, cols, i)


# Environment -------------------------------------------------------------

setwd("/storage")
SGE_ID <-  as.numeric(Sys.getenv('SLURM_JOB_ID'))
set.seed(SGE_ID)

# Monte Carlo Start -------------------------------------------------------

start_time <- Sys.time()
for (m in 1:M) {
  print(m)
  
  # Similarity matrix for first covariate
  Z_1_matrix <- matrix(nrow = n, ncol = n, data = rbeta(n*n, 2, 2) - 1/2)
  # Similarity values for first covariate in vector form
  Z_1 <- Z_1_matrix[lower.tri(Z_1_matrix)]
  # Similarity matrix for second covariate
  Z_2_matrix <- matrix(nrow = n, ncol = n, data = rbeta(n*n, 2, 2) - 1/2)
  # Similarity values for second covariate in vector form
  Z_2 <- Z_2_matrix[lower.tri(Z_2_matrix)]
  
  probabilities <- matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in i:n) {
      probabilities[j,i] <- p_ij(mu, betas[i], betas[j], i, j)
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
  rm(probabilities, Z_1_matrix, Z_2_matrix)
  # Record Network Summaries ------------------------------------------------
  network_summaries <- numeric()
  degrees_of_freedom <- numeric()
  log_lambda_matrix <- numeric()
  gamma_l1_error_matrix <- numeric()
  gamma_l2_error_matrix <- numeric()
  relative_l2_error_matrix <- numeric()
  BICs <- numeric()
  AICs <- numeric()
  GICs <- numeric()
  IC_matrix <- numeric()
  misclassification_matrix <- numeric()
  
  
  degrees <- rowSums(adj_matrix, na.rm = TRUE)
  network_summaries[1:3] <- c(s,
                              length(unique(degrees)),
                              sum(degrees)/2)
  
  # Create fit --------------------------------------------------------------
  
  
  # Design matrix
  X <- cbind(X_deterministic, Z_1, Z_2)
  
  # Values Y; I have to take lower.tri to get the right order of entries
  Y <- adj_matrix[lower.tri(adj_matrix)]
  # we can leave the default standardize = TRUE option, results are essentially the same
  fit <- glmnet::glmnet(X,Y, family = "binomial", alpha = 1, penalty.factor = c(rep(1,n),0,0),
                        lower.limits = c(rep(0,n), -Inf, -Inf))
  
  
  
  # Update Network summaries ------------------------------------------------
  
  # We only want to know number of active BETAS; also remember:
  # unpenalized variables, i.e. gamma, will always be active
  degrees_of_freedom <- (fit$df - p)
  
  log_lambda_matrix <- log(fit$lambda)
  covariate_error <- abs(as.matrix(fit$beta[(n+1):(n+p),]) - gamma_star)
  # Gamma errors
  gamma_l1_error_matrix <- colSums(covariate_error)
  gamma_l2_error_matrix <- sqrt(colSums(covariate_error^2))
  
  # Error Statistics
  # beta
  beta_abs_errors <- (abs(as.matrix(fit$beta[1:n,]) - betas))
  # l_1 errors
  beta_mae_errors <- colSums(beta_abs_errors)/n
  # l_2 errors: square all entries, build column sums, divide by n+1
  beta_l2_errors <- sqrt(colSums(beta_abs_errors^2))
  beta_mse_errors <- beta_l2_errors/n
  
  # mu
  mu_error <- fit$a0 - mu
  
  
  # log likelihood
  fitted_probabilities <- predict(fit, newx = X, type = "response")
  # This is the POSITIVE log likelihood for each fitted model.
  log_like_vector <- apply(fitted_probabilities, 
                           MARGIN = 2, 
                           function(i) log_like(adj = adj_matrix[lower.tri(adj_matrix)], probs = i))
  
  bic <- fit$df*log(choose(n,2)) - 2*log_like_vector
  aic <- 2*(fit$df - log_like_vector)
  gic <- fit$df*a_n - 2*log_like_vector
  BICs <- bic
  AICs <- aic
  GICs <- gic
  bic_best_pos <- which.min(bic)
  aic_best_pos <- which.min(aic)
  gic_best_pos <- which.min(gic)
  
  
  IC_matrix <- c(bic[bic_best_pos], # optimal bic value
                 bic_best_pos, # position of optimal bic value
                 fit$df[bic_best_pos] - p, # how many degrees of freedom (=active indices)?
                 sum(which(coef(fit, s = fit$lambda[bic_best_pos])[2:(n+1)] == 0) <= s), # how many actives were misclassified as inactive?
                 sum(which(coef(fit, s = fit$lambda[bic_best_pos])[2:(n+1)] > 0) > s), # how many inactives as active?
                 sum(coef(fit, s = fit$lambda[bic_best_pos])[(s+2):(n+1)])/(n-s), # mean l1 norm of misclas inactive betas
                 beta_mae_errors[bic_best_pos], # mae
                 beta_mse_errors[bic_best_pos], #mse
                 aic[aic_best_pos],
                 aic_best_pos,
                 fit$df[aic_best_pos] - p,
                 sum(which(coef(fit, s = fit$lambda[aic_best_pos])[2:(n+1)] == 0) <= s),
                 sum(which(coef(fit, s = fit$lambda[aic_best_pos])[2:(n+1)] > 0) > s),
                 sum(coef(fit, s = fit$lambda[aic_best_pos])[(s+2):(n+1)])/(n-s),
                 beta_mae_errors[aic_best_pos], 
                 beta_mse_errors[aic_best_pos], 
                 gic[gic_best_pos],
                 gic_best_pos,
                 fit$df[gic_best_pos] - p,
                 sum(which(coef(fit, s = fit$lambda[gic_best_pos])[2:(n+1)] == 0) <= s),
                 sum(which(coef(fit, s = fit$lambda[gic_best_pos])[2:(n+1)] > 0) > s),
                 sum(coef(fit, s = fit$lambda[gic_best_pos])[(s+2):(n+1)])/(n-s),
                 beta_mae_errors[gic_best_pos], 
                 beta_mse_errors[gic_best_pos], 
                 mu_error[bic_best_pos],
                 mu_error[aic_best_pos],
                 mu_error[gic_best_pos],
                 gamma_l1_error_matrix[bic_best_pos],
                 gamma_l2_error_matrix[bic_best_pos],
                 gamma_l1_error_matrix[aic_best_pos],
                 gamma_l2_error_matrix[aic_best_pos],
                 gamma_l1_error_matrix[gic_best_pos],
                 gamma_l2_error_matrix[gic_best_pos])
  
  # Misclassification errors
  misclassification_matrix <- c(apply(fit$beta[1:n,], MARGIN = 2, actives_as_inactive), 
                                apply(fit$beta[1:n,], MARGIN = 2, inactives_as_active))
  
  ### Inference
  
  # BIC
  hat_gamma <- fit$beta[(n+1):(n+p),bic_best_pos]
  Z <- cbind(Z_1, Z_2)
  W_hat <- sparseMatrix(i = 1:choose(n,2), j=1:choose(n,2), x = fitted_probabilities[,bic_best_pos]*(1-fitted_probabilities[,bic_best_pos]))
  sample_gram_matrix <- 1/choose(n,2)*t(Z)%*%W_hat%*%Z
  root_Theta_hat_diag <- sample_gram_matrix%>%as.matrix()%>%solve()%>%diag()%>%sqrt()
  
  standardized_gamma_matrix <- sqrt(choose(n,2))*(hat_gamma - gamma_star)/root_Theta_hat_diag
  CI_gamma <- rbind(hat_gamma + qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)),
                    hat_gamma - qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)))
  CI_lengths <- CI_gamma[1,] - CI_gamma[2,]
  gamma_1_coverage <- (gamma_star[1] <= CI_gamma[1,1]) && (gamma_star[1] >= CI_gamma[2,1])
  gamma_2_coverage <- (gamma_star[2] <= CI_gamma[1,2]) && (gamma_star[2] >= CI_gamma[2,2])
  confidence_intervals_gamma_matrix <- c(gamma_1_coverage, CI_lengths[1], gamma_2_coverage, CI_lengths[2])
  
  
  readr::write_csv(as.data.frame(t(c(standardized_gamma_matrix, confidence_intervals_gamma_matrix))), 
                   path = paste0(SGE_ID, "_", n, "_BIC_gamma_inference.csv"), append = TRUE)
  
  
  # Fixed lambda
  hat_theta <- coef(fit, s = lambda, exact = TRUE, x=X,y=Y, family = "binomial", alpha = 1, penalty.factor = c(rep(1,n),0,0),
                    lower.limits = c(rep(0,n), -Inf, -Inf))
  
  beta_l1_error_fl <- abs(hat_theta[2:(n+1)] - betas)%>%sum()
  mu_error_fl <- (hat_theta[1] - mu)
  gamma_l1_error_fl <- abs(hat_theta[(n+2):(n+1+p)] - gamma_star)%>%sum()
  
  actives_as_inactive_fl <- actives_as_inactive(hat_theta[2:(n+1)])
  inactives_as_active_fl <- inactives_as_active(hat_theta[2:(n+1)])  
  
  
  readr::write_csv(as.data.frame(t(beta_l1_error_fl)), path = paste0(SGE_ID, "_", n, "_FL_beta_l1_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(mu_error_fl)), path = paste0(SGE_ID, "_", n, "_FL_mu_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(gamma_l1_error_fl)), path = paste0(SGE_ID, "_", n, "_FL_gamma_l1_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(c(actives_as_inactive_fl, inactives_as_active_fl))), path = paste0(SGE_ID, "_", n, "_FL_misclassification.csv"), append = TRUE)
  
  rm(fitted_probabilities)
  fixed_lambda_probabilities <- predict(fit, newx = X, s=lambda, type = "response")%>%as.vector()
  
  hat_gamma <- hat_theta[(n+2):(n+1+p)]
  #Z <- cbind(Z_1, Z_2)
  W_hat <- sparseMatrix(i = 1:choose(n,2), j=1:choose(n,2), 
                        x = fixed_lambda_probabilities*(1-fixed_lambda_probabilities))
  sample_gram_matrix <- 1/choose(n,2)*t(Z)%*%W_hat%*%Z
  root_Theta_hat_diag <- sample_gram_matrix%>%as.matrix()%>%solve()%>%diag()%>%sqrt()
  
  standardized_gamma_matrix <- sqrt(choose(n,2))*(hat_gamma - gamma_star)/root_Theta_hat_diag
  CI_gamma <- rbind(hat_gamma + qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)),
                    hat_gamma - qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)))
  CI_lengths <- CI_gamma[1,] - CI_gamma[2,]
  gamma_1_coverage <- (gamma_star[1] <= CI_gamma[1,1]) && (gamma_star[1] >= CI_gamma[2,1])
  gamma_2_coverage <- (gamma_star[2] <= CI_gamma[1,2]) && (gamma_star[2] >= CI_gamma[2,2])
  confidence_intervals_gamma_matrix <- c(gamma_1_coverage, CI_lengths[1], gamma_2_coverage, CI_lengths[2])
  
  readr::write_csv(as.data.frame(t(c(standardized_gamma_matrix, confidence_intervals_gamma_matrix))), 
                   path = paste0(SGE_ID, "_", n, "_FL_gamma_inference.csv"), append = TRUE)
  
  
  
  # Other summaries
  
  readr::write_csv(as.data.frame(t(network_summaries)), 
                   path = paste0(SGE_ID, "_", n, "_network_summaries.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(IC_matrix)), 
                   path = paste0(SGE_ID, "_", n, "_IC_matrix.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(log_lambda_matrix)),
                   path = paste0(SGE_ID, "_", n, "_log_lambda_matrix.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(gamma_l1_error_matrix)), 
                   path = paste0(SGE_ID, "_", n, "_gamma_l1_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(gamma_l2_error_matrix)), 
                   path = paste0(SGE_ID, "_", n, "_gamma_l2_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(misclassification_matrix)), 
                   path = paste0(SGE_ID, "_", n, "_misclassification_matrix.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(fit$a0 - mu)), 
                   path = paste0(SGE_ID, "_", n, "_mu_error.csv"), append = TRUE)
  
  rm(fixed_lambda_probabilities, log_like_vector, adj_matrix, fit)
  
}
end_time <- Sys.time()
print(end_time - start_time)


# Clear workspace ---------------------------------------------------------

rm(list = ls())


# n=500 -------------------------------------------------------------------

# Number nodes
n <- 500
p <- 2
#s <- 7 #300
s <- 9 # 500
#s <- 12 # 1000

M <- 100

betas <- c(1.2,0.8,rep(1,s-2),rep(0,n-s))
mu <- -0.5*log(log(n))
gamma_star <- c(1,0.8)

# For fixed lambda
# Penaltiy
a <- sqrt(2*log(2*(n+p+1))/choose(n,2))
t <- 2
lambda_0 <- 4*a + 2*sqrt(2*t/choose(n,2) * (1+ sqrt(2*n) * a)) + t*sqrt(2*n)/(3*choose(n,2))
lambda_bar <- 8*lambda_0
lambda <- sqrt(2/n)*lambda_0

rm(a,t,lambda_0, lambda_bar)

# GIC parameter
a_n <- log(log(choose(n,2)))*log(n+1)

# Helper functions --------------------------------------------------------

p_ij <- function(mu, beta_i, beta_j, i, j){
  Z_ij <- c(Z_1_matrix[j,i], Z_2_matrix[j,i])%*%gamma_star
  1/(1 + exp(-(beta_i + beta_j + mu + Z_ij)))
}

log_like <- function(adj, probs){
  return(sum(adj%*%log(probs) + (1-adj)%*%log(1 - probs)))
}

actives_as_inactive <- function(vec){
  return(sum(which(vec == 0) <= s))
}
inactives_as_active <- function(vec){
  return(sum(which(vec > 0) > s))
}

# Create sparse design matrix ---------------------------------------------

# Sparse design matrix (same for each iteration, so created outside of loop)

# row indices: we take all rows (no empty rows) and two entries per row
rows <- rep(1:choose(n,2), each = 2)
# column indices
a <- 2:n
b <- rep(n, n-1)
c <- NULL
for (i in 1:(n-1)) {
  c <- c(c, seq(a[i],b[i]))
}
cols <- c(rbind(rep(1:(n-1), (n-1):1),c))

X_deterministic <- sparseMatrix(i = rows, j = cols, x = 1)
rm(a,b,c, rows, cols, i)


# Environment -------------------------------------------------------------

SGE_ID <-  as.numeric(Sys.getenv('SLURM_JOB_ID'))
set.seed(SGE_ID)

# Monte Carlo Start -------------------------------------------------------

start_time <- Sys.time()
for (m in 1:M) {
  print(m)
  
  # Similarity matrix for first covariate
  Z_1_matrix <- matrix(nrow = n, ncol = n, data = rbeta(n*n, 2, 2) - 1/2)
  # Similarity values for first covariate in vector form
  Z_1 <- Z_1_matrix[lower.tri(Z_1_matrix)]
  # Similarity matrix for second covariate
  Z_2_matrix <- matrix(nrow = n, ncol = n, data = rbeta(n*n, 2, 2) - 1/2)
  # Similarity values for second covariate in vector form
  Z_2 <- Z_2_matrix[lower.tri(Z_2_matrix)]
  
  probabilities <- matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in i:n) {
      probabilities[j,i] <- p_ij(mu, betas[i], betas[j], i, j)
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
  rm(probabilities, Z_1_matrix, Z_2_matrix)
  # Record Network Summaries ------------------------------------------------
  network_summaries <- numeric()
  degrees_of_freedom <- numeric()
  log_lambda_matrix <- numeric()
  gamma_l1_error_matrix <- numeric()
  gamma_l2_error_matrix <- numeric()
  relative_l2_error_matrix <- numeric()
  BICs <- numeric()
  AICs <- numeric()
  GICs <- numeric()
  IC_matrix <- numeric()
  misclassification_matrix <- numeric()
  
  
  degrees <- rowSums(adj_matrix, na.rm = TRUE)
  network_summaries[1:3] <- c(s,
                              length(unique(degrees)),
                              sum(degrees)/2)
  
  # Create fit --------------------------------------------------------------
  
  
  # Design matrix
  X <- cbind(X_deterministic, Z_1, Z_2)
  
  # Values Y; I have to take lower.tri to get the right order of entries
  Y <- adj_matrix[lower.tri(adj_matrix)]
  # we can leave the default standardize = TRUE option, results are essentially the same
  fit <- glmnet::glmnet(X,Y, family = "binomial", alpha = 1, penalty.factor = c(rep(1,n),0,0),
                        lower.limits = c(rep(0,n), -Inf, -Inf))
  
  
  
  # Update Network summaries ------------------------------------------------
  
  # We only want to know number of active BETAS; also remember:
  # unpenalized variables, i.e. gamma, will always be active
  degrees_of_freedom <- (fit$df - p)
  
  log_lambda_matrix <- log(fit$lambda)
  covariate_error <- abs(as.matrix(fit$beta[(n+1):(n+p),]) - gamma_star)
  # Gamma errors
  gamma_l1_error_matrix <- colSums(covariate_error)
  gamma_l2_error_matrix <- sqrt(colSums(covariate_error^2))
  
  # Error Statistics
  # beta
  beta_abs_errors <- (abs(as.matrix(fit$beta[1:n,]) - betas))
  # l_1 errors
  beta_mae_errors <- colSums(beta_abs_errors)/n
  # l_2 errors: square all entries, build column sums, divide by n+1
  beta_l2_errors <- sqrt(colSums(beta_abs_errors^2))
  beta_mse_errors <- beta_l2_errors/n
  
  # mu
  mu_error <- fit$a0 - mu
  
  
  # log likelihood
  fitted_probabilities <- predict(fit, newx = X, type = "response")
  # This is the POSITIVE log likelihood for each fitted model.
  log_like_vector <- apply(fitted_probabilities, 
                           MARGIN = 2, 
                           function(i) log_like(adj = adj_matrix[lower.tri(adj_matrix)], probs = i))
  
  bic <- fit$df*log(choose(n,2)) - 2*log_like_vector
  aic <- 2*(fit$df - log_like_vector)
  gic <- fit$df*a_n - 2*log_like_vector
  BICs <- bic
  AICs <- aic
  GICs <- gic
  bic_best_pos <- which.min(bic)
  aic_best_pos <- which.min(aic)
  gic_best_pos <- which.min(gic)
  
  
  IC_matrix <- c(bic[bic_best_pos], # optimal bic value
                 bic_best_pos, # position of optimal bic value
                 fit$df[bic_best_pos] - p, # how many degrees of freedom (=active indices)?
                 sum(which(coef(fit, s = fit$lambda[bic_best_pos])[2:(n+1)] == 0) <= s), # how many actives were misclassified as inactive?
                 sum(which(coef(fit, s = fit$lambda[bic_best_pos])[2:(n+1)] > 0) > s), # how many inactives as active?
                 sum(coef(fit, s = fit$lambda[bic_best_pos])[(s+2):(n+1)])/(n-s), # mean l1 norm of misclas inactive betas
                 beta_mae_errors[bic_best_pos], # mae
                 beta_mse_errors[bic_best_pos], #mse
                 aic[aic_best_pos],
                 aic_best_pos,
                 fit$df[aic_best_pos] - p,
                 sum(which(coef(fit, s = fit$lambda[aic_best_pos])[2:(n+1)] == 0) <= s),
                 sum(which(coef(fit, s = fit$lambda[aic_best_pos])[2:(n+1)] > 0) > s),
                 sum(coef(fit, s = fit$lambda[aic_best_pos])[(s+2):(n+1)])/(n-s),
                 beta_mae_errors[aic_best_pos], 
                 beta_mse_errors[aic_best_pos], 
                 gic[gic_best_pos],
                 gic_best_pos,
                 fit$df[gic_best_pos] - p,
                 sum(which(coef(fit, s = fit$lambda[gic_best_pos])[2:(n+1)] == 0) <= s),
                 sum(which(coef(fit, s = fit$lambda[gic_best_pos])[2:(n+1)] > 0) > s),
                 sum(coef(fit, s = fit$lambda[gic_best_pos])[(s+2):(n+1)])/(n-s),
                 beta_mae_errors[gic_best_pos], 
                 beta_mse_errors[gic_best_pos], 
                 mu_error[bic_best_pos],
                 mu_error[aic_best_pos],
                 mu_error[gic_best_pos],
                 gamma_l1_error_matrix[bic_best_pos],
                 gamma_l2_error_matrix[bic_best_pos],
                 gamma_l1_error_matrix[aic_best_pos],
                 gamma_l2_error_matrix[aic_best_pos],
                 gamma_l1_error_matrix[gic_best_pos],
                 gamma_l2_error_matrix[gic_best_pos])
  
  # Misclassification errors
  misclassification_matrix <- c(apply(fit$beta[1:n,], MARGIN = 2, actives_as_inactive), 
                                apply(fit$beta[1:n,], MARGIN = 2, inactives_as_active))
  
  ### Inference
  
  # BIC
  hat_gamma <- fit$beta[(n+1):(n+p),bic_best_pos]
  Z <- cbind(Z_1, Z_2)
  W_hat <- sparseMatrix(i = 1:choose(n,2), j=1:choose(n,2), x = fitted_probabilities[,bic_best_pos]*(1-fitted_probabilities[,bic_best_pos]))
  sample_gram_matrix <- 1/choose(n,2)*t(Z)%*%W_hat%*%Z
  root_Theta_hat_diag <- sample_gram_matrix%>%as.matrix()%>%solve()%>%diag()%>%sqrt()
  
  standardized_gamma_matrix <- sqrt(choose(n,2))*(hat_gamma - gamma_star)/root_Theta_hat_diag
  CI_gamma <- rbind(hat_gamma + qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)),
                    hat_gamma - qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)))
  CI_lengths <- CI_gamma[1,] - CI_gamma[2,]
  gamma_1_coverage <- (gamma_star[1] <= CI_gamma[1,1]) && (gamma_star[1] >= CI_gamma[2,1])
  gamma_2_coverage <- (gamma_star[2] <= CI_gamma[1,2]) && (gamma_star[2] >= CI_gamma[2,2])
  confidence_intervals_gamma_matrix <- c(gamma_1_coverage, CI_lengths[1], gamma_2_coverage, CI_lengths[2])
  
  
  readr::write_csv(as.data.frame(t(c(standardized_gamma_matrix, confidence_intervals_gamma_matrix))), 
                   path = paste0(SGE_ID, "_", n, "_BIC_gamma_inference.csv"), append = TRUE)
  
  
  # Fixed lambda
  hat_theta <- coef(fit, s = lambda, exact = TRUE, x=X,y=Y, family = "binomial", alpha = 1, penalty.factor = c(rep(1,n),0,0),
                    lower.limits = c(rep(0,n), -Inf, -Inf))
  
  beta_l1_error_fl <- abs(hat_theta[2:(n+1)] - betas)%>%sum()
  mu_error_fl <- (hat_theta[1] - mu)
  gamma_l1_error_fl <- abs(hat_theta[(n+2):(n+1+p)] - gamma_star)%>%sum()
  
  actives_as_inactive_fl <- actives_as_inactive(hat_theta[2:(n+1)])
  inactives_as_active_fl <- inactives_as_active(hat_theta[2:(n+1)])  
  
  
  readr::write_csv(as.data.frame(t(beta_l1_error_fl)), path = paste0(SGE_ID, "_", n, "_FL_beta_l1_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(mu_error_fl)), path = paste0(SGE_ID, "_", n, "_FL_mu_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(gamma_l1_error_fl)), path = paste0(SGE_ID, "_", n, "_FL_gamma_l1_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(c(actives_as_inactive_fl, inactives_as_active_fl))), path = paste0(SGE_ID, "_", n, "_FL_misclassification.csv"), append = TRUE)
  
  rm(fitted_probabilities)
  fixed_lambda_probabilities <- predict(fit, newx = X, s=lambda, type = "response")%>%as.vector()
  
  hat_gamma <- hat_theta[(n+2):(n+1+p)]
  #Z <- cbind(Z_1, Z_2)
  W_hat <- sparseMatrix(i = 1:choose(n,2), j=1:choose(n,2), 
                        x = fixed_lambda_probabilities*(1-fixed_lambda_probabilities))
  sample_gram_matrix <- 1/choose(n,2)*t(Z)%*%W_hat%*%Z
  root_Theta_hat_diag <- sample_gram_matrix%>%as.matrix()%>%solve()%>%diag()%>%sqrt()
  
  standardized_gamma_matrix <- sqrt(choose(n,2))*(hat_gamma - gamma_star)/root_Theta_hat_diag
  CI_gamma <- rbind(hat_gamma + qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)),
                    hat_gamma - qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)))
  CI_lengths <- CI_gamma[1,] - CI_gamma[2,]
  gamma_1_coverage <- (gamma_star[1] <= CI_gamma[1,1]) && (gamma_star[1] >= CI_gamma[2,1])
  gamma_2_coverage <- (gamma_star[2] <= CI_gamma[1,2]) && (gamma_star[2] >= CI_gamma[2,2])
  confidence_intervals_gamma_matrix <- c(gamma_1_coverage, CI_lengths[1], gamma_2_coverage, CI_lengths[2])
  
  readr::write_csv(as.data.frame(t(c(standardized_gamma_matrix, confidence_intervals_gamma_matrix))), 
                   path = paste0(SGE_ID, "_", n, "_FL_gamma_inference.csv"), append = TRUE)
  
  
  
  # Other summaries
  
  readr::write_csv(as.data.frame(t(network_summaries)), 
                   path = paste0(SGE_ID, "_", n, "_network_summaries.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(IC_matrix)), 
                   path = paste0(SGE_ID, "_", n, "_IC_matrix.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(log_lambda_matrix)),
                   path = paste0(SGE_ID, "_", n, "_log_lambda_matrix.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(gamma_l1_error_matrix)), 
                   path = paste0(SGE_ID, "_", n, "_gamma_l1_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(gamma_l2_error_matrix)), 
                   path = paste0(SGE_ID, "_", n, "_gamma_l2_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(misclassification_matrix)), 
                   path = paste0(SGE_ID, "_", n, "_misclassification_matrix.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(fit$a0 - mu)), 
                   path = paste0(SGE_ID, "_", n, "_mu_error.csv"), append = TRUE)
  
  rm(fixed_lambda_probabilities, log_like_vector, adj_matrix, fit)
  
}
end_time <- Sys.time()
print(end_time - start_time)



# Clear workspace ---------------------------------------------------------

rm(list = ls())

# 800 ---------------------------------------------------------------------



# Number nodes
n <- 800
p <- 2
#s <- 7 #300
#s <- 9 # 500
s <- 10 # 800
#s <- 12 # 1000

M <- 100

betas <- c(1.2,0.8,rep(1,s-2),rep(0,n-s))
mu <- -0.5*log(log(n))
gamma_star <- c(1,0.8)

# For fixed lambda
# Penaltiy
a <- sqrt(2*log(2*(n+p+1))/choose(n,2))
t <- 2
lambda_0 <- 4*a + 2*sqrt(2*t/choose(n,2) * (1+ sqrt(2*n) * a)) + t*sqrt(2*n)/(3*choose(n,2))
lambda_bar <- 8*lambda_0
lambda <- sqrt(2/n)*lambda_0

rm(a,t,lambda_0, lambda_bar)

# GIC parameter
a_n <- log(log(choose(n,2)))*log(n+1)

# Helper functions --------------------------------------------------------

p_ij <- function(mu, beta_i, beta_j, i, j){
  Z_ij <- c(Z_1_matrix[j,i], Z_2_matrix[j,i])%*%gamma_star
  1/(1 + exp(-(beta_i + beta_j + mu + Z_ij)))
}

log_like <- function(adj, probs){
  return(sum(adj%*%log(probs) + (1-adj)%*%log(1 - probs)))
}

actives_as_inactive <- function(vec){
  return(sum(which(vec == 0) <= s))
}
inactives_as_active <- function(vec){
  return(sum(which(vec > 0) > s))
}

# Create sparse design matrix ---------------------------------------------

# Sparse design matrix (same for each iteration, so created outside of loop)

# row indices: we take all rows (no empty rows) and two entries per row
rows <- rep(1:choose(n,2), each = 2)
# column indices
a <- 2:n
b <- rep(n, n-1)
c <- NULL
for (i in 1:(n-1)) {
  c <- c(c, seq(a[i],b[i]))
}
cols <- c(rbind(rep(1:(n-1), (n-1):1),c))

X_deterministic <- sparseMatrix(i = rows, j = cols, x = 1)
rm(a,b,c, rows, cols, i)


# Environment -------------------------------------------------------------

SGE_ID <-  as.numeric(Sys.getenv('SLURM_JOB_ID'))
set.seed(SGE_ID)

# Monte Carlo Start -------------------------------------------------------

start_time <- Sys.time()
for (m in 1:M) {
  print(m)
  
  # Similarity matrix for first covariate
  Z_1_matrix <- matrix(nrow = n, ncol = n, data = rbeta(n*n, 2, 2) - 1/2)
  # Similarity values for first covariate in vector form
  Z_1 <- Z_1_matrix[lower.tri(Z_1_matrix)]
  # Similarity matrix for second covariate
  Z_2_matrix <- matrix(nrow = n, ncol = n, data = rbeta(n*n, 2, 2) - 1/2)
  # Similarity values for second covariate in vector form
  Z_2 <- Z_2_matrix[lower.tri(Z_2_matrix)]
  
  probabilities <- matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in i:n) {
      probabilities[j,i] <- p_ij(mu, betas[i], betas[j], i, j)
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
  rm(probabilities, Z_1_matrix, Z_2_matrix)
  # Record Network Summaries ------------------------------------------------
  network_summaries <- numeric()
  degrees_of_freedom <- numeric()
  log_lambda_matrix <- numeric()
  gamma_l1_error_matrix <- numeric()
  gamma_l2_error_matrix <- numeric()
  relative_l2_error_matrix <- numeric()
  BICs <- numeric()
  AICs <- numeric()
  GICs <- numeric()
  IC_matrix <- numeric()
  misclassification_matrix <- numeric()
  
  
  degrees <- rowSums(adj_matrix, na.rm = TRUE)
  network_summaries[1:3] <- c(s,
                              length(unique(degrees)),
                              sum(degrees)/2)
  
  # Create fit --------------------------------------------------------------
  
  
  # Design matrix
  X <- cbind(X_deterministic, Z_1, Z_2)
  
  # Values Y; I have to take lower.tri to get the right order of entries
  Y <- adj_matrix[lower.tri(adj_matrix)]
  # we can leave the default standardize = TRUE option, results are essentially the same
  fit <- glmnet::glmnet(X,Y, family = "binomial", alpha = 1, penalty.factor = c(rep(1,n),0,0),
                        lower.limits = c(rep(0,n), -Inf, -Inf))
  
  
  
  # Update Network summaries ------------------------------------------------
  
  # We only want to know number of active BETAS; also remember:
  # unpenalized variables, i.e. gamma, will always be active
  degrees_of_freedom <- (fit$df - p)
  
  log_lambda_matrix <- log(fit$lambda)
  covariate_error <- abs(as.matrix(fit$beta[(n+1):(n+p),]) - gamma_star)
  # Gamma errors
  gamma_l1_error_matrix <- colSums(covariate_error)
  gamma_l2_error_matrix <- sqrt(colSums(covariate_error^2))
  
  # Error Statistics
  # beta
  beta_abs_errors <- (abs(as.matrix(fit$beta[1:n,]) - betas))
  # l_1 errors
  beta_mae_errors <- colSums(beta_abs_errors)/n
  # l_2 errors: square all entries, build column sums, divide by n+1
  beta_l2_errors <- sqrt(colSums(beta_abs_errors^2))
  beta_mse_errors <- beta_l2_errors/n
  
  # mu
  mu_error <- fit$a0 - mu
  
  
  # log likelihood
  fitted_probabilities <- predict(fit, newx = X, type = "response")
  # This is the POSITIVE log likelihood for each fitted model.
  log_like_vector <- apply(fitted_probabilities, 
                           MARGIN = 2, 
                           function(i) log_like(adj = adj_matrix[lower.tri(adj_matrix)], probs = i))
  
  bic <- fit$df*log(choose(n,2)) - 2*log_like_vector
  aic <- 2*(fit$df - log_like_vector)
  gic <- fit$df*a_n - 2*log_like_vector
  BICs <- bic
  AICs <- aic
  GICs <- gic
  bic_best_pos <- which.min(bic)
  aic_best_pos <- which.min(aic)
  gic_best_pos <- which.min(gic)
  
  
  IC_matrix <- c(bic[bic_best_pos], # optimal bic value
                 bic_best_pos, # position of optimal bic value
                 fit$df[bic_best_pos] - p, # how many degrees of freedom (=active indices)?
                 sum(which(coef(fit, s = fit$lambda[bic_best_pos])[2:(n+1)] == 0) <= s), # how many actives were misclassified as inactive?
                 sum(which(coef(fit, s = fit$lambda[bic_best_pos])[2:(n+1)] > 0) > s), # how many inactives as active?
                 sum(coef(fit, s = fit$lambda[bic_best_pos])[(s+2):(n+1)])/(n-s), # mean l1 norm of misclas inactive betas
                 beta_mae_errors[bic_best_pos], # mae
                 beta_mse_errors[bic_best_pos], #mse
                 aic[aic_best_pos],
                 aic_best_pos,
                 fit$df[aic_best_pos] - p,
                 sum(which(coef(fit, s = fit$lambda[aic_best_pos])[2:(n+1)] == 0) <= s),
                 sum(which(coef(fit, s = fit$lambda[aic_best_pos])[2:(n+1)] > 0) > s),
                 sum(coef(fit, s = fit$lambda[aic_best_pos])[(s+2):(n+1)])/(n-s),
                 beta_mae_errors[aic_best_pos], 
                 beta_mse_errors[aic_best_pos], 
                 gic[gic_best_pos],
                 gic_best_pos,
                 fit$df[gic_best_pos] - p,
                 sum(which(coef(fit, s = fit$lambda[gic_best_pos])[2:(n+1)] == 0) <= s),
                 sum(which(coef(fit, s = fit$lambda[gic_best_pos])[2:(n+1)] > 0) > s),
                 sum(coef(fit, s = fit$lambda[gic_best_pos])[(s+2):(n+1)])/(n-s),
                 beta_mae_errors[gic_best_pos], 
                 beta_mse_errors[gic_best_pos], 
                 mu_error[bic_best_pos],
                 mu_error[aic_best_pos],
                 mu_error[gic_best_pos],
                 gamma_l1_error_matrix[bic_best_pos],
                 gamma_l2_error_matrix[bic_best_pos],
                 gamma_l1_error_matrix[aic_best_pos],
                 gamma_l2_error_matrix[aic_best_pos],
                 gamma_l1_error_matrix[gic_best_pos],
                 gamma_l2_error_matrix[gic_best_pos])
  
  # Misclassification errors
  misclassification_matrix <- c(apply(fit$beta[1:n,], MARGIN = 2, actives_as_inactive), 
                                apply(fit$beta[1:n,], MARGIN = 2, inactives_as_active))
  
  ### Inference
  
  # BIC
  hat_gamma <- fit$beta[(n+1):(n+p),bic_best_pos]
  Z <- cbind(Z_1, Z_2)
  W_hat <- sparseMatrix(i = 1:choose(n,2), j=1:choose(n,2), x = fitted_probabilities[,bic_best_pos]*(1-fitted_probabilities[,bic_best_pos]))
  sample_gram_matrix <- 1/choose(n,2)*t(Z)%*%W_hat%*%Z
  root_Theta_hat_diag <- sample_gram_matrix%>%as.matrix()%>%solve()%>%diag()%>%sqrt()
  
  standardized_gamma_matrix <- sqrt(choose(n,2))*(hat_gamma - gamma_star)/root_Theta_hat_diag
  CI_gamma <- rbind(hat_gamma + qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)),
                    hat_gamma - qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)))
  CI_lengths <- CI_gamma[1,] - CI_gamma[2,]
  gamma_1_coverage <- (gamma_star[1] <= CI_gamma[1,1]) && (gamma_star[1] >= CI_gamma[2,1])
  gamma_2_coverage <- (gamma_star[2] <= CI_gamma[1,2]) && (gamma_star[2] >= CI_gamma[2,2])
  confidence_intervals_gamma_matrix <- c(gamma_1_coverage, CI_lengths[1], gamma_2_coverage, CI_lengths[2])
  
  
  readr::write_csv(as.data.frame(t(c(standardized_gamma_matrix, confidence_intervals_gamma_matrix))), 
                   path = paste0(SGE_ID, "_", n, "_BIC_gamma_inference.csv"), append = TRUE)
  
  
  # Fixed lambda
  hat_theta <- coef(fit, s = lambda, exact = TRUE, x=X,y=Y, family = "binomial", alpha = 1, penalty.factor = c(rep(1,n),0,0),
                    lower.limits = c(rep(0,n), -Inf, -Inf))
  
  beta_l1_error_fl <- abs(hat_theta[2:(n+1)] - betas)%>%sum()
  mu_error_fl <- (hat_theta[1] - mu)
  gamma_l1_error_fl <- abs(hat_theta[(n+2):(n+1+p)] - gamma_star)%>%sum()
  
  actives_as_inactive_fl <- actives_as_inactive(hat_theta[2:(n+1)])
  inactives_as_active_fl <- inactives_as_active(hat_theta[2:(n+1)])  
  
  
  readr::write_csv(as.data.frame(t(beta_l1_error_fl)), path = paste0(SGE_ID, "_", n, "_FL_beta_l1_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(mu_error_fl)), path = paste0(SGE_ID, "_", n, "_FL_mu_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(gamma_l1_error_fl)), path = paste0(SGE_ID, "_", n, "_FL_gamma_l1_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(c(actives_as_inactive_fl, inactives_as_active_fl))), path = paste0(SGE_ID, "_", n, "_FL_misclassification.csv"), append = TRUE)
  
  rm(fitted_probabilities)
  fixed_lambda_probabilities <- predict(fit, newx = X, s=lambda, type = "response")%>%as.vector()
  
  hat_gamma <- hat_theta[(n+2):(n+1+p)]
  #Z <- cbind(Z_1, Z_2)
  W_hat <- sparseMatrix(i = 1:choose(n,2), j=1:choose(n,2), 
                        x = fixed_lambda_probabilities*(1-fixed_lambda_probabilities))
  sample_gram_matrix <- 1/choose(n,2)*t(Z)%*%W_hat%*%Z
  root_Theta_hat_diag <- sample_gram_matrix%>%as.matrix()%>%solve()%>%diag()%>%sqrt()
  
  standardized_gamma_matrix <- sqrt(choose(n,2))*(hat_gamma - gamma_star)/root_Theta_hat_diag
  CI_gamma <- rbind(hat_gamma + qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)),
                    hat_gamma - qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)))
  CI_lengths <- CI_gamma[1,] - CI_gamma[2,]
  gamma_1_coverage <- (gamma_star[1] <= CI_gamma[1,1]) && (gamma_star[1] >= CI_gamma[2,1])
  gamma_2_coverage <- (gamma_star[2] <= CI_gamma[1,2]) && (gamma_star[2] >= CI_gamma[2,2])
  confidence_intervals_gamma_matrix <- c(gamma_1_coverage, CI_lengths[1], gamma_2_coverage, CI_lengths[2])
  
  readr::write_csv(as.data.frame(t(c(standardized_gamma_matrix, confidence_intervals_gamma_matrix))), 
                   path = paste0(SGE_ID, "_", n, "_FL_gamma_inference.csv"), append = TRUE)
  
  
  
  # Other summaries
  
  readr::write_csv(as.data.frame(t(network_summaries)), 
                   path = paste0(SGE_ID, "_", n, "_network_summaries.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(IC_matrix)), 
                   path = paste0(SGE_ID, "_", n, "_IC_matrix.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(log_lambda_matrix)),
                   path = paste0(SGE_ID, "_", n, "_log_lambda_matrix.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(gamma_l1_error_matrix)), 
                   path = paste0(SGE_ID, "_", n, "_gamma_l1_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(gamma_l2_error_matrix)), 
                   path = paste0(SGE_ID, "_", n, "_gamma_l2_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(misclassification_matrix)), 
                   path = paste0(SGE_ID, "_", n, "_misclassification_matrix.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(fit$a0 - mu)), 
                   path = paste0(SGE_ID, "_", n, "_mu_error.csv"), append = TRUE)
  
  rm(fixed_lambda_probabilities, log_like_vector, adj_matrix, fit)
  
}
end_time <- Sys.time()
print(end_time - start_time)



# Clear workspace ---------------------------------------------------------

rm(list = ls())


# n=1000 ------------------------------------------------------------------

# Number nodes
n <- 1000
p <- 2
#s <- 7 #300
#s <- 9 # 500
s <- 12 # 1000

M <- 100

betas <- c(1.2,0.8,rep(1,s-2),rep(0,n-s))
mu <- -0.5*log(log(n))
gamma_star <- c(1,0.8)

# For fixed lambda
# Penaltiy
a <- sqrt(2*log(2*(n+p+1))/choose(n,2))
t <- 2
lambda_0 <- 4*a + 2*sqrt(2*t/choose(n,2) * (1+ sqrt(2*n) * a)) + t*sqrt(2*n)/(3*choose(n,2))
lambda_bar <- 8*lambda_0
lambda <- sqrt(2/n)*lambda_0

rm(a,t,lambda_0, lambda_bar)

# GIC parameter
a_n <- log(log(choose(n,2)))*log(n+1)

# Helper functions --------------------------------------------------------

p_ij <- function(mu, beta_i, beta_j, i, j){
  Z_ij <- c(Z_1_matrix[j,i], Z_2_matrix[j,i])%*%gamma_star
  1/(1 + exp(-(beta_i + beta_j + mu + Z_ij)))
}

log_like <- function(adj, probs){
  return(sum(adj%*%log(probs) + (1-adj)%*%log(1 - probs)))
}

actives_as_inactive <- function(vec){
  return(sum(which(vec == 0) <= s))
}
inactives_as_active <- function(vec){
  return(sum(which(vec > 0) > s))
}

# Create sparse design matrix ---------------------------------------------

# Sparse design matrix (same for each iteration, so created outside of loop)

# row indices: we take all rows (no empty rows) and two entries per row
rows <- rep(1:choose(n,2), each = 2)
# column indices
a <- 2:n
b <- rep(n, n-1)
c <- NULL
for (i in 1:(n-1)) {
  c <- c(c, seq(a[i],b[i]))
}
cols <- c(rbind(rep(1:(n-1), (n-1):1),c))

X_deterministic <- sparseMatrix(i = rows, j = cols, x = 1)
rm(a,b,c, rows, cols, i)


# Environment -------------------------------------------------------------

SGE_ID <-  as.numeric(Sys.getenv('SLURM_JOB_ID'))
set.seed(SGE_ID)

# Monte Carlo Start -------------------------------------------------------

start_time <- Sys.time()
for (m in 1:M) {
  print(m)
  
  # Similarity matrix for first covariate
  Z_1_matrix <- matrix(nrow = n, ncol = n, data = rbeta(n*n, 2, 2) - 1/2)
  # Similarity values for first covariate in vector form
  Z_1 <- Z_1_matrix[lower.tri(Z_1_matrix)]
  # Similarity matrix for second covariate
  Z_2_matrix <- matrix(nrow = n, ncol = n, data = rbeta(n*n, 2, 2) - 1/2)
  # Similarity values for second covariate in vector form
  Z_2 <- Z_2_matrix[lower.tri(Z_2_matrix)]
  
  probabilities <- matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in i:n) {
      probabilities[j,i] <- p_ij(mu, betas[i], betas[j], i, j)
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
  rm(probabilities, Z_1_matrix, Z_2_matrix)
  # Record Network Summaries ------------------------------------------------
  network_summaries <- numeric()
  degrees_of_freedom <- numeric()
  log_lambda_matrix <- numeric()
  gamma_l1_error_matrix <- numeric()
  gamma_l2_error_matrix <- numeric()
  relative_l2_error_matrix <- numeric()
  BICs <- numeric()
  AICs <- numeric()
  GICs <- numeric()
  IC_matrix <- numeric()
  misclassification_matrix <- numeric()
  
  
  degrees <- rowSums(adj_matrix, na.rm = TRUE)
  network_summaries[1:3] <- c(s,
                              length(unique(degrees)),
                              sum(degrees)/2)
  
  # Create fit --------------------------------------------------------------
  
  
  # Design matrix
  X <- cbind(X_deterministic, Z_1, Z_2)
  
  # Values Y; I have to take lower.tri to get the right order of entries
  Y <- adj_matrix[lower.tri(adj_matrix)]
  # we can leave the default standardize = TRUE option, results are essentially the same
  fit <- glmnet::glmnet(X,Y, family = "binomial", alpha = 1, penalty.factor = c(rep(1,n),0,0),
                        lower.limits = c(rep(0,n), -Inf, -Inf))
  
  
  
  # Update Network summaries ------------------------------------------------
  
  # We only want to know number of active BETAS; also remember:
  # unpenalized variables, i.e. gamma, will always be active
  degrees_of_freedom <- (fit$df - p)
  
  log_lambda_matrix <- log(fit$lambda)
  covariate_error <- abs(as.matrix(fit$beta[(n+1):(n+p),]) - gamma_star)
  # Gamma errors
  gamma_l1_error_matrix <- colSums(covariate_error)
  gamma_l2_error_matrix <- sqrt(colSums(covariate_error^2))
  
  # Error Statistics
  # beta
  beta_abs_errors <- (abs(as.matrix(fit$beta[1:n,]) - betas))
  # l_1 errors
  beta_mae_errors <- colSums(beta_abs_errors)/n
  # l_2 errors: square all entries, build column sums, divide by n+1
  beta_l2_errors <- sqrt(colSums(beta_abs_errors^2))
  beta_mse_errors <- beta_l2_errors/n
  
  # mu
  mu_error <- fit$a0 - mu
  
  
  # log likelihood
  fitted_probabilities <- predict(fit, newx = X, type = "response")
  # This is the POSITIVE log likelihood for each fitted model.
  log_like_vector <- apply(fitted_probabilities, 
                           MARGIN = 2, 
                           function(i) log_like(adj = adj_matrix[lower.tri(adj_matrix)], probs = i))
  
  bic <- fit$df*log(choose(n,2)) - 2*log_like_vector
  aic <- 2*(fit$df - log_like_vector)
  gic <- fit$df*a_n - 2*log_like_vector
  BICs <- bic
  AICs <- aic
  GICs <- gic
  bic_best_pos <- which.min(bic)
  aic_best_pos <- which.min(aic)
  gic_best_pos <- which.min(gic)
  
  
  IC_matrix <- c(bic[bic_best_pos], # optimal bic value
                 bic_best_pos, # position of optimal bic value
                 fit$df[bic_best_pos] - p, # how many degrees of freedom (=active indices)?
                 sum(which(coef(fit, s = fit$lambda[bic_best_pos])[2:(n+1)] == 0) <= s), # how many actives were misclassified as inactive?
                 sum(which(coef(fit, s = fit$lambda[bic_best_pos])[2:(n+1)] > 0) > s), # how many inactives as active?
                 sum(coef(fit, s = fit$lambda[bic_best_pos])[(s+2):(n+1)])/(n-s), # mean l1 norm of misclas inactive betas
                 beta_mae_errors[bic_best_pos], # mae
                 beta_mse_errors[bic_best_pos], #mse
                 aic[aic_best_pos],
                 aic_best_pos,
                 fit$df[aic_best_pos] - p,
                 sum(which(coef(fit, s = fit$lambda[aic_best_pos])[2:(n+1)] == 0) <= s),
                 sum(which(coef(fit, s = fit$lambda[aic_best_pos])[2:(n+1)] > 0) > s),
                 sum(coef(fit, s = fit$lambda[aic_best_pos])[(s+2):(n+1)])/(n-s),
                 beta_mae_errors[aic_best_pos], 
                 beta_mse_errors[aic_best_pos], 
                 gic[gic_best_pos],
                 gic_best_pos,
                 fit$df[gic_best_pos] - p,
                 sum(which(coef(fit, s = fit$lambda[gic_best_pos])[2:(n+1)] == 0) <= s),
                 sum(which(coef(fit, s = fit$lambda[gic_best_pos])[2:(n+1)] > 0) > s),
                 sum(coef(fit, s = fit$lambda[gic_best_pos])[(s+2):(n+1)])/(n-s),
                 beta_mae_errors[gic_best_pos], 
                 beta_mse_errors[gic_best_pos], 
                 mu_error[bic_best_pos],
                 mu_error[aic_best_pos],
                 mu_error[gic_best_pos],
                 gamma_l1_error_matrix[bic_best_pos],
                 gamma_l2_error_matrix[bic_best_pos],
                 gamma_l1_error_matrix[aic_best_pos],
                 gamma_l2_error_matrix[aic_best_pos],
                 gamma_l1_error_matrix[gic_best_pos],
                 gamma_l2_error_matrix[gic_best_pos])
  
  # Misclassification errors
  misclassification_matrix <- c(apply(fit$beta[1:n,], MARGIN = 2, actives_as_inactive), 
                                apply(fit$beta[1:n,], MARGIN = 2, inactives_as_active))
  
  ### Inference
  
  # BIC
  hat_gamma <- fit$beta[(n+1):(n+p),bic_best_pos]
  Z <- cbind(Z_1, Z_2)
  W_hat <- sparseMatrix(i = 1:choose(n,2), j=1:choose(n,2), x = fitted_probabilities[,bic_best_pos]*(1-fitted_probabilities[,bic_best_pos]))
  sample_gram_matrix <- 1/choose(n,2)*t(Z)%*%W_hat%*%Z
  root_Theta_hat_diag <- sample_gram_matrix%>%as.matrix()%>%solve()%>%diag()%>%sqrt()
  
  standardized_gamma_matrix <- sqrt(choose(n,2))*(hat_gamma - gamma_star)/root_Theta_hat_diag
  CI_gamma <- rbind(hat_gamma + qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)),
                    hat_gamma - qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)))
  CI_lengths <- CI_gamma[1,] - CI_gamma[2,]
  gamma_1_coverage <- (gamma_star[1] <= CI_gamma[1,1]) && (gamma_star[1] >= CI_gamma[2,1])
  gamma_2_coverage <- (gamma_star[2] <= CI_gamma[1,2]) && (gamma_star[2] >= CI_gamma[2,2])
  confidence_intervals_gamma_matrix <- c(gamma_1_coverage, CI_lengths[1], gamma_2_coverage, CI_lengths[2])
  
  
  readr::write_csv(as.data.frame(t(c(standardized_gamma_matrix, confidence_intervals_gamma_matrix))), 
                   path = paste0(SGE_ID, "_", n, "_BIC_gamma_inference.csv"), append = TRUE)
  
  
  # Fixed lambda
  hat_theta <- coef(fit, s = lambda, exact = TRUE, x=X,y=Y, family = "binomial", alpha = 1, penalty.factor = c(rep(1,n),0,0),
                    lower.limits = c(rep(0,n), -Inf, -Inf))
  
  beta_l1_error_fl <- abs(hat_theta[2:(n+1)] - betas)%>%sum()
  mu_error_fl <- (hat_theta[1] - mu)
  gamma_l1_error_fl <- abs(hat_theta[(n+2):(n+1+p)] - gamma_star)%>%sum()
  
  actives_as_inactive_fl <- actives_as_inactive(hat_theta[2:(n+1)])
  inactives_as_active_fl <- inactives_as_active(hat_theta[2:(n+1)])  
  
  
  readr::write_csv(as.data.frame(t(beta_l1_error_fl)), path = paste0(SGE_ID, "_", n, "_FL_beta_l1_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(mu_error_fl)), path = paste0(SGE_ID, "_", n, "_FL_mu_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(gamma_l1_error_fl)), path = paste0(SGE_ID, "_", n, "_FL_gamma_l1_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(c(actives_as_inactive_fl, inactives_as_active_fl))), path = paste0(SGE_ID, "_", n, "_FL_misclassification.csv"), append = TRUE)
  
  rm(fitted_probabilities)
  fixed_lambda_probabilities <- predict(fit, newx = X, s=lambda, type = "response")%>%as.vector()
  
  hat_gamma <- hat_theta[(n+2):(n+1+p)]
  #Z <- cbind(Z_1, Z_2)
  W_hat <- sparseMatrix(i = 1:choose(n,2), j=1:choose(n,2), 
                        x = fixed_lambda_probabilities*(1-fixed_lambda_probabilities))
  sample_gram_matrix <- 1/choose(n,2)*t(Z)%*%W_hat%*%Z
  root_Theta_hat_diag <- sample_gram_matrix%>%as.matrix()%>%solve()%>%diag()%>%sqrt()
  
  standardized_gamma_matrix <- sqrt(choose(n,2))*(hat_gamma - gamma_star)/root_Theta_hat_diag
  CI_gamma <- rbind(hat_gamma + qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)),
                    hat_gamma - qnorm(0.975)*root_Theta_hat_diag/sqrt(choose(n,2)))
  CI_lengths <- CI_gamma[1,] - CI_gamma[2,]
  gamma_1_coverage <- (gamma_star[1] <= CI_gamma[1,1]) && (gamma_star[1] >= CI_gamma[2,1])
  gamma_2_coverage <- (gamma_star[2] <= CI_gamma[1,2]) && (gamma_star[2] >= CI_gamma[2,2])
  confidence_intervals_gamma_matrix <- c(gamma_1_coverage, CI_lengths[1], gamma_2_coverage, CI_lengths[2])
  
  readr::write_csv(as.data.frame(t(c(standardized_gamma_matrix, confidence_intervals_gamma_matrix))), 
                   path = paste0(SGE_ID, "_", n, "_FL_gamma_inference.csv"), append = TRUE)
  
  
  
  # Other summaries
  
  readr::write_csv(as.data.frame(t(network_summaries)), 
                   path = paste0(SGE_ID, "_", n, "_network_summaries.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(IC_matrix)), 
                   path = paste0(SGE_ID, "_", n, "_IC_matrix.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(log_lambda_matrix)),
                   path = paste0(SGE_ID, "_", n, "_log_lambda_matrix.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(gamma_l1_error_matrix)), 
                   path = paste0(SGE_ID, "_", n, "_gamma_l1_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(gamma_l2_error_matrix)), 
                   path = paste0(SGE_ID, "_", n, "_gamma_l2_error.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(misclassification_matrix)), 
                   path = paste0(SGE_ID, "_", n, "_misclassification_matrix.csv"), append = TRUE)
  readr::write_csv(as.data.frame(t(fit$a0 - mu)), 
                   path = paste0(SGE_ID, "_", n, "_mu_error.csv"), append = TRUE)
  
  rm(fixed_lambda_probabilities, log_like_vector, adj_matrix, fit)
  
}
end_time <- Sys.time()
print(end_time - start_time)
