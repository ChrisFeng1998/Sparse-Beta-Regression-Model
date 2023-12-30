
# Friendship network, symmetric friendship -------------------------------

# Corresponding website:
# https://www.stats.ox.ac.uk/~snijders/siena/Lazega_lawyers_data.htm

library(glmnet)
library(matrixStats)
library(tidyverse)
library(igraph)

# import attributes
attributes <- read.delim("data/LazegaLawyers/ELattr.dat", 
                         header = FALSE, sep = " ")
names(attributes) <- c("seniority",
                       "status (partner/assoc)",
                       "gender (man/ woman)",
                       "office",
                       "years_with_firm",
                       "age",
                       "practice",
                       "law_school")
# remove unnecessary columns
attributes <- attributes[, -c(1,9)]

# import lawyer adjacency matrix.
adj <- read.delim("data/LazegaLawyers/ELfriend.dat",
                  header = FALSE, sep = " ")
adj <- as.matrix(adj)

# Create covariate vectors ------------------------------------------------

# We want to use the covariates as in Jochmas: same status, same gender, same office,
# difference in tenure and difference in age (p.9)

# For each covariate we want to use the same structure from the simulations
n <- nrow(adj)
p <- 7
# same status
status_matrix <- matrix(nrow = n, ncol = n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    status_matrix[j,i] <- as.numeric(attributes$`status (partner/assoc)`[i] == attributes$`status (partner/assoc)`[j])
  }
}
status <- status_matrix[lower.tri(status_matrix)]
rm(status_matrix)

# same gender
gender_matrix <- matrix(nrow = n, ncol = n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    gender_matrix[j,i] <- as.numeric(attributes$`gender (man/ woman)`[i] == attributes$`gender (man/ woman)`[j])
  }
}
gender <- gender_matrix[lower.tri(gender_matrix)]
rm(gender_matrix)

# same office
office_matrix <- matrix(nrow = n, ncol = n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    office_matrix[j,i] <- as.numeric(attributes$office[i] == attributes$office[j])
  }
}
office <- office_matrix[lower.tri(office_matrix)]
rm(office_matrix)

# difference in tenure (diff. years with firm)
tenure_matrix <- matrix(nrow = n, ncol = n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    tenure_matrix[j,i] <- abs(attributes$years_with_firm[i] - attributes$years_with_firm[j])
  }
}
tenure <- tenure_matrix[lower.tri(tenure_matrix)]
rm(tenure_matrix)

# difference in age
age_matrix <- matrix(nrow = n, ncol = n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    age_matrix[j,i] <- abs(attributes$age[i] == attributes$age[j])
  }
}
age <- age_matrix[lower.tri(age_matrix)]
rm(age_matrix)

# same practice
practice_matrix <- matrix(nrow = n, ncol = n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    practice_matrix[j,i] <- as.numeric(attributes$practice[i] == attributes$practice[j])
  }
}
practice <- practice_matrix[lower.tri(practice_matrix)]
rm(practice_matrix)

# same law school
school_matrix <- matrix(nrow = n, ncol = n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    school_matrix[j,i] <- as.numeric(attributes$law_school[i] == attributes$law_school[j])
  }
}
school <- school_matrix[lower.tri(school_matrix)]
rm(school_matrix)
rm(i,j)


# Deterministic desgin matrix ---------------------------------------------

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
rm(a,b,c, i, rows, cols)
X <- cbind(X_deterministic, status, gender, office, tenure, age, practice, school)
rm(X_deterministic, status, gender, office, tenure, age, practice, school)


# Generate symmetric adjacency matrix of mutual friendships
sym_adj <- graph_from_adjacency_matrix(adjmatrix = adj)%>%
  as.undirected(mode = "mutual")%>%
  get.adjacency()%>%
  as.matrix()
degrees <- rowSums(sym_adj)
print(sprintf("Maximum degree: %s", max(degrees))) # 16
print("Degree distribution")
print(summary(degrees)) # median = 5
# edge density of 0.07
print("Edge density:")
print((sum(sym_adj)/2)/choose(n,2))

Y <- sym_adj[lower.tri(sym_adj)]
fit <- glmnet(X,Y, family = "binomial", alpha = 1, penalty.factor = c(rep(1,n), rep(0,p)),
              lower.limits = c(rep(0,n), rep(-Inf, p)))

# Calculate BIC
log_like <- function(adj, probs){
  return(sum(adj%*%log(probs) + (1-adj)%*%log(1 - probs)))
}
# log likelihood
fitted_probabilities <- predict(fit, newx = X, type = "response")
# This is the POSITIVE log likelihood for each fitted model.
log_like_vector <- apply(fitted_probabilities, 
                         MARGIN = 2, 
                         function(i) log_like(adj = sym_adj[lower.tri(sym_adj)], probs = i))

bic <- (fit$df+1)*log(choose(n,2)) - 2*log_like_vector

bic_best_pos <- which.min(bic) # 9
fitted_values <- coef(fit, s = fit$lambda[bic_best_pos])
print(fitted_values)
# remove intercept and covariates
betas <- as.vector(coef(fit, s = fit$lambda[bic_best_pos]))[-c(1,73:79)]

print(sprintf("Number of non-zeor betas: %s", length(betas[betas!=0]))) # 6
print("Non-zero beta summary")
print(summary(betas[betas!=0]))

# Inference ---------------------------------------------------------------
N <- n*(n-1)/2
Z <- X[,72:78]
hat_gamma <- fit$beta[(n+1):(n+p),bic_best_pos]
names(hat_gamma) <- c("same_status",
                      "same_gender",
                      "same_office",
                      "diff_years_with_firm",
                      "diff_age",
                      "same_practice",
                      "same_law_school")

W_hat <- sparseMatrix(i = 1:N, j=1:N, x = fitted_probabilities[,bic_best_pos]*(1-fitted_probabilities[,bic_best_pos]))
sample_gram_matrix <- 1/N*t(Z)%*%W_hat%*%Z
root_Theta_hat_diag <- sample_gram_matrix%>%as.matrix()%>%solve()%>%diag()%>%sqrt()
CI_hat <- rbind(hat_gamma + qnorm(0.975)*root_Theta_hat_diag/sqrt(N),
                hat_gamma - qnorm(0.975)*root_Theta_hat_diag/sqrt(N))
CI_lengths <- CI_hat[1,] - CI_hat[2,]

print("Coefficients")
print(hat_gamma)

print("Confidence interval")
print(t(CI_hat))

# Plotting ----------------------------------------------------------------

g <- graph_from_adjacency_matrix(adjmatrix = sym_adj, mode = "undirected")
set.seed(123)
g_force <- layout_with_fr(g)
sizes <- ifelse(degree(g) > 5, degree(g), 5)  

colors <- attributes%>%
  transmute(color = ifelse(office ==1,  "#56B4E9",
                        ifelse(office == 2, "#E69F00", "black")))%>%
  pull()
  
  
# Network by office
plot(g, vertex.size = sizes, vertex.label = NA, layout = g_force,
     vertex.color = colors, edge.width = 0.8)

colors_status <- attributes%>%
  transmute(color_status = ifelse(`status (partner/assoc)` ==1, "firebrick", "forestgreen"))%>%
  pull()
# Network by status
plot(g, vertex.size = sizes, vertex.label = NA, layout = g_force,
     vertex.color = colors_status, edge.width = 0.8)
