library(glmnet)
library(matrixStats)
library(tidyverse)
library(igraph)

library(readxl)

# From "The log of Gravity" paper, http://personal.lse.ac.uk/tenreyro/LGW.html

Log_of_Gravity <- read_excel("data/trade_network/Log of Gravity.xls")
countrycodes <- read_excel("data/trade_network/countrycodes.xls")
names(countrycodes)[1] <- "code"


# Significant trade partnerships -------------------------------------------------

# Find signigicant trade partnerships based on trade volume
# Calculate GDP per country
gdps <- countrycodes%>%
  right_join(Log_of_Gravity%>%
               select(s1_im, lypim)%>%
               group_by(s1_im)%>%
               summarise(log_gdp = mean(lypim))%>%
               mutate(importer_GDP = exp(log_gdp))%>%
               select(s1_im, importer_GDP),
             by = c("code" = "s1_im"))%>%
  select(-(3:4))

Total_import <- Log_of_Gravity%>%
  group_by(s1_im)%>%
  summarise(Total_import = sum(trade))
Total_export <- Log_of_Gravity%>%
  group_by(s2_ex)%>%
  summarise(Total_export = sum(trade))
Total_import_export <- Total_import%>%
  right_join(Total_export, by = c("s1_im" = "s2_ex"))%>%
  right_join(countrycodes, by = c("s1_im" = "code"))%>%
  select(c(1,4,2,3))
names(Total_import_export)[1] <- "code"

# generate edge list: an edge is placed between countries with
# a significant trade partnership as defined in the paper
Network_data <- Log_of_Gravity%>%
  select(s1_im, s2_ex, trade)%>%
  right_join(countrycodes[,1:2], by = c("s1_im" = "code"))%>%
  rename(importer_country = country)%>%
  right_join(countrycodes[,1:2], by = c("s2_ex" = "code"))%>%
  rename(exporter_country = country)%>%
  right_join(Total_import, by = "s1_im")%>%
  right_join(Total_export, by = "s2_ex")%>%
  select(s1_im, s2_ex, importer_country, exporter_country,
         trade, Total_import, Total_export)%>%
  mutate(relevant_import = ifelse(trade >= 0.03*Total_import, 1, 0),
         relevant_export = ifelse(trade >= 0.03*Total_export, 1, 0))%>%
  filter((relevant_import == 1) | (relevant_export == 1))%>%
  select(importer_country, exporter_country)%>%
  arrange(importer_country)


n <- nrow(countrycodes)
library(ggraph)
g <- graph_from_data_frame(Network_data, directed = FALSE, 
                           vertices = countrycodes$country)%>%
  igraph::simplify()
print(sprintf("Total number of edges: %s", length(E(g)))) # 1279
print(sprintf("Edge density: %s", edge_density(g))) # 0.1393246

adj_matrix <- as_adjacency_matrix(g, type = "both")%>%as.matrix()
isSymmetric(adj_matrix)
degrees <- rowSums(adj_matrix)

# Covariates --------------------------------------------------------------

# use same covariates as Jochmans: log-distance, common border, common language,
# colonial ties, preferential trade agreement

# Nodes in g are ordered in alphabetical order, just as in countrycodes

covariates <- Log_of_Gravity%>%select(s1_im, s2_ex, ldist, border, comlang, colony, comfrt)%>%
  filter(FALSE)
for (i in 1:(n-1)) {
  if(i %% 50 == 0){print(i)}
  country_i <- countrycodes$code[i]
  for (j in (i+1):n) {
    country_j <- countrycodes$code[j]
    covariates <- rbind(covariates,
                        Log_of_Gravity%>%select(s1_im, s2_ex, ldist, border, comlang, colony, comfrt)%>%
                          filter(s1_im == country_i, s2_ex == country_j))
  }
}
rm(country_i, country_j, i, j)

# Model fit ---------------------------------------------------------------

p <- 5
log_like <- function(adj, probs){
  return(sum(adj%*%log(probs) + (1-adj)%*%log(1 - probs)))
}
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

# model fit
Y <- adj_matrix[lower.tri(adj_matrix)]
X <- cbind(X_deterministic, as.matrix(covariates[,-c(1:2)]))
fit <- glmnet(X,Y, family = "binomial", alpha = 1, penalty.factor = c(rep(1,n), rep(0,p)),
              lower.limits = c(rep(0,n), rep(-Inf, p)))


fitted_probabilities <- predict(fit, newx = X, type = "response")
# This is the POSITIVE log likelihood for each fitted model.
log_like_vector <- apply(fitted_probabilities, 
                         MARGIN = 2, 
                         function(i) log_like(adj = adj_matrix[lower.tri(adj_matrix)], probs = i))

bic <- (fit$df+1)*log(choose(n,2)) - 2*log_like_vector

bic_best_pos <- which.min(bic)

betas <- fit$beta[, bic_best_pos]
print(sprintf("Number of non-zero betas: %s", sum(betas[1:n] != 0))) # 32 betas unequal to zero for BIC
# covariate coefficients
print("Estimated covariate coefficients")
print(betas[(n+1):(n+p)])
# data frame containing country name, estimated beta, degree and GDP
data.frame(countrycodes, beta = betas[1:n], degree = degrees)%>%
  right_join(gdps, by = "country")%>%
  select(country, beta, degree, importer_GDP)%>%
  arrange(desc(beta))%>%View()

# Inference ---------------------------------------------------------------
N <- n*(n-1)/2
Z <- X[,137:141]
hat_gamma <- fit$beta[(n+1):(n+p),bic_best_pos]
names(hat_gamma) <- c("log_distance",
                      "common_border",
                      "common_language",
                      "colonial_ties",
                      "preferential_trade_agreement")

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
