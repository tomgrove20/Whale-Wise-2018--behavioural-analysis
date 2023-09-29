# https://github.com/fabsig/GPBoost/blob/master/R-package/demo/classification_non_Gaussian_data.R


install.packages("gpboost")
library(gpboost)

# specifying likelihood
likelihood <- "bernoulli_logit"

params <- list(learning_rate = 0.1, min_data_in_leaf = 20,
               objective = "binary", monotone_constraints = c(0,1))

nrounds <- 50

# defining random variables
group = as.matrix(sapply(tot[,c("surface.feeding")], as.numeric))
group = as.vector(tot$folnum.unique)

# defining data
X = as.matrix(sapply(tot[,vars[c(7,8)]], as.numeric))  
y = as.vector(tot$surface.active)
data_train = gpb.Dataset(X, y)

gp_model <- GPModel(group_data = group, likelihood = likelihood)
bst <- gpboost(data = X, label = y, 
               group_data = gp_model, 
               nrounds = nrounds, 
               params = params)
summary(bst) # Trained random effects model (true variance = 0.5)

cvbst <- gpb.cv(params = params,
                data = X, label = y,
                gp_model = gp_model,
                nrounds = 200,
                nfold = 4,
                verbose = 1,
                early_stopping_rounds = 5,
                use_gp_model_for_validation = TRUE)
print(paste0("Optimal number of iterations: ", cvbst$best_iter))


#--------------------Choosing tuning parameters----------------
param_grid = list("learning_rate" = c(1,0.1,0.01), "min_data_in_leaf" = c(1,10,100),
                  "max_depth" = c(1,3,5,10))
gp_model <- GPModel(group_data = group, likelihood = likelihood)
dataset <- gpb.Dataset(data = X, label = y)
set.seed(10)
opt_params <- gpb.grid.search.tune.parameters(param_grid = param_grid,
                                              params = params,
                                              num_try_random = NULL,
                                              nfold = 4,
                                              data = dataset,
                                              gp_model = gp_model,
                                              verbose_eval = 1,
                                              nrounds = 1000,
                                              early_stopping_rounds = 10)
print(paste0("Best parameters: ",paste0(unlist(lapply(seq_along(opt_params$best_params), function(y, n, i) { paste0(n[[i]],": ", y[[i]]) }, y=opt_params$best_params, n=names(opt_params$best_params))), collapse=", ")))
print(paste0("Best number of iterations: ", opt_params$best_iter))
print(paste0("Best score: ", round(opt_params$best_score, digits=3)))





## Examples for combining tree-boosting with Gaussian process and random effects models
##    for several non-Gaussian likelihoods
## See the examples in GPBoost_algorithm.R for more functionality
## Author: Fabio Sigrist

library(gpboost)

## Choose likelihood: either "bernoulli_probit" (=default for binary data), "bernoulli_logit",
##                       "poisson", or "gamma"
likelihood <- "bernoulli_probit"

# Non-linear prior mean function for simulation in examples below
f1d <- function(x) 1/(1+exp(-(x-0.5)*10)) - 0.5
sim_non_lin_f <- function(n){
  X <- matrix(runif(2*n),ncol=2)
  f <- f1d(X[,1])
  return(list(X=X,f=f))
}

# Parameters for gpboost in examples below 
# Note: the tuning parameters are by no means optimal for all situations considered here
params <- list(learning_rate = 0.1, min_data_in_leaf = 20,
               objective = likelihood, monotone_constraints = c(1,0))
nrounds <- 25
if (likelihood=="bernoulli_logit") nrounds <- 50
if (likelihood %in% c("bernoulli_probit","bernoulli_logit")) params$objective="binary"

#--------------------Combine tree-boosting and grouped random effects model----------------
# Simulate data
n <- 5000 # number of samples
m <- 500 # number of groups
set.seed(1)
# Simulate random and fixed effects
group <- rep(1,n) # grouping variable
for(i in 1:m) group[((i-1)*n/m+1):(i*n/m)] <- i
b1 <- sqrt(0.5) * rnorm(m)
eps <- b1[group]
eps <- eps - mean(eps)
sim_data <- sim_non_lin_f(n=n)
f <- sim_data$f
X <- sim_data$X
# Simulate response variable
if (likelihood == "bernoulli_probit") {
  probs <- pnorm(f+eps)
  y <- as.numeric(runif(n) < probs)
} else if (likelihood == "bernoulli_logit") {
  probs <- 1/(1+exp(-(f+eps)))
  y <- as.numeric(runif(n) < probs)
} else if (likelihood == "poisson") {
  mu <- exp(f+eps)
  y <- qpois(runif(n), lambda = mu)
} else if (likelihood == "gamma") {
  mu <- exp(f+eps)
  y <- qgamma(runif(n), scale = mu, shape = 1)
}
hist(y,breaks=50)# visualize response variable

#--------------------Training----------------
# Define random effects model
gp_model <- GPModel(group_data = group, likelihood = likelihood)
bst <- gpboost(data = X, label = y, verbose = 0,
               gp_model = gp_model,
               nrounds = nrounds, 
               params = params)
summary(gp_model) # Trained random effects model (true variance = 0.5)
