library(coda)
# library(mvtnorm)
library(mvnfast)
library(MASS)
source(here::here("R", "update_tuning_example.R"))

n = 100
p = 4
truebeta <- rnorm(p, 0, 2)
sigma2 = 1

## very slow # x <-  mvrnorm(n, rep(0, p), diag(1, p))
x <-  rmvn(n, rep(0, p), diag(1, p))
y =   x %*% truebeta + rnorm(n)

likelihood = function(param, x, y){
  beta = param
  # sigma2 = param[2]
  ## univariate likelihood
  ll <- sum(dnorm(y, mean = x %*% beta, sd = sqrt(sigma2), log = T))
  return(ll)
}

prior = function(param){
  beta = param
  p <- length(beta)
  ## slow # prior_beta <- dmvnorm(beta, rep(0, p), diag(10, p), log = T)
  ## ok
  # prior_beta <- dmvn(beta, rep(0, p), diag(10, p), log = T)
  ## faster
  prior_beta <- dmvn(beta, rep(0, p), diag(sqrt(10), p), isChol = T, log = T)
  return(prior_beta)
}

posterior <- function(param, x, y){
  return(prior(param) + likelihood(param, x, y))
}

Adapt_MH_MCMC_lm <- function(K, Sigma_tune, Sigma_tune_chol){
  bathch_idx <- seq(1, K, by = 50)
  mcmc_sample <- matrix(0, K, p)
  mcmc_sample[1,] = rep(2, p)
  accept <- 0
  accept_adapt <- 0
  lambda <- 5
  
  for (i in 2:K){
    #s <- rep(1, p)*lambda
    current_mcmc_sample = mcmc_sample[i-1,]
    
    # proposal <-  mvrnorm(1, current_mcmc_sample, Sigma_tune_chol ) 
    proposal <- as.vector(rmvn(1, current_mcmc_sample, sigma = lambda * Sigma_tune_chol, isChol = TRUE))
    
    
    #proposal <-  mvrnorm(1, current_mcmc_sample, diag(lambda, p) %*% Sigma_tune_chol %*% diag(lambda, p) ) 

   

    # if (proposal[p] < 0){
    #   proposal[p] <- -proposal[p]
    # } else {
    #   proposal <- proposal
    # }
    

    probab = exp(posterior(proposal, x, y) - posterior(current_mcmc_sample, x, y))
    
    if (runif(1) < probab){
      mcmc_sample[i,] = proposal
      accept <- accept + 1/K
      accept_adapt <- accept_adapt + 1/ 50
    }else{
      mcmc_sample[i,] = current_mcmc_sample
    }
      # update step
    
    if(i %in% bathch_idx){
      batch_samples <- mcmc_sample[(i-50) : i, ]
      #accept <- sum(accept[(i-50) : i])/nrow(batch_samples)
      param_update <- update_tuning_mv(i, accept_adapt, lambda, batch_samples,
                                       Sigma_tune, Sigma_tune_chol)

      Sigma_tune_chol <- param_update$Sigma_tune_chol
      lambda          <- param_update$lambda
      Sigma_tune      <- param_update$Sigma_tune
      accept_adapt    <- param_update$accept
    }
  }
  
  message("Acceptance rate = ", accept)
  
  return(mcmc(mcmc_sample))
}

# Initial values

Sig_chol <- matrix(0, p, p)
Sig_chol[lower.tri(Sig_chol, diag = TRUE)] <- 1

mcmc_out1 = Adapt_MH_MCMC_lm(K = 5000, Sigma_tune = diag(1,p), Sig_chol)
str(mcmc_out1)
matplot(mcmc_out1, type = "l", col = c("black", "red",  "blue", "green"))
abline(h = truebeta, col = c("black", "red",  "blue", "green"))

plot(apply(mcmc_out1, 2, mean), truebeta)
abline(0, 1, col = "red")

effectiveSize(mcmc_out1)

summary(mcmc_out1)




plot(chain)
plot(chain1)
## Test update_tuning_mv()

m <- matrix(0, 50, 4)
m[, 1:4] <- rnorm(50)
m
lambda <- rnorm(1)
update_tuning_mv(51, 0.5, lambda, m, diag(1, 4),  diag(1, 4))




mat <- matrix(1:9, 3, 3)
mat[lower.tri(mat)] <- 0
mat
