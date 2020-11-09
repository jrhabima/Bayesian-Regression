library(mvnfast)
library(ggplot2)
library(truncnorm)
library(coda)

set.seed(111)
n = 100
p <- 5 # length of truebeta
truesigma2 = 4
truebeta <- c(rmvn(1, rep(2, p), diag(20, nrow = p, ncol = p)))
X <- cbind(1, rmvn(n, rep(50, p-1), diag(1, nrow = p-1, ncol = p-1)))
y <- X %*% truebeta + rnorm(n, mean = 0,  sd = sqrt(truesigma2))

plot(y, X %*% truebeta, main = " Response vs expected response", ylab = "expecte response", cex.lab = 1, cex.main=1)

#mcmc_lm_mh <- function(niter, beta_tune, simga2_tune){

# Prior distribution

prior = function(param){
  beta   = param[1:p]
  sigma2 = param[p + 1]
  
  prior_beta <-  dmvn(beta, rep(0, p), diag(10, nrow = p, ncol = p), log = T) 
  prior_sigma2 <- dunif(sigma2, min=0, max=30, log = T)
  return(prior_beta +  prior_sigma2)
}

# Likelihood function

likelihood = function(X, y, param){
  beta   = param[1:p]
  sigma2 = param[p + 1]
  
  l = sum(dnorm(y, mean = X %*% beta, sd = sqrt(sigma2), log = T))
  return(l)
}

# Posterior distribution

posterior <- function(X, y, param){
  return(prior(param) + likelihood(X, y,param))
}

mcmc_multi_lm <- function(X, y, K, proposal_sd){
  chain = array(dim = c(K, p+1))
  chain[1,] = rep(2, p+1)
  accept <- 0
  for (i in 2:K){
    current_chain = chain[i-1,]
    #proposal <- rmvn(p+2, mean = current_chain, sd = proposal_sd)
    proposal <- as.vector(rmvn(1, current_chain, sigma = diag(proposal_sd, p+1), isChol = TRUE))
    
    if (proposal[p+1] < 0){
      proposal[p+1] <- -proposal[p+1]
    } else {
      proposal <- proposal
    }
    probab = exp(posterior(X, y, proposal) - posterior(X, y, current_chain))
    
    if (runif(1) < probab){
      chain[i,] = proposal
      accept <- accept + 1 / K
    }else{
      chain[i,] = current_chain
    }
  }
  message("acceptance ratio is", accept)
  return(mcmc(chain))
}


chain1 = mcmc_multi_lm(X = X, y = y, 10000, rep(10, p+1))
chain2 = mcmc_multi_lm(X = X, y = y, 10000, rep(1, p+1))
chain3 = mcmc_multi_lm(X = X, y = y, 10000, rep(0.1, p+1))

png(file = here::here("figures", "MH_tuning_choices.png"),  #Creating a plotting window in folder
    res = 400, height = 6, width = 3, units = "in") # This number can be adjusted for a better view

par(mfrow =c(3, 1))

colors <- c("blue", "red", "green", "purple", "darkblue", "black", "violet")
matplot(chain1, type = "l", main = "proposal sd = 10", col = colors,
        cex.main = 1, xlab = "iterations", ylab = "parameter value")
#abline(h = c(truebeta0,  truebeta,  truesigma2), col = colors)

matplot(chain2, type = "l", main = "proposal sd = 1", col = colors,
        cex.main = 1, xlab = "iterations", ylab = "parameter value")
#abline(h = c(truebeta0,  truebeta,  truesigma2), col = colors)

matplot(chain3, type = "l", main = "proposal sd = 0.1", col = colors,
        cex.main = 1, xlab = "iterations", ylab = "parameter value")
#abline(h = c(truebeta0,  truebeta,  truesigma2), col = colors)

dev.off() #Save the graph to the window 

# colors <- c("blue", "red", "green", "purple", "darkblue", "black")
# chain4 = mcmc_multi_lm(X = X, y = y, 100000, c(2, rep(0.25, p+1)))
# matplot(chain4, type = "l", main = "proposal sd = 0.1", col = colors,
#         cex.main = 1, xlab = "iterations", ylab = "parameter value", ylim = c(-10, 10))
# abline(h = c(truebeta,  truesigma2), col = colors)


