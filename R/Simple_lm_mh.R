
# Generate the data

truebeta0 = 0
truebeta1 = 5
truesigma2 = 10
n = 100

x <- rnorm(n, 5, 10)
y =  truebeta0 + truebeta1* x + rnorm(n, mean = 0, sqrt(truesigma2))


likelihood = function(param){
  beta0  = param[1]
  beta1  = param[2] 
  sigma2 = param[3]
  
  mu = beta0 + beta1*x 
  likelihood = sum(dnorm(y, mean = mu, sd = sqrt(sigma2), log = T))
  return(likelihood)
}

likelihood(startvalue)

prior = function(param){
  beta0  = param[1]
  beta1  = param[2] 
  sigma2 = param[3]
  
  prior_beta0 = dnorm(beta0, sd = 5, log = T)
  prior_beta1 = dunif(beta1, min=0, max=10, log = T)
  prior_sigma2 = dunif(sigma2, min=0, max = 30, log = T)
  return(prior_beta0 + prior_beta1 + prior_sigma2)
}

prior(startvalue)

posterior <- function(param){
  return(likelihood(param) + prior(param))
}

posterior(startvalue)

mcmc_lm = function(startvalue, tuning, K){
  chain = matrix(0, nrow = K, ncol = 3)
  chain[1,] = startvalue
  for (i in 2:K){
    current_chain <- chain[i-1, ]
    
    proposed_chain  = rnorm(3, mean = current_chain, sd = tuning )
    
    if (proposed_chain[3] < 0){
      proposed_chain[3] <- -proposed_chain[3]
    } else {
      proposed_chain <-  proposed_chain
    }
    
    a = exp(posterior(proposed_chain) - posterior(current_chain))
    
    if (runif(1) < a){
      chain[i,] = proposed_chain
    }else{
      chain[i,] = current_chain
    }
  }
  return(chain)
}

tuning <- c(0.5, 0.5, 0.5)
startvalue = rep(4, 3)

chain = mcmc_lm(startvalue, tuning, 10000)
matplot(chain, type = "l", col = c("red", "blue", 'orange'))
abline(h = c(truebeta0 , truebeta1,truesigma2), col = c("red", "blue", 'orange'))
