# Metropolis Hastings
library(MASS)
library(mvtnorm)
library(invgamma)
library(ggplot2)
library(here)
library(coda)
library(tidyverse)
library(tidybayes)
library(bayesplot)
source(here::here("R", "convert-to-coda.R"))

# Gibbs sampling

##create the data for a simple linear regression

set.seed(111)
k <- 5
n <- 100
mu_0 <- rep(5,k)
a <- 1
b <- 1
sigma_sq <- rinvgamma(a,b)
X <- cbind(1,mvrnorm(n, mu_0, diag(sqrt(sigma_sq), nrow = k, ncol = k)))
p = ncol(X)
mean_true_beta <- rep(0, 6)
true_beta <- mvrnorm(1, mean_true_beta, diag(1, nrow = p, ncol = p))
y <- X %*% true_beta + rnorm(n)

plot(y, X %*% true_beta)

# Metropolis Hastings algorithm

mcmc_MH <- function(niter, beta_tune, sigma2_tune){
    
    set.seed(222)
    
    # Prior distribution
    
    prior <- function(beta, sigma2){  
        
        prior_beta <- dmvnorm(beta, rep(0, p), diag(1, nrow = p, ncol = p), log = T) 
        prior_sigma2  <- dunif(sigma2 , min = 0, max = 30, log = T)  
        return(prior_beta + prior_sigma2)
    }
    
    # Likelihood function
    
    likelihood <- function(beta, sigma2){   
        
        mean_y = X %*% beta
        log_likelihood = dnorm(y, mean_y, sd = sqrt(sigma2), log = T)
        return(sum(log_likelihood))   
    }
    
    
    # Posterior
    
    posterior <- function(beta, sigma2){
        
        return(prior(beta, sigma2) + likelihood(beta, sigma2))
        
    }
    
   # Initialize parameters
    
    beta_chain        <- matrix(0,  nrow = niter, ncol = p)
    sigma2_chain      <- rep(0, niter)
    beta_chain[1,]    <- mvrnorm(1, rep(0, p), diag(1, nrow = p, ncol = p)) 
    sigma2_chain[1]   <- rinvgamma(1, shape = 1, rate =1)

    
    # Metropolis-Hastings step
    
    for (i in 2 : niter){ 
        
        proposed_beta <- mvrnorm(1, beta_chain[i-1,], diag(beta_tune))
        proposed_simga2 <- rnorm(1, sigma2_chain[i-1], sigma2_tune)
        
        if(proposed_simga2 < 0){
            proposed_simga2 <- sigma2_chain[i-1]

        } else {

            proposed_simga2 <- proposed_simga2
        }
        
        #alpha  <- min(1, exp(posterior(proposed_beta, proposed_simga2)  - posterior(beta_chain[i-1,], sigma2_chain[i-1])))
        alpha  <- exp(posterior(proposed_beta, proposed_simga2)  - posterior(beta_chain[i-1,], sigma2_chain[i-1]))
        # Include the acceptance ratio
        
        a <- logical(length = niter)
        
        if(runif(1) < alpha){ 
            beta_chain[i,] <- proposed_beta 
            sigma2_chain[i] <- proposed_simga2
            a[i] <- TRUE
        }else{
            
            beta_chain[i,] <- beta_chain[i-1,] 
            sigma2_chain[i] <- sigma2_chain[i-1] 
        }
        
    } 
    return(list(
        beta_chain = beta_chain,
        sigma2_chain = sigma2_chain,
        a = a
    ))
}


# mcmc_MH results

beta_tune <- rep(0.01, p)
sigma2_tune <- 0.01
niter <- 5000
out_MH = mcmc_MH(niter, beta_tune, sigma2_tune)

# Trace plots

out_MH$beta_chain
post_sigma2 <- out_MH$sigma2_chain
post_beta <- out_MH$beta_chain

plot(out_MH$sigma2_chain, type = "l" )
abline(h = sigma_sq)
matplot(out_MH$beta_chain, type = "l",  xlab = "iterations", col = c("red",  "blue", "green", "purple", "black"))
abline(h = true_beta, col = c("red",  "blue", "green", "purple", "black"))

# Accepted samples

acceptance <- out_MH$a
accepted_samples <- length(acceptance[acceptance[TRUE]])
accepted_samples

accept_ratio <- 100*(accepted_samples/niter)
accept_ratio

# Or

accept_ratio1 <- 100*mean(as.numeric(out_MH$a))
accept_ratio1




