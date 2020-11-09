source(here::here("R", "libraries.R"))
library(truncnorm)
set.seed(111)
n <- 1000
k <- 5
true_sigma_sq <- 4
X <- cbind(1,mvrnorm(n, rep(50,k), diag(50, nrow = k, ncol = k)))
p = ncol(X)
true_beta <- mvrnorm(1, rep(0, p), diag(4, nrow = p, ncol = p))
y <- X %*% true_beta + rnorm(n, mean = 0,  sd= sqrt(true_sigma_sq))
plot(y, X %*% true_beta, main = " Response vs expected response", ylab = "expecte response", cex.lab = 1, cex.main=1)

set.seed(222)

mcmc_lm_mh <- function(niter, beta_tune, simga2_tune){
    
    # Prior distribution
    
    prior <- function(beta, sigma2){  
        
        prior_beta <- dmvnorm(beta, rep(0, p), diag(10, nrow = p, ncol = p), log = T) 
        prior_sigma2  <- dunif(sigma2 , min = 0, max = 50, log = T)  
        return(prior_beta + prior_sigma2)
      
    }

    # Likelihood function
    
    likelihood <- function(beta, sigma2){   
        
        log_likelihood = dnorm(y, mean = X %*% true_beta, sd = sqrt(sigma2), log = T)
        return(sum(log_likelihood))   
    }
    
   # Posterior distribution
    
    posterior <- function(beta, sigma2){
        
        return(prior(beta, sigma2) + likelihood(beta, sigma2))
        
    }
    
    # Initialize parameters
    
    beta_chain        <- matrix(0,  nrow = niter, ncol = p)
    sigma2_chain      <- rep(0, niter)
    beta_chain[1,]    <- mvrnorm(1, rep(0, p), diag(1, nrow = p, ncol = p)) 
    sigma2_chain[1]   <- 1
    
    
    # Metropolis-Hastings step
    
    for (i in 2 : niter){ 
        proposed_beta <- mvrnorm(1, beta_chain[i-1,], diag(beta_tune))
        proposed_simga2 <- rnorm(1, sigma2_chain[i-1], sigma2_tune)
        #proposed_simga2 <- rtruncnorm(1, a = 0, b=Inf, mean = sigma2_chain[i-1], sigma2_tune)

         if(proposed_simga2 < 0){
             proposed_simga2 <- sigma2_chain[i-1]

         } else {

             proposed_simga2 <- proposed_simga2
         }
        accep_ratio  <- exp(posterior(proposed_beta, proposed_simga2) - posterior(beta_chain[i-1,], sigma2_chain[i-1]))
       
        if (runif(1) <  accep_ratio) {
            beta_chain[i,] <- proposed_beta
            sigma2_chain[i] <- proposed_simga2
            
        } else {
            beta_chain[i,] <- beta_chain[i-1,]
            sigma2_chain[i] <- sigma2_chain[i-1]
        }
        
    } 
    return(list(
        beta_chain = beta_chain,
        sigma2_chain = sigma2_chain
        #a = a
    ))
}


# mcmc_lm_mh results

beta_tune <- rep(0.2, p)
sigma2_tune <- 0.5
niter <- 50000

out_mh = mcmc_lm_mh(niter, beta_tune)

matplot(out_mh$beta_chain, type = "l")
plot(out_mh$sigma2_chain, type = "l")
abline(h = true_sigma_sq,  col = "red")
