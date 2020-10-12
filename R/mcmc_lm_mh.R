
mcmc_lm_mh <- function(niter, beta_tune, sigma2_tune){
    
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
        
        accep_ratio  <- exp(posterior(proposed_beta, proposed_simga2) - posterior(beta_chain[i-1,], sigma2_chain[i-1]))
              alpha  <- min(1, accep_ratio)
              
        # Include the acceptance ratio
        
        a <- logical(length = niter)
        
        
        if (runif(1) < alpha) {
            beta_chain[i,] <- proposed_beta
            sigma2_chain[i] <- proposed_simga2
            a[i] <- TRUE
  } else {
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

