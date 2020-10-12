
mcmc_lm <- function(niter) {
    
    set.seed(100)
    #Initial values
    
    beta <- matrix(0, niter, p)
    sigma_sq <- rep(0, niter)
    sigma_sq[1] <- 1
    
    # sample beta 
    
    for(i in 2 : niter){
        
        mu <- solve(t(X) %*% (X)) %*% (t(X) %*% (y)) # mean of beta
        Sigma <- solve(t(X) %*% (X)) * sigma_sq[i-1] # covariance matrix for beta
        beta[i,] <- mvrnorm(1, mu, Sigma)
        
        # sample sigma_sq
        
        b_n <- 0.5 * t(y - X %*% beta[i,]) %*% (y - X %*% beta[i,]) + 1
        a_n <- (n/2) + 1
        sigma_sq[i] <- rinvgamma(1, a, rate = b)
        
    }
    return(
        list(
            beta = beta, 
            sigma_sq = sigma_sq
        )
        
    )
    
    
    
}
