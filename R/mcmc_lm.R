library(mvnfast)
library(ggplot2)
library(truncnorm)
library(coda)
library(invgamma)

set.seed(100)
# n = 100
# p <- 5 # length of truebeta
# truesigma2 = 4
# truebeta <- c(rmvn(1, rep(2, p), diag(20, nrow = p, ncol = p)))
# X <- cbind(1, rmvn(n, rep(0, p-1), diag(1, nrow = p-1, ncol = p-1)))
# y <- X %*% truebeta + rnorm(n, mean = 0,  sd = sqrt(truesigma2))
# 
# plot(y, X %*% truebeta, main = " Response vs expected response", ylab = "expecte response", cex.lab = 1, cex.main=1)

mcmc_lm <- function(niter) {
    
    #Initial values
    
    beta <- matrix(0, niter, p)
    sigma2 <- rep(0, niter)
    sigma2[1] <- 1
    
    # sample beta 
    
    for(i in 2 : niter){
        
        mu <- solve(t(X) %*% (X)) %*% (t(X) %*% (y)) # mean of beta
        Sigma <- solve(t(X) %*% (X)) * sigma2[i-1] # covariance matrix for beta
        beta[i,] <- rmvn(1, mu, Sigma)
        
        # sample sigma_sq
        
        b_n <- 0.5 * t(y - X %*% beta[i,]) %*% (y - X %*% beta[i,]) + 1
        a_n <- (n/2) + 1
        sigma2[i] <- rinvgamma(1, a_n, rate = b_n)
        
    }
    return(
        list(
            beta = beta, 
            sigma2 = sigma2
        )
        
    )
    
    
    
}

# 
# 
# out_lm <- mcmc_lm(niter = 10000) 
# 
# colors <- c("blue", "red", "green", "purple", "darkblue", "black")
# matplot(out_lm[[1]], type  = "l", col = colors)
# abline(h = truebeta, col = colors)
# 
# 
# save(out_lm,  file = here::here("results", "out_lm.RData"))
# 
# 



