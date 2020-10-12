# # Load packages
# 
# library(MASS)
# library(mvtnorm)
# library(invgamma)
# library(ggplot2)
# library(here)
# library(coda)
# library(tidyverse)
# library(tidybayes)
# library(bayesplot)
# source(here::here("R", "convert-to-coda.R"))
# 
# 
# # Gibbs sampling

## Create the data to train the model





# Codes


 
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

# # mcmc output
# 
# N = 5000
# 
# if(!file.exists(here::here("results", "out_lm"))){
#     
#    ptm <- proc.time() # timing
#    
#     out_lm <- mcmc_lm(N) 
#     save(out_lm,  file = here::here("results", "out_lm"))
#     
#    proc.time() - ptm
# }
# 
# load(here::here("results", "out_lm"))
# 
# 
# 
# # Trace plots
# 
# out_coda %>%
#     mcmc_trace(regex_pars = c("beta"),
#                facet_args = list(nrow = 3,  ncol = 3))
# abline(h = true_beta, col = c("red",  "blue", "green", "purple", "black", "green"))
# 
# out_coda %>%
#     mcmc_trace(regex_pars = c("sigma_sq"))
#             
# ### Other plots illustrating posterior distribution for our model parameters.
# 
# plot_title1 <- ggtitle("Credible intervals for model parameter")
# plot_title2 <- ggtitle("Credible intervals for regression coefficints")
# title3 <- ggtitle("Posterior density for regression coefficients")
# 
# out_coda %>%
#     mcmc_intervals(regex_pars = c("beta", "sigma_sq")) + plot_title1
# out_coda %>%
#     mcmc_areas(regex_pars = c("beta"),
#                prob = 0.8, # 80% intervals
#                prob_outer = 0.99, # 99%
#                point_est = "mean") + plot_title2
# 
# color_scheme_set("blue")
# mcmc_dens(out_coda, regex_pars = c("beta")) + title3
# 
# ## Regression model diagnostics
# 
# beta_post_mean <- apply(out_lm$beta, 2, mean) # Posterior mean of beta
# 
# y_pred <- X %*% beta_post_mean
# 
# residuals <- y - y_pred
# 
# # Fitted density vs true density
# 
# plot(density(y), col = "red", main = "Fitted density vs observed density")
# lines(density(y_pred), col = "blue")
# legend("topright", legend = c("true density",  "fitted density"),
#        col = c("red", "blue"), lty=1:2, cex=0.8)
# 
# # Residual and fitted values plots
# 
# par(mfrow = c(2,2))
# hist(residuals, main = "residual histogram")
# qqnorm(residuals, main = "residual qqplot")
# abline(c(0, 1), col ="blue")
# plot( y_pred, residuals, main = "resuduals vs fitted values")
# plot(y, y_pred, main = "Observed values vs predicted values")
# abline(c(0, 1), col ="blue") ## a line of 45 degree angle
