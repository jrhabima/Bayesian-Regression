library(coda)
library(mvnfast)
library(here)
source(here::here("R", "update_tuning_example.R"))

n = 100
p <- 5
#truebeta0 = 2
truesigma2 = 4

truebeta <- c(rmvn(1, rep(0, p), diag(10, p)))
X <-  cbind(1, rmvn(n, rep(0, p-1), diag(1, p-1)))
y =  X %*% truebeta + rnorm(n, mean = 0,  sd = sqrt(truesigma2))

likelihood = function(param, X, y){
  beta  = param[1:p]
  sigma2 = param[p+1]
  ## univariate likelihood
  ll <- sum(dnorm(y, mean = X %*% beta, sd = sqrt(sigma2), log = T))
  return(ll)
}

prior = function(param){
  beta  = param[1:p]
  sigma2 = param[p+1]
  p <- length(beta)
  ## slow # prior_beta <- dmvnorm(beta, rep(0, p), diag(10, p), log = T)
  ## ok
  # prior_beta <- dmvn(beta, rep(0, p), diag(10, p), log = T)
  ## faster
  prior_beta <- dmvn(beta, rep(0, p), diag(sqrt(10), p), isChol = T, log = T)
  return(prior_beta)
}

posterior <- function(param, X, y){
  return(prior(param) + likelihood(param, X, y))
}

#posterior(rep(2, 6), X, y)
Adapt_MH_MCMC_lm <- function(X = X, y = y, K, Sigma_tune, Sigma_tune_chol, 
                             n_burn = K / 2, update_tune = FALSE){
  bathch_idx <- seq(1, K, by = 50)
  mcmc_sample <- matrix(0, K, p+1)
  mcmc_sample[1,] = rep(2, p+1)
  accept <- 0
  accept_adapt <- 0
  lambda <- 5
  
  for (k in 2:K){
    #s <- rep(1, p)*lambda
    current_mcmc_sample = mcmc_sample[k-1,]
    
    # proposal <-  mvrnorm(1, current_mcmc_sample, Sigma_tune_chol ) 
    proposal <- as.vector(rmvn(1, current_mcmc_sample, lambda * Sigma_tune_chol, isChol = TRUE))
    
    if (proposal[p+1] < 0){
      proposal[p+1] <- -proposal[p+1]
    } else {
      proposal <- proposal
    }
    
    
    probab = exp(posterior(proposal, X, y) - posterior(current_mcmc_sample, X, y))
    
    if (runif(1) < probab){
      mcmc_sample[k,] = proposal
      if (k > n_burn) {
        accept <- accept + 1 / (K - n_burn)
      }
      accept_adapt <- accept_adapt + 1/ 50
    }else{
      mcmc_sample[k,] = current_mcmc_sample
    }
    # update step
    if (update_tune & (k < n_burn)){ ## check whether to update the tuning
      if(k %in% bathch_idx){
        if(k %in% bathch_idx){
          batch_samples <- mcmc_sample[(k-50) : k, ]
          #accept <- sum(accept[(i-50) : i])/nrow(batch_samples)
          param_update <- update_tuning_mv(k, accept_adapt, lambda, batch_samples,
                                           Sigma_tune, Sigma_tune_chol)
          
          Sigma_tune_chol <- param_update$Sigma_tune_chol
          lambda          <- param_update$lambda
          Sigma_tune      <- param_update$Sigma_tune
          accept_adapt    <- param_update$accept
        }
      }
    }
  }
  message("Acceptance rate = ", accept)
  # if (update_tune) {
  #   message("Final tuning parameter is ", Sigma_tune_chol)
  # }
  return(mcmc(mcmc_sample))
}

# Initial values

Sig_chol <- matrix(0, p+1, p+1)
Sig_chol[lower.tri(Sig_chol, diag = TRUE)] <- 1

# non_adapt_out10 <- Adapt_MH_MCMC_lm(X = X, y = y, K = 10000, Sigma_tune = diag(10, p+1), Sig_chol, update_tune = FALSE)
# non_adapt_out1  <- Adapt_MH_MCMC_lm(X = X, y = y, K = 10000, Sigma_tune = diag(1,  p+1), Sig_chol, update_tune = FALSE)
# non_adapt_out0.1 <- Adapt_MH_MCMC_lm(X = X, y = y, K = 10000, Sigma_tune = diag(0.1, p+1), Sig_chol, update_tune = FALSE)
#adap_out <- Adapt_MH_MCMC_lm(X = X, y = y, K = 10000, Sigma_tune = diag(1, p+1), Sig_chol, update_tune = TRUE)

# save(non_adapt_out10, file = here::here("results", "non_adapt_out10.RData"))
# save(non_adapt_out1 , file = here::here("results", "non_adapt_out1.RData"))
# save(non_adapt_out0.1, file = here::here("results", "non_adapt_out0.1.RData"))
# save(adap_out, file = here::here("results", "adap_out.RData"))
# 
# png(file = here::here("figures", "compare-adap-non-adapt.png"),  #Creating a plotting window in folder
#     res = 400, height = 6, width = 6, units = "in") # This number can be adjusted for a better view
# 
# par(mfrow =c(2, 2))
# 
# colors <- c("blue", "red", "green", "purple", "darkblue", "black", "violet")
# matplot(non_adapt_out10, type = "l", col = colors, main = "prop sd = 10", ylim = c(-8, 10))
# matplot(non_adapt_out1, type = "l", col = colors, main = "prop sd = 1", ylim = c(-8, 10))
# matplot(non_adapt_out0.1, type = "l", col = colors, main = "prop sd  = 0.1", ylim = c(-8, 10))
#matplot(adap_out, type = "l",  main = "adaptive mcmc", ylim = c(-8, 10))
# 
# dev.off() #Save the graph to the window
# 
# png(file = here::here("figures", "adap-non-adapt-post-mean.png"),  #Creating a plotting window in folder
#     res = 400, height = 6, width = 6, units = "in") # This number can be adjusted for a better view
# 
# par(mfrow =c(2, 2))
# plot(apply(non_adapt_out10, 2, mean), c(truebeta, truesigma2), main = "Estimates vs true parameters sd = 10",
#      cex.main = 0.7, xlab = "parameter posterior mean", ylab = "true parameter values", cex.lab = 0.7)
# abline(0, 1, col = "red")
# 
# plot(apply(non_adapt_out1, 2, mean), c(truebeta, truesigma2), main = "Estimates vs true parameters sd = 1",
#      cex.main = 0.7, xlab = "parameter posterior mean", ylab = "true parameter values",  cex.lab = 0.7)
# abline(0, 1, col = "red")
# 
# plot(apply(non_adapt_out0.1, 2, mean), c(truebeta, truesigma2), main = "Estimates vs true parameters sd = 0.1",
#      cex.main = 0.7, xlab = "parameter posterior mean", ylab = "true parameter values", cex.lab = 0.7)
# abline(0, 1, col = "red")
# plot(apply(adap_out, 2, mean), c(truebeta, truesigma2), main = "adaptive mcmc estimates vs true parameters",
#      cex.main = 0.7, xlab = "parameter posterior mean", ylab = "true parameter values", cex.lab = 0.7)
# abline(0, 1, col = "red")
# 
# dev.off() 
# 
# # Effective sample sizes
# 
# knitr::kable(cbind(
#   effectiveSize(non_adapt_out10), 
#   effectiveSize(non_adapt_out1), 
#   effectiveSize(non_adapt_out0.1),
#   effectiveSize(adap_out)), 
#   col.names = c('ESS with sd = 10', 'ESS with sd = 1', 'ESS with sd = 0.1', 'adaptive ESS'), "rst")
# 




