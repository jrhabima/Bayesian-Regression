#' this function updates the univariate Gaussian random walk proposal Metropolis-Hastings tuning parameters
#' @param k is a positive integer that is the current MCMC iteration.
#' @param accept is a number between 0 and 1 that represents the batch accpetance rate. The target of this adaptive tuning algorithm is an acceptance rate of 0.44.
#' @param tune is a positive number that is the univariate Gaussian random walk proposal standard deviation.
#' 
#' @keywords internal

update_tuning <- function(k, accept, tune) {
    delta = 1.0 / sqrt(k)
    tune_out <- 0
    if(accept > 0.44){
        tune_out <- exp(log(tune) + delta)
    } else {
        tune_out <- exp(log(tune) - delta)
    }
    accept_out <- 0.0
    return(
        list(
            accept = accept_out,
            tune   = tune_out
        )
    )
}


#' this function updates a block Gaussian random walk proposal Metropolis-Hastings tuning parameters
#' @param k is a positive integer that is the current MCMC iteration.
#' @param accept is a positive number between 0 and 1 that represents the batch accpetance rate. The target of this adaptive tuning algorithm is an acceptance rate of between 0.44 (a univariate update) and 0.234 (5+ dimensional upate).
#' @param lambda is a positive scalar that scales the covariance matrix Sigma_tune and diminishes towards 0 as k increases.
#' @param batch_samples is a \eqn{50 \times d}{50 x d} dimensional matrix that consists of the 50 batch samples of the d-dimensional parameter being sampled.
#' @param Sigma_tune is a \eqn{d \times d}{d x d} positive definite covariance matrix of the batch samples used to generate the multivariate Gaussian random walk proposal.
#' @param Sigma_tune_chol is the \eqn{d \times d}{d x d}  Cholesky decomposition of Sigma_tune.
#' 
#' @keywords internal

update_tuning_mv <- function(k, accept, lambda, batch_samples,
                             Sigma_tune, Sigma_tune_chol) {
    arr <- c(0.44, 0.35, 0.32, 0.25, 0.234)
    # std::vector<double> acceptance_rates (arr, arr + sizeof(arr) / sizeof(arr[0]))
    dimension <- nrow(batch_samples)
    if (dimension >= 5) {
        dimension <- 5
    }
    d <- ncol(batch_samples)
    batch_size <- nrow(batch_samples)
    optimal_accept <- arr[dimension]
    times_adapted <- floor(k / 50)
    gamma1 <- 1.0 / ((times_adapted + 3.0)^0.8)
    gamma2 <- 10.0 * gamma1
    adapt_factor <- exp(gamma2 * (accept - optimal_accept))
    ## update the MV scaling parameter
    lambda_out <- lambda * adapt_factor
    ## center the batch of MCMC samples
    batch_samples_tmp <- batch_samples
    for (j in 1:d) {
        mean_batch = mean(batch_samples[, j])
        for (i in 1:batch_size) {
            batch_samples_tmp[i, j] <- batch_samples[i, j] - mean_batch
        }
    }
    ## 50 is an MCMC batch size, can make this function more general later...
    Sigma_tune_out <- Sigma_tune + gamma1 *
        (t(batch_samples) %*% batch_samples / (50.0-1.0) - Sigma_tune)
    Sigma_tune_chol_out <- tryCatch(
        chol(Sigma_tune),
        error = function(e) {
            chol(Sigma_tune + 1e-8 * diag(nrow(Sigma_tune)))                    
        }
    )
    accept_out <- 0.0
    batch_samples_out <- matrix(0, batch_size, d)
    return(
        list(
            accept          = accept_out,
            lambda          = lambda_out,
            batch_samples   = batch_samples_out,
            Sigma_tune      = Sigma_tune_out,
            Sigma_tune_chol = Sigma_tune_chol_out
        )
    )
}
