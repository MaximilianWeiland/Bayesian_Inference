# define function to calculate the log posterior 
log_posterior <- function(beta, X, y) {
  eta <- X %*% beta
  p <- 1 / (1 + exp(-eta))
  loglik <- sum(dbinom(y, size = 1, prob = p, log = TRUE))
  logprior <- sum(dnorm(beta, mean = 0, sd = 10, log = TRUE))
  log_posterior <- loglik + logprior
  log_posterior
}

# define the Metropolis-Hastings sampler
mh <- function(X, y, n_samples = 1500, warmup = 500,
               proposal_sd_vec = rep(0.1, ncol(X)),
               start_beta = rep(1, ncol(X))) {
  
  # initialize objects to store the sampling results
  beta_samples <- matrix(NA, nrow = n_samples, ncol = ncol(X))
  ll_samples <- numeric(n_samples)
  accept <- logical(n_samples)
  
  # initialize the starting value and calculate first ll
  beta_samples[1, ] <- start_beta
  ll_samples[1] <- log_posterior(start_beta, X, y)
  accept[1] <- TRUE
  
  # after start, loop over the number of samples
  for (t in 2:n_samples) {
    
    # look at current beta
    beta_curr <- beta_samples[t - 1, ]
    
    # propose new beta values
    Sigma_prop <- diag(proposal_sd_vec^2, ncol(X))
    beta_prop <- mvrnorm(1, mu = beta_curr, Sigma = Sigma_prop)

    # calculate the acceptance ratio
    log_r <- log_posterior(beta_prop, X, y) - log_posterior(beta_curr, X, y)
    
    if (log(runif(1)) < log_r) { # acceptance
      beta_samples[t, ] <- beta_prop
      accept[t] <- TRUE
      ll_samples[t] <- log_posterior(beta_prop, X, y)
    } else {                      # rejection
      beta_samples[t, ] <- beta_curr
      accept[t] <- FALSE
      ll_samples[t] <- log_posterior(beta_curr, X, y)
    }
  }
  
  # remove the first 500 iterations due to warmup
  beta_samples <- beta_samples[(warmup + 1):n_samples, , drop = FALSE]
  ll_samples <- ll_samples[(warmup + 1):n_samples]
  accept <- accept[(warmup + 1):n_samples]
  
  # generate replicate data for the ppc
  y_rep_samples <- matrix(0, nrow = nrow(beta_samples), ncol = nrow(X))
  for (i in 1:nrow(beta_samples)) {
    p <- 1/(1+exp(-(X %*% beta_samples[i, ])))
    y_rep_samples[i, ] <- rbinom(nrow(X), size = 1, prob = p)
  }
  
  list(beta = beta_samples, ll = ll_samples,
       y_rep = y_rep_samples, accepted = accept)
}
