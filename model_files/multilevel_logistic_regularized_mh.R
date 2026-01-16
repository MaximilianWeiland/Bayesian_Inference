log_posterior <- function(beta, sigma_beta, gamma, X, U, y, country, J, K, L) {
  
  # hyperprior on sigma-beta
  log_prior_sigma_beta <- sum(dnorm(sigma_beta, mean = 0, sd = .5, log = TRUE))
  
  # regularized double-exponential prior on gamma
  lambda <- 1 / 0.1
  log_prior_gamma <- sum(log(lambda) - log(2) - lambda * abs(gamma))
  
  # prior on beta-values
  log_prior_beta <- 0
  for (j in 1:J) {
    mu_j <- as.numeric(U[j, ] %*% gamma)
    for (k in 1:K) {
      log_prior_beta <- log_prior_beta +
        dnorm(beta[j, k], mu_j[k], sigma_beta[k], log = TRUE)
    }
  }
  
  # likelihood
  log_lik <- 0
  for (i in 1:length(y)) {
    eta <- sum(X[i, ] * beta[country[i], ])
    p <- 1/(1+exp(-eta))
    log_lik <- log_lik + dbinom(y[i], size = 1, prob = p, log = TRUE)
  }
  
  # combine all to log posterior
  log_post <- log_lik + log_prior_gamma + log_prior_sigma_beta + log_prior_beta
  log_post
}

mh <- function(X, U, y, country,
               n_samples = 2000, warmup = 500, proposal_sd = 0.01,
               beta_start = NULL, gamma_start = NULL, sigma_beta_start = NULL) {
  
  N <- nrow(X)
  K <- ncol(X)
  J <- max(country)
  L <- ncol(U)
  
  # initialize the starting parameters
  if (is.null(beta_start)) beta_start <- matrix(0, J, K)
  if (is.null(gamma_start)) gamma_start <- matrix(0, L, K)
  if (is.null(sigma_beta_start)) sigma_beta_start <- rep(1, K)
  
  # initialize objects to store the sampling results
  beta_samples <- array(NA, dim = c(n_samples, J, K))
  gamma_samples <- array(NA, dim = c(n_samples, L, K))
  sigma_beta_samples <- matrix(NA, nrow = n_samples, ncol = K)
  ll_samples <- numeric(n_samples)
  accept <- logical(n_samples)
  
  # initialize the starting value and calculate first ll
  beta_samples[1, , ] <- beta_start
  gamma_samples[1, , ] <- gamma_start
  sigma_beta_samples[1, ] <- sigma_beta_start
  ll_samples[1] <- log_posterior(beta_start, sigma_beta_start, gamma_start,
                                 X, U, y, country, J, K, L)
  accept[1] <- TRUE
  
  # after start, loop over the number of samples
  for (t in 2:n_samples) {
    
    # extract current parameter values
    beta_curr <- beta_samples[t - 1, , ]
    gamma_curr <- gamma_samples[t - 1, , ]
    sigma_beta_curr <- sigma_beta_samples[t - 1, ]
    
    # flatten the matrices to vectors
    beta_vec_curr <- as.vector(beta_curr)
    gamma_vec_curr <- as.vector(gamma_curr)
    
    # get proposal values as vectors
    beta_vec_prop <- mvrnorm(1, mu = beta_vec_curr,
                             Sigma = diag(proposal_sd^2, length(beta_vec_curr)))
    gamma_vec_prop <- mvrnorm(1, mu = gamma_vec_curr,
                              Sigma = diag(proposal_sd^2, length(gamma_vec_curr)))
    
    # reshape back to matrices:
    beta_prop <- matrix(beta_vec_prop, nrow = J, ncol = K)
    gamma_prop <- matrix(gamma_vec_prop, nrow = L, ncol = K)
    
    # get proposal for sigma beta (propose on log to keep sd positive)
    log_sigma_beta_curr <- log(sigma_beta_curr)
    log_sigma_beta_prop <- rnorm(length(log_sigma_beta_curr), mean = log_sigma_beta_curr, sd = proposal_sd)
    sigma_beta_prop <- exp(log_sigma_beta_prop)
    
    # calculate the acceptance ratio
    log_r <- log_posterior(beta_prop, sigma_beta_prop, gamma_prop,
                           X, U, y, country, J, K, L) -
      log_posterior(beta_curr, sigma_beta_curr, gamma_curr,
                    X, U, y, country, J, K, L)
    
    if (log(runif(1)) < log_r) { # acceptance
      beta_samples[t, , ] <- beta_prop
      gamma_samples[t, , ] <- gamma_prop
      sigma_beta_samples[t, ] <- sigma_beta_prop
      accept[t] <- TRUE
      ll_samples[t] <- log_posterior(beta_prop, sigma_beta_prop, gamma_prop,
                                     X, U, y, country, J, K, L)
    } else {                      # rejection
      beta_samples[t, , ] <- beta_curr
      gamma_samples[t, , ] <- gamma_curr
      sigma_beta_samples[t, ] <- sigma_beta_curr
      accept[t] <- FALSE
      ll_samples[t] <- log_posterior(beta_curr, sigma_beta_curr, gamma_curr,
                                     X, U, y, country, J, K, L)
    }
  }
  
  # remove the first 500 iterations due to warmup
  beta_samples <- beta_samples[(warmup + 1):n_samples, , ]
  gamma_samples <- gamma_samples[(warmup + 1):n_samples, , ]
  sigma_beta_samples <- sigma_beta_samples[(warmup + 1):n_samples, ]
  ll_samples <- ll_samples[(warmup + 1):n_samples]
  accept <- accept[(warmup + 1):n_samples]
  
  # generate replicate data for the ppc
  y_rep_samples <- matrix(0, nrow = dim(beta_samples)[1], ncol = nrow(X))
  for (i in 1:nrow(beta_samples)) {
    eta <- numeric(nrow(X))
    for (n in 1:nrow(X)) {
      eta[n] <- sum(X[n, ] * beta_samples[i, country[n], ])
    }
    p <- 1 / (1 + exp(-eta))
    y_rep_samples[i, ] <- rbinom(nrow(X), size = 1, prob = p)
  }
  
  # collect all results inside a list
  list(beta = beta_samples, gamma = gamma_samples,
       sigma_beta = sigma_beta_samples, ll = ll_samples,
       y_rep = y_rep_samples, accepted = accept)
}
