// data
data {
  int<lower=1> N;  // number of observations
  int<lower=1> J;  // number of groups
  int<lower=1> K;  // number of individual-level predictors + 1 
  int<lower=1> L;  // number of group-level predictors + 1 
  array[N] int<lower=1,upper=J> country;  // country for each observation
  array[J] row_vector[L] U;  // group-level predictor, J*L matrix
  array[N] row_vector[K] X;  // individual-level predictor, N*K matrix
  array[N] int<lower=0, upper=1> y;  // binary outcome y
}

// parameters
parameters {
  matrix[J,K]  beta;  // matrix of coefficients of individual-level regression J * K  
  matrix[L,K] gamma;  // matrix of coefficients of group-level regression L*K
  cholesky_factor_corr[K] LKJ;  // cholesky factor of a correlation matrix
  vector<lower=0>[K] D;  // scale vector
}

// model
model {
  // priors
  // create covariance matrix with LKJ distribution and cholesky decomposition
  LKJ ~ lkj_corr_cholesky(1);  // prior for the correlation 
  D ~ normal(0,1);  // prior on paramater scale (sd)
  matrix[K,K] Sigma;  // matrix to store covariance matrix of beta
  Sigma = diag_pre_multiply(D, LKJ);  // prior for Sigma  D*L
  
  // weakly informative normal prior on all gamma values
  to_vector(gamma) ~ normal(0, 10);

  matrix[J, K] mu;  // matrix to store prior means
  // loop over all groups and sample beta values
  for (j in 1:J){
     mu[j] = U[j,] * gamma;  // calculate the group level mean vector
     beta[j] ~ multi_normal_cholesky(mu[j], Sigma);  // prior on individual-level coefficients
  }

  // likelihood
  vector[N] y_hat;  // vector to store all predictions
  // loop over all observations
  for (i in 1:N){
    y_hat[i] = dot_product(X[i,], beta[country[i],]);  // get prediction
    y[i] ~ bernoulli_logit(y_hat[i]);  // calculate likelihood of prediction
  }  
}
generated quantities {
  array[N] int<lower=0, upper=1> y_rep;
  for (i in 1:N) {
    y_rep[i] = bernoulli_logit_rng(dot_product(X[i,], beta[country[i],]));
  }
}
