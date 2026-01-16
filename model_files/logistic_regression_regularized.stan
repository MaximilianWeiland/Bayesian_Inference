data {
  int<lower=0> N;                 // number of observations
  int<lower=0> K;                 // number of predictors
  matrix[N, K] X;                 // predictor matrix
  array[N] int<lower=0,upper=1> y;   
}
parameters {
  vector[K] beta;                 // coefficients
}
model {
  // Priors
  beta ~ double_exponential(0, 0.1);
  // Likelihood
  y ~ bernoulli_logit(X * beta); // For logit
}
generated quantities {
  array[N] int<lower=0, upper=1> y_rep;
  for (i in 1:N) {
    y_rep[i] = bernoulli_logit_rng(dot_product(X[i], beta));
  }
}
