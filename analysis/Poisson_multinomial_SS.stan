// Poisson-multinomial model for Poisson state-space process

data {
  int<lower=0> T;      // length of time series
  int<lower=0> y[T];   // time series of counts
}

transformed data {
  int<lower=0> N = sum(y);  // sum of all counts  
}

parameters {
  real<lower=0> sigma; // process error SD
  vector[T] z;         // process error innovations (Z-scored)
}

transformed parameters {
  vector[T] x;         // log-means (states)
  
  // Random walk process model
  // use hierarchical non-centered parameterization aka "Matt trick"
  x[1] = sigma*z[1];
  for(t in 2:T)
    x[t] = x[t-1] + sigma*z[t];
}

model {
  vector[T] lambda = exp(x);         // means
  real sum_lambda = sum(lambda);     // sum of means
  vector[T] pi = lambda/sum_lambda;  // cell probabilities
  
  // Priors
  sigma ~ normal(0,3);
  
  // Process model
  z ~ std_normal();         // implies x[t] ~ N(x[t-1], sigma)
  
  // Observation model
  N ~ poisson(sum_lambda);  // marginal distribution of total count
  y ~ multinomial(pi);      // conditional distribution of counts
}

