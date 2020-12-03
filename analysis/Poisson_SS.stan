// Poisson state-space model

data {
  int<lower=0> T;      // length of time series
  int<lower=0> y[T];   // time series of counts
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
  // Priors
  sigma ~ normal(0,3);
  
  // Process model
  z ~ std_normal();   // implies x[t] ~ N(x[t-1], sigma)
  
  // Observation model
  y ~ poisson_log(x);
}

