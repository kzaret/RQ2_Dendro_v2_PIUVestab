// Multinomial model for distribution of times in a Poisson state-space model,
// conditioned on total count

data {
  int<lower=0> T;      // length of time series
  int<lower=0> y[T];   // time series of counts
}

transformed data {
  int<lower=0> N = sum(y);  // sum of all counts  
}

parameters {
  real<lower=0> sigma; // process errror SD
  vector[T] z;         // process error innovations (Z-scored)
}

transformed parameters {
  vector[T] psi;       // additive log-ratios of yearly probabilities
  vector[T] pi;        // yearly probabilities (states)
  
  // Random walk process model
  // use hierarchical non-centered parameterization aka "Matt trick"
  psi[1] = 0;          // year 1 is reference class
  for(t in 2:T)
    psi[t] = psi[t-1] + sigma*z[t];
  pi = softmax(psi);   // inverse log-ratio transform
}

model {
  // Priors
  sigma ~ normal(0,3);
  
  // Process model
  z ~ std_normal();    // implies psi[t] ~ N(psi[t-1], sigma)
  
  // Observation model
  y ~ multinomial(pi); // distribution of counts by year
}

generated quantities {
  vector[T] lambda = pi*N;  // predicted counts by year 
}
