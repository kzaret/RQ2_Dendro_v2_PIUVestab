// Poisson-multinomial model for Poisson state-space process, 
// with observation error in times

functions {
  // Truncated geometric PMF, evaluated over domain x 
  vector geometric_pmf(vector x, real r) {
    vector[rows(x)] px;
    px = r * exp(x * log1m(r));
    return(px/sum(px));
  }
}

data {
  int<lower=0> T;      // length of time series
  int<lower=0> y[T];   // time series of observed counts
  real<lower=0> r;     // probability parameter for geometric obs error in time (known)
}

transformed data {
  int<lower=0> N = sum(y);  // sum of all counts  
  vector<lower=1>[T] t1T;   // 1:T
  matrix<lower=0>[T,T] gamma = rep_matrix(0,T,T); // conditional likelihood P(t | tau, r) 
  
  for(t in 1:T) 
    t1T[t] = t;
  
  // Geometric observation error likelihood of observed time t (rows) 
  // given discrete latent state tau = 1, ..., T (cols)
  for(tau in 1:T)
    gamma[tau:T,tau] = geometric_pmf(t1T[tau:T] - tau, r);
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
  vector[T] lambda = exp(x);          // means
  real sum_lambda = sum(lambda);      // sum of means
  vector[T] pi = lambda/sum_lambda;   // cell probabilities

  // Priors
  sigma ~ normal(0,3);
  
  // Process model
  z ~ std_normal();         // implies x[t] ~ N(x[t-1], sigma)
  
  // Observation model
  // Marginal likelihood of t (t = 1, ..., T), summing over discrete latent states tau:
  // sum(prior * likelihood), where the geometric observation error likelihood is
  // weighted by the process-model prior
  y ~ multinomial(gamma * pi);   // distribution of counts by observed year
  N ~ poisson(sum_lambda);       // marginal distribution of total count
}
