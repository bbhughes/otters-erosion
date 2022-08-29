data {
  int<lower=0> N; // number of data points
  real width[N];  // observed mean width
  real otters[N]; // observed number of otters, with near 0 error
}
parameters {
  real<lower=0> r; // base rate of width erosion without otters
  real<lower=0> a; // initial width
  real<lower=0> b; // reduction in erosion per otter at max effect (i.e. near otters = 0)
  real<lower=0> sigma; // observation error SD on width
}
transformed parameters {
  vector[N] eta; // predicted rate of erosion
  vector[N] zeta; // predicted rate of erosion without otters
  vector[N] reff; // effective rate of erosion
  eta[1] = a; // start the time series
  zeta[1] = a; // start the time series
  reff[1] = r; // start the time series
  for (i in 2:N) {
    reff[i] = r -  // base rate of erosion
              r * (1 - exp(-b * otters[i])); // minus otter effect
    eta[i] = eta[i-1] * exp(reff[i]); // last year width x erosion rate
    zeta[i] = zeta[i-1] * exp(r);
  }
}
model {
  target += student_t_lpdf(sigma | 3, 0, 1); // prior
  target += normal_lpdf(b | 0, 1); // prior
  target += normal_lpdf(a | 0, 3); // prior
  target += normal_lpdf(r | 0, 1); // prior
  target += lognormal_lpdf(width | log(eta), sigma); // data likelihood
}
