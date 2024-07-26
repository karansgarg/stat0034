// A standard implementation of a 10D Neal's funnel in Stan

parameters {
  real y;
  vector[9] x;
}
model {
  y ~ normal(0, 3);
  x ~ normal(0, exp(y/2));
}
