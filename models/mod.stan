data {
  int<lower=0> n_levels;
  int<lower=0> n;
  int variable[n];
  int era[n];
  vector[n] y;
  vector[n] x;
}
parameters {
  vector[n_levels] beta;
  vector[n_levels] ratio;
  real<lower=0> sigma_resid;
  real<lower=0> sigma_ratio;
  real<lower=0> sigma_beta;
  real mu_ratio;
  real mu_beta;
}
transformed parameters {
  vector[n_levels] exp_ratio;
  real exp_mu_ratio;
  vector[n] pred;
  exp_mu_ratio = exp(mu_ratio);
  for(i in 1:n_levels) {exp_ratio[i] = exp(ratio[i]);}
  for(i in 1:n) {
    if(era[i]==1) {
      pred[i] = x[i]*beta[variable[i]];
    } else {
      pred[i] = x[i]*(beta[variable[i]] * ratio[variable[i]]);
    }
  }
}
model {
  mu_beta ~ normal(0,1);
  sigma_beta ~ student_t(3,0,2);
  mu_ratio ~ normal(0,1);
  sigma_ratio ~ student_t(3,0,2);
  sigma_resid ~ student_t(3,0,2);
  beta ~ normal(mu_beta,sigma_beta);
  ratio ~ normal(mu_ratio,sigma_ratio);
  y ~ normal(pred, sigma_resid);
}
