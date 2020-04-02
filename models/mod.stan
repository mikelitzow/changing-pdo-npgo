// Stan Code for Version 6 of Salmon SR Era Regression Model

//Notes:
  //  matrix[3, 3] m[6, 7] - m to be a two-dimensional array of size 6 x 7,
//    containing values that are 3 x 3 matrices.

data {

  int<lower=0> S; //number of stocks
  int N[S];  //number of SR observations for each stock
  int<lower=0> maxN;  //maximum number of observations across stocks
  //int years[S, maxN]; //Pointer vector for brood year references.
  matrix[S,maxN] ln_rps;
  matrix[S,maxN] spawn;
  int era[S,maxN]; //Indicator for before/after 89

  // Covariates
  matrix[S,maxN] covar;

  int<lower=0> R; //number of regions
  int region[S]; //region pointer

}

parameters {
  // Coefficients: Stock Level
  vector[S] beta;

  // Coef. Ratios: Region Level
  vector[R] ratio;

  //Ricker Params
  real ricker_mu_alpha;
  real<lower=0> ricker_sigma_alpha;
  real ricker_alpha[S]; 
  real<lower=0> ricker_beta[S]; 
  //Variances
  real<lower=0> sigma_resid[S];
  //Autoregressive term for errors
  vector[S] phi;
}

transformed parameters {
  // Predited Log-R/S
  vector[maxN] pred[S]; // [S,maxN] - vector of length S, each element of which is a vector of length maxN

  //Iniital Residual for autocorrelated error structure
  vector[maxN] residual[S];

  // Convert Ratios
  // vector[S] exp_ratio;

  // Convert phi to be (-1, 1)
  vector[S] phi_trans;
  
  //Beta in the second era
  vector[S] beta2; 
  
  //Generate Predicted values
  for(s in 1:S) {
    beta2[s] = beta[s] * ratio[region[s]];
    
    phi_trans[s] = 2*exp(phi[s])/(1+exp(phi[s])) - 1;
    //Alternatively, we may be able to do this with an if statement.... But it makes the JAGSer inside of me cry :(
      for(n in 1:maxN) {
        if(n<=N[s]) {
          // Calculate Residual for AR(1) Error
          if(n==1) {
            residual[s,n] = 0;
          }else {
            residual[s,n] = ln_rps[s,n-1] - pred[s,n-1];
          }

          // Calculate predicted ln(R/S)
          if(era[s,n]==1) {
            pred[s,n] = ricker_alpha[s] - ricker_beta[s]*spawn[s,n] + covar[s,n]*beta[s];//  +
              // phi_trans[s] * residual[s,n];
          }else {
            pred[s,n] = ricker_alpha[s] - ricker_beta[s]*spawn[s,n] + covar[s,n]*beta2[s];//  +
              // phi_trans[s] * residual[s,n];
          }
        }else {
          pred[s,n] = 0;
          residual[s,n] = 0;
        }
      }// next n
  }//next s
}

model {
  //Priors
  phi ~ normal(0,2); //Autoregressive term for errors

  // Ricker alpha hyperparameters
  ricker_mu_alpha ~ normal(0,5);
  ricker_sigma_alpha ~ normal(0,5);

  for(s in 1:S) { // Stocks
    ricker_alpha[s] ~ normal(ricker_mu_alpha,ricker_sigma_alpha);
    ricker_beta[s] ~ normal(0,0.001);
    sigma_resid[s] ~ normal(0,1);
    
    //Coefficients
    beta[s] ~ normal(0,1);
  }//next s
  



  for(r in 1:R) { //Regions
    ratio[r] ~ normal(1,1);
  }


  //Likelihood
  for(s in 1:S) {
    for(n in 1:N[s]) {
      ln_rps[s,n] ~ normal(pred[s,n] + phi_trans[s] * residual[s,n], sigma_resid[s]);
    }//next n
  }//next s
}

generated quantities {
  vector[maxN] log_lik[S];
  vector[S] beta_ratio; //Beta2/Beta1
  
  for(s in 1:S) {
    // for(n in 1:N[s]){
      for(n in 1:maxN) {
        if(n<=N[s]) {
          log_lik[s,n] = normal_lpdf(ln_rps[s,n] | pred[s,n] + phi_trans[s] * residual[s,n], sigma_resid[s]);
        }else {
          log_lik[s,n] = 0;
        }
      }//next n
      
      //Ratio of Betas
      beta_ratio[s] = beta2[s]/beta[s];
    }//next s

  }


