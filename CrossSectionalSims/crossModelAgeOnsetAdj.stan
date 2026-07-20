data {
  int<lower=1> n_subj;

  // subject-level parental onset covariate
  vector[n_subj] parentOnset;

  // biomarker measurement age and value, one per subject
  vector[n_subj] age;
  vector[n_subj] y;

  // number with onset observed
  int<lower=0> n_obs_onset;
  // number with onset censored
  int<lower=0> n_cens_onset;

  // ids of those with onset observed
  array[n_obs_onset] int<lower=1, upper=n_subj> subj_obs_onset;
  // ids of those with onset right censored
  array[n_cens_onset] int<lower=1, upper=n_subj> subj_cens_onset;

  // observed onset values
  vector[n_obs_onset] onset_obs;
  // right-censored onset values (i.e. age at censoring)
  vector[n_cens_onset] censor_age;
}

parameters {
  // age-at-onset model
  real ageAtOnsetInt;
  real ageAtOnsetParent;
  real<lower=1e-6> sigma_ageAtOnset;

  // latent censored onset times
  vector<lower=censor_age>[n_cens_onset] onset_cens;

  // biomarker model
  real alpha_0;
  real alpha_t;
  real alpha_o;
  real<lower=1e-6> sdYRes;
}

transformed parameters {
  vector[n_subj] ageAtOnset;

  for (j in 1:n_obs_onset) {
    ageAtOnset[subj_obs_onset[j]] = onset_obs[j];
  }

  for (j in 1:n_cens_onset) {
    ageAtOnset[subj_cens_onset[j]] = onset_cens[j];
  }
}

model {
  // priors
  ageAtOnsetInt ~ normal(0, 100);
  ageAtOnsetParent ~ normal(0, 100);
  sigma_ageAtOnset ~ normal(0, 10);

  alpha_0 ~ normal(0, 100);
  alpha_t ~ normal(0, 100);
  alpha_o ~ normal(0, 100);
  sdYRes ~ normal(0, 5);

  // observed onset times
  for (j in 1:n_obs_onset) {
    int i = subj_obs_onset[j];
    onset_obs[j] ~ normal(ageAtOnsetInt
                          + ageAtOnsetParent * parentOnset[i],
                          sigma_ageAtOnset);
  }

  // censored onset times
  for (j in 1:n_cens_onset) {
    int i = subj_cens_onset[j];
    onset_cens[j] ~ normal(ageAtOnsetInt
                           + ageAtOnsetParent * parentOnset[i],
                           sigma_ageAtOnset);
  }

  // biomarker model
  for (i in 1:n_subj) {
    real t_i = age[i] - ageAtOnset[i];
    real mu_i = alpha_0 + alpha_t * t_i + alpha_o * ageAtOnset[i];
    y[i] ~ normal(mu_i, sdYRes);
  }
}

generated quantities {
  real ageAtOnsetIntOrig;
  ageAtOnsetIntOrig = ageAtOnsetInt+40*(1-ageAtOnsetParent);
  
  real alpha_0Orig, alpha_tOrig, alpha_oOrig, sdYResOrig;
  alpha_tOrig = alpha_t*10;
  alpha_oOrig = alpha_o*10;

  alpha_0Orig = 100+10*(alpha_0-40*alpha_o);
  sdYResOrig = sdYRes*10;
  
  real alpha_t_minus_o;
  alpha_t_minus_o = alpha_tOrig-alpha_oOrig;
  
  real fitted_mean;
  fitted_mean = alpha_0Orig -6*alpha_tOrig + 46*alpha_oOrig;
}
