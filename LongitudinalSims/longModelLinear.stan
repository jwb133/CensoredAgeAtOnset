data {
  int<lower=1> n_subj;
  int<lower=1> n_parents;
  int<lower=1> n_obs;
  int<lower=1> n_mutation;

  vector[n_parents] parentOnset;

  array[n_subj] int<lower=1, upper=n_mutation> mutation_by_sub;
  array[n_parents] int<lower=1, upper=n_mutation> mutation_by_parent;
  array[n_obs] int<lower=1, upper=n_subj> id;
  array[n_obs] int<lower=1, upper=n_mutation> mutation;

  vector[n_obs] age;
  vector[n_obs] y;

  int<lower=0> n_obs_onset;
  int<lower=0> n_cens_onset;

  array[n_obs_onset] int<lower=1, upper=n_subj> subj_obs_onset;
  array[n_cens_onset] int<lower=1, upper=n_subj> subj_cens_onset;

  vector[n_obs_onset] onset_obs;
  vector[n_cens_onset] censor_age;
}

parameters {
  // age-at-onset model
  real ageAtOnsetInt;
  real<lower=1e-6> sigma_ageAtOnset;

  // latent censored onset times
  vector<lower=censor_age>[n_cens_onset] onset_cens;

  // subject random effects
  vector<lower=1e-6>[2] sd_alpha;
  cholesky_factor_corr[2] Lcorr;
  matrix[2, n_subj] z_alpha;

  // mutation effects
  real<lower=1e-6> sd_mutation;
  real<lower=1e-6> sd_aoo_mutation;
  vector[n_mutation] z_mut;
  vector[n_mutation] z_aoo_mut;

  // fixed effects
  real alpha0;
  real b;

  // residual SD
  real<lower=1e-6> sdYRes;
}

transformed parameters {
  vector[n_mutation] m;
  vector[n_mutation] m_aoo;
  matrix[2, n_subj] alpha;
  vector[n_subj] ageAtOnset;

  m = sd_mutation * z_mut;
  m_aoo = sd_aoo_mutation * z_aoo_mut;
  alpha = diag_pre_multiply(sd_alpha, Lcorr) * z_alpha;

  // fill in subject-specific onset ages
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
  sigma_ageAtOnset ~ normal(0, 10);

  sd_alpha[1] ~ normal(0, 5);
  sd_alpha[2] ~ normal(0, 2);
  Lcorr ~ lkj_corr_cholesky(1);
  to_vector(z_alpha) ~ normal(0, 1);

  sd_mutation ~ normal(0, 5);
  sd_aoo_mutation ~ normal(0, 5);
  z_mut ~ normal(0, 1);
  z_aoo_mut ~ normal(0, 1);

  alpha0 ~ normal(0, 100);
  b ~ normal(0, 100);
  sdYRes ~ normal(0, 5);

  // observed onsets
  for (j in 1:n_obs_onset) {
    int i = subj_obs_onset[j];
    onset_obs[j] ~ normal(ageAtOnsetInt + m_aoo[mutation_by_sub[i]],
                          sigma_ageAtOnset);
  }

  // right-censored onsets
  for (j in 1:n_cens_onset) {
    int i = subj_cens_onset[j];
    onset_cens[j] ~ normal(ageAtOnsetInt + m_aoo[mutation_by_sub[i]],
                           sigma_ageAtOnset);
  }

  // parent onsets
  for (i in 1:n_parents) {
    parentOnset[i] ~ normal(ageAtOnsetInt + m_aoo[mutation_by_parent[i]],
                            sigma_ageAtOnset);
  }

  // longitudinal model
  for (i in 1:n_obs) {
    real t_i = age[i] - ageAtOnset[id[i]];
    real mu_i = alpha0
                + m[mutation[i]]
                + b * t_i
                + alpha[1, id[i]]
                + alpha[2, id[i]] * t_i;

    y[i] ~ normal(mu_i, sdYRes);
  }
}

generated quantities {
  corr_matrix[2] Cor_alpha;
  Cor_alpha = multiply_lower_tri_self_transpose(Lcorr);
}