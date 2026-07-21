functions {
  // Restricted cubic spline basis (3 terms) with knots
  vector spline_basis(real t, vector knot) {
    real term1 = fmax((t - knot[1])^3, 0);
    real term2 = fmax((t - knot[2])^3, 0);
    real term3 = fmax((t - knot[3])^3, 0);
    real term4 = fmax((t - knot[4])^3, 0);
    real u1 = t;
    real u2 = (term1 - (term3*(knot[4]-knot[1]) - term4*(knot[3]-knot[1]))/(knot[4]-knot[3])) / (knot[4]-knot[1])^2;
    real u3 = (term2 - (term3*(knot[4]-knot[2]) - term4*(knot[3]-knot[2]))/(knot[4]-knot[3])) / (knot[4]-knot[1])^2;
    return [u1, u2, u3]';
  }
}

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
  
  // --- spline-related data
  vector[4] knot;
  vector[3] u_means;
  vector[3] u_sds;
  matrix[3,3] splineRot;

  int<lower=0> n_grid;
  vector[n_grid] t_grid;
}

parameters {
  real ageAtOnsetInt;
  real<lower=1e-6> sigma_ageAtOnset;
  vector<lower=censor_age>[n_cens_onset] onset_cens;
  vector<lower=1e-6>[2] sd_alpha;
  cholesky_factor_corr[2] Lcorr;
  matrix[2, n_subj] z_alpha;
  real<lower=1e-6> sd_mutation;
  real<lower=1e-6> sd_aoo_mutation;
  vector[n_mutation] z_mut;
  vector[n_mutation] z_aoo_mut;
  real alpha0;
  vector[3] b;
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
  for (j in 1:n_obs_onset) {
    ageAtOnset[subj_obs_onset[j]] = onset_obs[j];
  }
  for (j in 1:n_cens_onset) {
    ageAtOnset[subj_cens_onset[j]] = onset_cens[j];
  }
}

model {
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

  for (j in 1:n_obs_onset) {
    int i = subj_obs_onset[j];
    onset_obs[j] ~ normal(ageAtOnsetInt + m_aoo[mutation_by_sub[i]], sigma_ageAtOnset);
  }
  for (j in 1:n_cens_onset) {
    int i = subj_cens_onset[j];
    onset_cens[j] ~ normal(ageAtOnsetInt + m_aoo[mutation_by_sub[i]], sigma_ageAtOnset);
  }
  for (i in 1:n_parents) {
    parentOnset[i] ~ normal(ageAtOnsetInt + m_aoo[mutation_by_parent[i]], sigma_ageAtOnset);
  }

  // longitudinal model with cubic spline fixed effect
  for (i in 1:n_obs) {
    real t_i = age[i] - ageAtOnset[id[i]];
    vector[3] u_raw = spline_basis(t_i, knot);
    vector[3] u_std = (u_raw - u_means) ./ u_sds;
    vector[3] v_i = splineRot' * u_std;
    real mu_i = alpha0
                + m[mutation[i]]
                + dot_product(v_i, b)
                + alpha[1, id[i]]
                + alpha[2, id[i]] * t_i;
    y[i] ~ normal(mu_i, sdYRes);
  }
}

generated quantities {
  corr_matrix[2] Cor_alpha;
  Cor_alpha = multiply_lower_tri_self_transpose(Lcorr);

  // back-transform spline coefficients to original scale
  vector[3] b_scaled = splineRot * b;
  vector[3] b_orig = b_scaled ./ u_sds;
  real alpha0_orig = alpha0 - dot_product(b_orig, u_means);

  // population-average predicted curve vs time since onset (no mutation/subject REs)
  vector[n_grid] y_pred;
  for (i in 1:n_grid) {
    vector[3] u_raw_g = spline_basis(t_grid[i], knot);
    vector[3] u_std_g = (u_raw_g - u_means) ./ u_sds;
    vector[3] v_g = splineRot' * u_std_g;
    y_pred[i] = alpha0 + dot_product(v_g, b);
  }
}
