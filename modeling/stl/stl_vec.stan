//vectorized stl model from zhou et al
functions {
  real partial_log_lik(array[] int slice_subj_idx,
                       int start,
                       int end,
                       array[,] int outcome,
                       array[,] int npumps,
                       array[,] int opportunity,
                       array[,] int nmax,
                       array[,,] int d,
                       vector vwin,
                       vector vloss,
                       vector beta,
                       vector omegaone) {

    real ll = 0;
    int ntrial = dims(outcome)[2];

    for (i in 1:size(slice_subj_idx)) {
      int subj = slice_subj_idx[i];
      vector[ntrial] omega;
      int last = 0;

      for (k in 1:ntrial) {
        if (last == 0) {
          omega[k] = nmax[subj, k] * omegaone[subj];
        } else {
          real vloss_k = vloss[subj];
          real vwin_k = vwin[subj];
          if (outcome[subj, last] == 1) {
            omega[k] = nmax[subj,k] *
              (omega[last] / nmax[subj,last]) *
              (1 - vloss_k * (1 - npumps[subj, last] * 1.0 / nmax[subj,last]));
          } else {
            omega[k] = nmax[subj,k] *
              (omega[last] / nmax[subj,last]) *
              (1 + vwin_k  * (     npumps[subj, last] * 1.0 / nmax[subj,last]));
          }
        }
        last = k;
      }

      for (k in 1:ntrial) {
        real beta_k = beta[subj];
        int opp = opportunity[subj, k];

        if (opp > 0) {
          row_vector[opp] n_idx = linspaced_row_vector(opp, 1, opp);
          if (dims(d)[3] >= opp) {
            ll += bernoulli_logit_lpmf(d[subj, k, 1:opp] | -beta_k * (n_idx - omega[k]));
          }
        }
      }
    }

    return ll;
  }
}

data {
  int<lower=1> nsub;
  int<lower=1> ntrial;
  array[nsub, ntrial] int<lower=0, upper=1> outcome;
  array[nsub, ntrial] int<lower=0> npumps;
  array[nsub, ntrial] int<lower=1> opportunity;
  array[nsub, ntrial] int nmax;
  int<lower=1> maxpump;
  array[nsub, ntrial, maxpump] int d;
}

parameters {
  real<lower=0> mu_vwin;
  real<lower=0> mu_vloss;
  real<lower=0> mu_beta;
  real<lower=0, upper=1> mu_omegaone;
  
  vector<lower=0>[4] sigma;

  vector<lower=0, upper=1>[nsub] vwin;
  vector<lower=0, upper=1>[nsub] vloss;
  vector<lower=0, upper=3>[nsub] beta;
  vector<lower=0, upper=1>[nsub] omegaone;
}

model {
  mu_vwin ~ normal(0, 1);
  mu_vloss ~ normal(0, 1);
  mu_beta ~ normal(0, 1);
  mu_omegaone ~ normal(0, 1);
  
  sigma ~ inv_gamma(1, 1);

  for (i in 1:nsub) {
    vwin[i] ~ normal(mu_vwin, sigma[1]);
    vloss[i] ~ normal(mu_vloss, sigma[2]);
    beta[i] ~ normal(mu_beta, sigma[3]);
    omegaone[i] ~ normal(mu_omegaone, sigma[4]);
  }

  array[nsub] int subj_idx = linspaced_int_array(nsub, 1, nsub);
  target += reduce_sum(partial_log_lik, subj_idx, 1,
                       outcome, npumps, opportunity, nmax, d,
                       vwin, vloss, beta, omegaone);
}

generated quantities {
  array[nsub, ntrial] real omega_out;
  array[nsub, ntrial] real log_lik;
  array[nsub, ntrial] int npumps_pred;  // <-- added prediction array

  for (i in 1:nsub) {
    vector[ntrial] omega;
    int last = 0;

    for (k in 1:ntrial) {
      if (last == 0) {
        omega[k] = nmax[i, k] * omegaone[i];
      } else {
        if (outcome[i, last] == 1) {
          omega[k] = nmax[i,k] *
            (omega[last] / nmax[i,last]) *
            (1 - vloss[i] * (1 - npumps[i, last] * 1.0 / nmax[i,last]));
        } else {
          omega[k] = nmax[i,k] *
            (omega[last] / nmax[i,last]) *
            (1 + vwin[i]  * (     npumps[i, last] * 1.0 / nmax[i,last]));
        }
      }
      last = k;
    }

    for (k in 1:ntrial) {
      omega_out[i, k] = omega[k];
      int opp = opportunity[i, k];

      // simulate predicted number of pumps
      int pred_pumps = 0;
      for (n in 1:opp) {
        real p_continue = inv_logit(-beta[i] * (n - omega[k]));
        if (bernoulli_rng(p_continue) == 1) {
          pred_pumps += 1;
        } else {
          break;
        }
      }
      npumps_pred[i, k] = pred_pumps;

      if (opp > 0 && dims(d)[3] >= opp) {
        row_vector[opp] n_idx = linspaced_row_vector(opp, 1, opp);
        log_lik[i, k] = bernoulli_logit_lpmf(d[i, k, 1:opp] | -beta[i] * (n_idx - omega[k]));
      } else {
        log_lik[i, k] = 0.0;
      }
    }
  }
}
