functions {
  real partial_log_lik(array[] int slice_subj_idx,   // <-- rename to show it's a slice
                       int start, int end,           // (required args, but you won’t use them)
                       array[] int Tsubj,
                       array[,] int npumps,
                       array[,] int outcome,
                       array[,,] int d,
                       vector phi,
                       vector eta,
                       vector gam,
                       vector tau) {
    real ll = 0;
    
    //------------------------------------------------------------------
    // iterate over *slice_subj_idx* itself, not start:end
    //------------------------------------------------------------------
    for (i in 1:size(slice_subj_idx)) {
      int j = slice_subj_idx[i];
      int n_succ = 0;
      int n_pump = 0;

      for (k in 1:Tsubj[j]) {
        // Improved numerical stability in p_burst calculation
        real belief_component;
        real p_burst;
        real denom;
        real omega;
        
        // Prevent division by zero when n_pump is 0
        if (n_pump == 0) {
          belief_component = phi[j];
        } else {
          belief_component = (phi[j] + eta[j] * n_succ) / (1 + eta[j] * n_pump);
        }
        
        // Ensure belief_component is in valid range [0.001, 0.999]
        belief_component = fmin(0.999, fmax(0.001, belief_component));
        p_burst = 1 - belief_component;
        
        // Ensure p_burst is in valid range [0.001, 0.999]
        p_burst = fmin(0.999, fmax(0.001, p_burst));
        
        // Compute log(1-p_burst) with safety bounds
        denom = log1m(p_burst);
        denom = fmin(-0.001, denom);  // Stricter bound to avoid division issues
        
        // Calculate omega with bounds check
        omega = -gam[j] / denom;
        omega = fmin(1000, fmax(-1000, omega));  // Prevent extreme values
        
        int npump_eff = npumps[j, k] + 1 - outcome[j, k];

        if (npump_eff > 0) {
          array[npump_eff] int d_slice;
          vector[npump_eff] linpred;

          d_slice = d[j, k, 1:npump_eff];
          linpred = tau[j] * (omega - linspaced_vector(npump_eff, 1, npump_eff));
          linpred = fmin(100, fmax(-100, linpred));

          ll += bernoulli_logit_lpmf(d_slice | linpred);
        }

        n_succ += npumps[j, k] - outcome[j, k];
        n_pump += npumps[j, k];
      }
    }
    return ll;
  }
}

data {
  int<lower=1> nsub;            // Number of subjects
  int<lower=1> ntrial;          // Maximum number of trials
  array[nsub] int<lower=0> Tsubj;
  int<lower=1> maxpump;         // Maximum number of possible pumps
  array[nsub, ntrial] int<lower=0> npumps;
  array[nsub, ntrial] int<lower=0,upper=1> outcome;

}


transformed data {
  array[nsub, ntrial, maxpump] int d;
  array[nsub] int subj_idx = linspaced_int_array(nsub, 1, nsub);
  vector[maxpump] pump_indices = linspaced_vector(maxpump, 1, maxpump);

  for (j in 1:nsub) {
    for (k in 1:Tsubj[j]) {
      for (l in 1:maxpump) {
        d[j, k, l] = (l <= npumps[j, k]) ? 1 : 0;
      }
    }
  }
}

parameters {
  // Group-level parameters
  vector[4] mu_pr;
  vector<lower=0>[4] sigma;

  // Normally distributed error for Matt trick
  vector[nsub] phi_pr;
  vector[nsub] eta_pr;
  vector[nsub] gam_pr;
  vector[nsub] tau_pr;
}

transformed parameters {
  // Subject-level parameters with Matt trick and tighter bounds
  vector<lower=0.01,upper=0.99>[nsub] phi;  // Tighter bounds to avoid edge cases
  vector<lower=0.01>[nsub] eta;             // Lower bound to prevent eta=0 edge case
  vector<lower=0.01>[nsub] gam;             // Lower bound to prevent gam=0 edge case
  vector<lower=0.01>[nsub] tau;             // Lower bound to prevent tau=0 edge case

  for (j in 1:nsub) {
    // Apply constraints more carefully
    phi[j] = 0.01 + 0.98 * Phi_approx(mu_pr[1] + sigma[1] * phi_pr[j]);  // Range [0.01, 0.99]
    eta[j] = 0.01 + exp(mu_pr[2] + sigma[2] * eta_pr[j]);
    gam[j] = 0.01 + exp(mu_pr[3] + sigma[3] * gam_pr[j]);
    tau[j] = 0.01 + exp(mu_pr[4] + sigma[4] * tau_pr[j]);
  }
}

model {
  // Prior distributions
  mu_pr  ~ normal(0, 1);
  sigma ~ normal(0, 0.2);
  phi_pr ~ normal(0, 1);
  eta_pr ~ normal(0, 1);
  gam_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);

  // Log-likelihood calculation using reduce_sum
  target += reduce_sum(
    partial_log_lik,
    subj_idx,    
    1,                        // grain size
    Tsubj, npumps, outcome, d,
    phi, eta, gam, tau
  );
}

// Generated quantities with improved stability
generated quantities {
  // Group-level means (reparameterized)
  real<lower=0.01, upper=0.99> mu_phi = 0.01 + 0.98 * Phi_approx(mu_pr[1]);
  real<lower=0.01> mu_eta = 0.01 + exp(mu_pr[2]);
  real<lower=0.01> mu_gam = 0.01 + exp(mu_pr[3]);
  real<lower=0.01> mu_tau = 0.01 + exp(mu_pr[4]);
  
  vector[nsub] phi_out = phi;   // subject-level burst prior
  vector[nsub] eta_out = eta;   // learning-rate
  vector[nsub] gam_out = gam;   // risk coefficient
  vector[nsub] tau_out = tau;   // inverse-temperature

  // Trial-level log-likelihoods for LOO/WAIC
  array[nsub, ntrial] real log_lik;

  // -------- trial-level model predictions (new) --------
  array[nsub, ntrial] int npumps_pred;   // number of pumps the model would make
  // -----------------------------------------------------

  for (j in 1:nsub) {
    int n_succ = 0;
    int n_pump = 0;

    for (k in 1:Tsubj[j]) {
      real belief_component;
      real p_burst;
      real omega;

      // Compute belief about burst probability with improved stability
      if (n_pump == 0) {
        belief_component = phi[j];
      } else {
        belief_component = (phi[j] + eta[j] * n_succ) / (1 + eta[j] * n_pump);
      }
      belief_component = fmin(0.999, fmax(0.001, belief_component));
      p_burst = 1 - belief_component;
      p_burst = fmin(0.999, fmax(0.001, p_burst));

      real denom = log1m(p_burst);
      denom = fmin(-0.001, denom);           // avoid division issues
      omega  = -gam[j] / denom;
      omega  = fmin(1000, fmax(-1000, omega));

      // ---------- log-likelihood (unchanged) ------------
      log_lik[j, k] = 0;
      for (l in 1:(npumps[j, k] + 1 - outcome[j, k])) {
        real raw_pred = tau[j] * (omega - l);
        real bounded_pred = fmin(100, fmax(-100, raw_pred));
        log_lik[j, k] += bernoulli_logit_lpmf(d[j, k, l] | bounded_pred);
      }
      // --------------------------------------------------

      // ---------- simulate model-predicted pumps --------
      {
        int pred = 0;
        for (l in 1:maxpump) {                       // pump until stop or maxpump
          real raw_pred = tau[j] * (omega - l);
          real bounded_pred = fmin(100, fmax(-100, raw_pred));
          if (bernoulli_rng(inv_logit(bounded_pred)) == 1) {
            pred += 1;
          } else {
            break;
          }
        }
        npumps_pred[j, k] = pred;                    // save prediction
      }
      // --------------------------------------------------

      // Update counters for next trial
      n_succ += npumps[j, k] - outcome[j, k];
      n_pump += npumps[j, k];
    }

    // Fill unused trials with zeros
    for (k in (Tsubj[j] + 1):ntrial) {
      log_lik[j, k]   = 0;
      npumps_pred[j, k] = 0;                         // no prediction for non-existent trial
    }
  }
}
