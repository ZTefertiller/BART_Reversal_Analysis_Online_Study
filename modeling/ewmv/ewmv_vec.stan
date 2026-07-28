functions {
  real partial_log_lik(array[] int slice_subj_idx,
                       int start,
                       int end,
                       // Other data and parameters
                       int ntrial,
                       int maxpump,
                       array[,] int npumps,
                       array[,] int outcome,
                       array[,,] int d,
                       vector phi,
                       vector eta,
                       vector rho,
                       vector tau,
                       vector lambda) {
    real ll = 0;

    for (i in 1:size(slice_subj_idx)) {
      int j = slice_subj_idx[i];  // Use the index from the slice
      int n_succ = 0;
      int n_pump = 0;
      real p_burst = phi[j];

      for (k in 1:ntrial) {
        real u_gain = 1;
        real u_stop = 0;
        real u_loss, u_pump, delta_u;

        int max_decision = npumps[j, k] + 1 - outcome[j, k];

        for (l in 1:max_decision) {
          if (d[j, k, l] != 75) {
            u_loss = l - 1;
            u_pump = (1 - p_burst) * u_gain
                      - lambda[j] * p_burst * u_loss
                      + rho[j] * p_burst * (1 - p_burst) * square(u_gain + lambda[j] * u_loss);
            delta_u = u_pump - u_stop;

            ll += bernoulli_logit_lpmf(d[j, k, l] | tau[j] * delta_u);
          }
        }

        n_succ += npumps[j, k] - outcome[j, k];
        n_pump += npumps[j, k];

        if (n_pump > 0) {
          p_burst = phi[j] + (1 - exp(-n_pump * eta[j])) * ((n_pump - n_succ) * 1.0 / n_pump - phi[j]);
        }
      }
    }

    return ll;
  }
}

data {
  int<lower=1> nsub;
  array[nsub] int<lower=1, upper=nsub> subj_idx;  // Make sure indexes are valid
  int<lower=1> ntrial;
  int<lower=1> maxpump;
  array[nsub, ntrial] int<lower=0> npumps;
  array[nsub, ntrial] int<lower=0, upper=1> outcome;
  array[nsub, ntrial, maxpump] int<lower=-1> d;
}

parameters {
  vector[5] mu_pr;
  vector<lower=0>[5] sigma;

  vector[nsub] phi_pr;
  vector[nsub] eta_pr;
  vector[nsub] rho_pr;
  vector[nsub] tau_pr;
  vector[nsub] lambda_pr;
}

model {
  mu_pr ~ normal(0, 1);
  sigma ~ normal(0, 0.2);

  phi_pr ~ normal(0, 1);
  eta_pr ~ normal(0, 1);
  rho_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  lambda_pr ~ normal(0, 1);

  {
    vector[nsub] phi = Phi_approx(mu_pr[1] + sigma[1] .* phi_pr);
    vector[nsub] eta = Phi_approx(mu_pr[2] + sigma[2] .* eta_pr);
    vector[nsub] rho = 0.5 - Phi_approx(mu_pr[3] + sigma[3] .* rho_pr);
    vector[nsub] tau = exp(mu_pr[4] + sigma[4] .* tau_pr);
    vector[nsub] lambda = exp(mu_pr[5] + sigma[5] .* lambda_pr);
    
    // Create array of indices 1:nsub for reduce_sum
    array[nsub] int indices;
    for (i in 1:nsub) {
      indices[i] = i;
    }
    
    target += reduce_sum(partial_log_lik, indices, 1, 
                         ntrial, maxpump, npumps, outcome, d,
                         phi, eta, rho, tau, lambda);
  }
}

generated quantities {
  real<lower=0> mu_phi = Phi(mu_pr[1]);
  real<lower=0> mu_eta = Phi(mu_pr[2]);
  real<lower=-0.5, upper=0.5> mu_rho = 0.5 - Phi(mu_pr[3]);
  real<lower=0> mu_tau = exp(mu_pr[4]);
  real<lower=0> mu_lambda = exp(mu_pr[5]);

  vector<lower=0, upper=1>[nsub] phi;
  vector<lower=0>[nsub] eta;
  vector<lower=-0.5, upper=0.5>[nsub] rho;
  vector<lower=0>[nsub] tau;
  vector<lower=0>[nsub] lambda;

  array[nsub, ntrial] real log_lik;

  for (j in 1:nsub) {
    phi[j] = Phi_approx(mu_pr[1] + sigma[1] * phi_pr[j]);
    eta[j] = Phi_approx(mu_pr[2] + sigma[2] * eta_pr[j]);
    rho[j] = 0.5 - Phi_approx(mu_pr[3] + sigma[3] * rho_pr[j]);
    tau[j] = exp(mu_pr[4] + sigma[4] * tau_pr[j]);
    lambda[j] = exp(mu_pr[5] + sigma[5] * lambda_pr[j]);

    int n_succ = 0;
    int n_pump = 0;
    real p_burst = phi[j];

    for (k in 1:ntrial) {
      real u_gain = 1;
      real u_stop = 0;
      real u_loss, u_pump, delta_u;
      log_lik[j, k] = 0;
      int max_decision = npumps[j, k] + 1 - outcome[j, k];

      for (l in 1:max_decision) {
        if (d[j, k, l] != 75) {
          u_loss = l - 1;
          u_pump = (1 - p_burst) * u_gain
                   - lambda[j] * p_burst * u_loss
                   + rho[j] * p_burst * (1 - p_burst) * square(u_gain + lambda[j] * u_loss);
          delta_u = u_pump - u_stop;

          log_lik[j,k] += bernoulli_logit_lpmf(d[j, k, l] | tau[j] * delta_u);
        }
      }

      n_succ += npumps[j, k] - outcome[j, k];
      n_pump += npumps[j, k];

      if (n_pump > 0) {
        p_burst = phi[j] + (1 - exp(-n_pump * eta[j])) * ((n_pump - n_succ) * 1.0 / n_pump - phi[j]);
      }
    }
  }
}
