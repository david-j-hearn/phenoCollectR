functions {

  real log_p_stage_counts_given_t(int S, int ind, array[,] int counts, vector t_std, vector O_std, real sigma_std) {

    //log multinomial: remove constant terms -> sum of count * log(p) for each stage
    real logp = 0.0;
    real z = 0.0;
    real z1 = 0.0;
    real p;
    z = (t_std[ind] - O_std[2] ) / sigma_std;
    logp += counts[ind,1] * normal_lccdf(z | 0, 1); //Need to update this with fixed
    for(i in 2:S-1) {
        z = (t_std[ind] - O_std[i] ) / sigma_std;
        z1 = (t_std[ind] - O_std[i+1]) / sigma_std; 
    logp += counts[ind,i] * log_diff_exp(normal_lcdf(z | 0, 1), normal_lcdf(z1 | 0, 1));
    }
    z = (t_std[ind] - O_std[S] ) / sigma_std;
    logp += counts[ind,S] * normal_lcdf(z | 0, 1);
    return logp;
  }
}

//In this coding, stage 1 is the stage before the first full stage in the time period
//Stage 1 in this coding doesn't have an onset (or it is constant 0), but it does have a duration
//Stage 2 onset is the onset of the first full stage
data {
  int<lower=1> N;
  int<lower=2> S;	//Number of stages + 1 to include the times before the stage 2 onset as stage 1
  int<lower=1> K; //Number of covariates
  matrix[N, K] X_raw; //Covariate data

  vector[N] t_raw;  //observed collection times

  array[N, S] int stage_counts;  //int array for each individual of counts per stage

  real T_max; //Maximum possible observed time (minimum assumed to be 0)

  //Hyperparameters
  matrix[S-1,K] betaMeans;
  matrix<lower=0>[S-1,K] betaSDs;

  vector[S-1] anchorMeans;
  vector<lower=0>[S-1] anchorSDs;

  real<lower=0> sigmaMean;
  real<lower=0> sigmaSD;

}

transformed data {

  vector[N] t_std;
  matrix[N, K] X_std;

  real mean_t = mean(t_raw);
  real sd_t = sd(t_raw);

  vector[K] mean_X;
  vector[K] sd_X;

  real T_min_std = (0.0 - mean_t) / sd_t;	//times start at m_raw = 0 for stage 1, but this might not always be true... so probabilities calcualted for "below stage 2" for stage 1.

  // Standardize time
  for (n in 1:N)
    t_std[n] = (t_raw[n] - mean_t) / sd_t;

  // Standardize covariates
  for (k in 1:K) {
    mean_X[k] = mean(col(X_raw, k));
    sd_X[k]   = sd(col(X_raw, k));

    for (n in 1:N) {
      X_std[n, k] = (X_raw[n, k] - mean_X[k]) / sd_X[k];
    }
  }

  // Standardize hyperparameters

  matrix[S-1,K] betaMeans_std;
  matrix<lower=0>[S-1,K] betaSDs_std;

  vector[S-1] anchorMeans_std;
  vector<lower=0>[S-1] anchorSDs_std;

  real<lower=0> sigmaMean_std;
  real<lower=0> sigmaSD_std;

  for (s in 1:(S-1)) {
    anchorMeans_std[s] = anchorMeans[s] / sd_t;
    anchorSDs_std[s] = anchorSDs[s] / sd_t;
    for (k in 1:K) {
      betaMeans_std[s,k] = betaMeans[s,k] * sd_X[k] / sd_t;
      betaSDs_std[s,k] = betaSDs[s,k] * sd_X[k] / sd_t;
    }
  }

  sigmaMean_std = sigmaMean / sd_t;
  sigmaSD_std = sigmaSD / sd_t;
  print("beta hyper means: ", betaMeans_std);
  print("beta hyper sds: ", betaSDs_std);
  print("anchor hyper means: ", anchorMeans_std);
  print("anchor hyper SDs: ", anchorSDs_std);
  print("sigma hyper mean: ", sigmaMean_std);
  print("anchor hyper SD: ", sigmaSD_std);
}

parameters {

  // Duration models (standardized time units)
  vector<lower=0.0>[S-1] alpha_d_std;		//intercepts w/o last stage 
  matrix[S-1, K] beta_d_std;			//slopes for duration models  - does not include last stage 
  real<lower=0> sigma_std; 	//single sigma flower-level variability
  //vector[N] gamma;  //individual random effects (shifting onset per individual)
  //real<lower=0> sigma_gamma;  //st dev of individual random effects
}

model {

  //epsilon
  real epsilon = 1e-12;

  // Priors tuned for standardized space - centered so intercepts alpha are anchors
  alpha_d_std ~ normal(anchorMeans_std, anchorSDs_std);
  to_vector(beta_d_std) ~ normal(to_vector(betaMeans_std), to_vector(betaSDs_std));
  //sigma_std ~ normal(sigmaMean_std, sigmaSD_std);	
  sigma_std ~ normal(0, 1);	

  //Likelihoods
  vector[S] O_mean_std;
  vector[S-1] D_mean_std;

  for (n in 1:N) {
    for(s in 1:S-1) {

      D_mean_std[s] = alpha_d_std[s] + dot_product(to_vector(beta_d_std[s,]), X_std[n,]);
      if(D_mean_std[s] < 0.0) {
        D_mean_std[s] = epsilon;	//Hack to make non-0 probabilities during initiaton and warmup - stops after parameter values get closer to true ones for well-separated stages
      //better would be to scale response times to a reasonable level so that softplus is not in non-linear regime for durations
      }

//Might be worth reparameterizing durations in terms of proportion of total time
      if(s == 1) 
        O_mean_std[1] = T_min_std;  //somewhat arbitrary
      else 
        O_mean_std[s] = O_mean_std[s-1] + D_mean_std[s-1];
    }
    O_mean_std[S] = O_mean_std[S-1] + D_mean_std[S-1];

    target += log_p_stage_counts_given_t(S, n, stage_counts, t_std, O_mean_std, sigma_std);
  }
}

generated quantities {

  // Back-transformed parameters (original scale)

  //Duration models 	//don't need one for final stage, as it is determined by the other stages
  vector[S-1] alpha_d;
  matrix[S-1, K] beta_d;
  vector[S-1] anchor_d;

  //Onset models	//onset 1 is at baseline 0
  vector[S] alpha_o;
  matrix[S, K] beta_o;
  vector[S] anchor_o; //no stage 1 onset 

  // Transform sigma
  real sigma = sd_t * sigma_std;

  for (s in 1:(S-1)) {
    for (k in 1:K) {

      beta_d[s, k] = (sd_t / sd_X[k]) * beta_d_std[s,k];

      if(s==1) {
        beta_o[s, k] = 0.0;
      }
      else {
        beta_o[s, k] = beta_o[s-1,k] + beta_d[s-1,k];
      }
    }

    alpha_d[s] = sd_t * alpha_d_std[s] -  dot_product(beta_d[s,], mean_X) ;
    anchor_d[s] = alpha_d[s];

    if(s == 1) {
      anchor_o[1] = 0.0;
      alpha_o[1] = 0.0;
    }
    else {
      anchor_o[s] = anchor_o[s-1] + anchor_d[s-1];
      alpha_o[s] = alpha_o[s-1] + alpha_d[s-1];
    }
  }

  anchor_o[S] = anchor_o[S-1] + anchor_d[S-1];
  alpha_o[S] = alpha_o[S-1] + alpha_d[S-1];
  for (k in 1:K) {
    beta_o[S, k] = beta_o[S-1,k] + beta_d[S-1,k];
  }
}
