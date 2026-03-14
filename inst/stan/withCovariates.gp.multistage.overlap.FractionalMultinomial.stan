functions {

  real log_p_stage_counts_given_t(
      int S,
      int ind,
      matrix counts,
      vector t_std,
      vector O_std,
      real sigma_std) {

    real z_lo;
    real z_hi;
    real logp = 0.0;

    // 
    // P(stage=1 | t, theta) = P(t<O2)
    if(counts[ind,1]>0) {
      z_hi = (t_std[ind] - O_std[2]) / sigma_std; //sigma index 1 is stage 2, S-1 is S, stage-1 is stage
      logp += counts[ind,1] * normal_lccdf(z_hi | 0, 1);
    }

    for(i in 2:(S-1)) {
      // Internal stage 
      // P(stage=i | t, theta) = P(O_i <= t < O_{i+1})
      // use de morgan's law: P(stage=i) = CDF(t|O_i) - CDF(t|O_{i+1}) 
      if(counts[ind,i]>0) {
          z_lo = (t_std[ind] - O_std[i]) / sigma_std;
          z_hi = (t_std[ind] - O_std[i+1]) / sigma_std;

          logp += counts[ind,i] * log_diff_exp(
              normal_lcdf(z_lo | 0, 1),
              normal_lcdf(z_hi | 0, 1)
              );
      }
    }

          // Final stage - final stage is the "after" part
          // Stage S is the final stage after the onset of the last stage and before the end of the time period
          // P(stage=S | t, theta) = P(t>OS)
     if(counts[ind,S]>0) {
       z_lo = (t_std[ind] - O_std[S]) / sigma_std;
       logp += counts[ind,S] * normal_lcdf(z_lo | 0, 1);
     }
  return(logp);
}

}

//In this coding, stage 1 is the stage before the first full stage in the time period
//Stage 1 in this coding doesn't have an onset, but it does have a duration
//Stage 2 onset is the onset of the first full stage
data {
  int<lower=1> N;
  int<lower=2> S;	//this is the number of stages + 1 to include the times before the stage 2 onset as stage 1
  int<lower=1> K;

  vector[N] t_raw;
  matrix[N, K] X_raw;

  array[N, S] int stage_counts;

  real T_max;

  //Hyperparameters
  matrix[S-1,K] betaMeans;
  matrix<lower=0>[S-1,K] betaSDs;

  vector[S-1] anchorMeans;
  vector<lower=0>[S-1] anchorSDs;

  real<lower=0> sigmaMean;
  real<lower=0> sigmaSD;

  //int<lower=0,upper=1> ppd;
  //int<lower=1> nXs;
  //int<lower=1> nReps;
  //matrix[nXs*nReps,K] xPPD;
}

transformed data {

  vector[N] t_std;
  matrix[N, K] X_std;

  real mean_t = mean(t_raw);
  real sd_t = sd(t_raw);

  vector[K] mean_X;
  vector[K] sd_X;

  matrix<lower=0>[N, S] stage_proportions;
  vector[S] row;

  real T_min_std = (0.0 - mean_t) / sd_t;	//times start at m_raw = 0 for stage 1

  // Standardize time
  for (n in 1:N)
    t_std[n] = (t_raw[n] - mean_t) / sd_t;

  for (n in 1:N) {
    row = to_vector(stage_counts[n,]);
    stage_proportions[n,] = (row / sum(row))';
    //print("t_raw: ", t_raw[n], "; t_std: ", t_std[n], "; orig counts: ", row, "; proportions: ", stage_proportions[n,], "; t mean: ", mean_t, "; sd t: ", sd_t);
  }


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
  print("beta hyperparameter means: ", betaMeans_std);
  print("beta hyperparameter sds: ", betaSDs_std);
  print("anchor hyperparameter means: ", anchorMeans_std);
  print("anchor hyperparameter SDs: ", anchorSDs_std);
  print("sigma hyperparameter mean: ", sigmaMean_std);
  print("sigma hyperparameter SD: ", sigmaSD_std);

}

parameters {
  // Duration models (standardized time units)
  vector<lower=0.0>[S-1] alpha_d_std;		//does not include last stage 
  matrix[S-1, K] beta_d_std;			//does not include last stage 
                                  //vector<lower=0.001>[S-1] sigma_std; 	//does not include stage 1, which is fixed at m
  real<lower=0> sigma_std; 	//single sigma
  //vector[N] gamma; //individual random effects (shifting onset per individual)
  //real<lower=0> sigma_gamma; //st dev of individual random effects
}

model {

  //epsilon
  real epsilon = 1e-12;

  // Priors tuned for standardized space
  alpha_d_std ~ normal(anchorMeans_std, anchorSDs_std);
  to_vector(beta_d_std) ~ normal(to_vector(betaMeans_std), to_vector(betaSDs_std));
  sigma_std ~ normal(sigmaMean_std, sigmaSD_std);	

  //Likelihoods
  vector[S] O_mean_std;
  vector[S-1] D_mean_std;

  for (n in 1:N) {
    for(s in 1:S-1) {

      D_mean_std[s] = alpha_d_std[s] + dot_product(to_vector(beta_d_std[s,]), X_std[n,]);
      if(D_mean_std[s] < 0.0)
        D_mean_std[s] = epsilon;	//Hack to make non-0 probabilities during initiaton and warmup - stops after parameter values get closer to true ones for well-separated stages

      if(s == 1) 
        O_mean_std[1] = T_min_std;
      else 
        O_mean_std[s] = O_mean_std[s-1] + D_mean_std[s-1];

    }
    O_mean_std[S] = O_mean_std[S-1] + D_mean_std[S-1];

    target += log_p_stage_counts_given_t(S,
        n,
        stage_proportions,
        t_std,
        O_mean_std,
        sigma_std);
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
  vector[S] anchor_o; 

  // Transform sigma
  //vector[S] sigma; 	//one for each onset, stage 1 has sigma 0
  real sigma = sd_t * sigma_std;

  //matrix[nXs*nReps, S] y_pred = rep_matrix(0.0, nXs*nReps, S); 

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

    //sigma[s] = sd_t*sigma_std[s-1]
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

  //sigma[S] = sd_t * sigma_std[S-1];
  anchor_o[S] = anchor_o[S-1] + anchor_d[S-1];
  alpha_o[S] = alpha_o[S-1] + alpha_d[S-1];
  for (k in 1:K) {
    beta_o[S, k] = beta_o[S-1,k] + beta_d[S-1,k];
  }

  //if(ppd == 1) {
    //for (n in 1:(nXs*nReps) ) {
      //y_pred[n,1] = 0.0;	//Should be identically 0
      //for(s in 2:S) {
        ////if(s == 1) {
        ////y_pred[n,s] = alpha_o[s] + dot_product(beta_o[s,], xPPD[n,]);
        ////}
        ////else {
        //y_pred[n,s] = normal_rng(alpha_o[s] + dot_product(beta_o[s,], xPPD[n,]),sigma);
        ////}
      //}
    //}
  //}
}
