//NOTES:
//  Case I: If individual has sequential stages and no phenological units that develop on indiviudal, then all counts will be 0 except for the stage the individual is in, which will be 1
//  Case II: If individual has stages with 0 or more phenological units, 0 or more pre stages before phenological units have developed and 0 or more post stages after which phenological units have abscised:
//    If individual is in a pre or post stage, that stage gets a 1 and all other stages get a 0  - this might be problematic, as the missing units doesn't impact the likelihood at all -> drifting sampling
//      Set lower and upper bounds on expected units so that prior is proper
//    If individual has 1 or more phenological units, each stage with units gets a count of the units in that stage and all pre and post stages get a 0
//  If either the true number of units is known a priori, or if Case I is in action and expected counts are all set to 1, then prior variance is 0 using a fixed parameter:
//  parameters {
//real<lower=0> theta_free;  // sampled only if needed
//}

//transformed parameters {
//real theta;
//if (fix_theta == 1)
//theta = fixed_value;      // fixed
//else
//theta = theta_free;       // sampled
//} 
//MAKE SURE TO PLACE WELL-DEFINED UPPER AND LOWER ON SUCH PARAMETERS AND USE A WELL-DEFINED PRIOR, NOT A FLAT PRIOR - THIS IS LIKE PRIOR PREDICTIVE CHECKS but can lead to sampler trouble when flat

functions {

  //fractional multinomial giving equal total weight of 1 to each individual
  real log_p_stage_counts_given_t(
      int ind,            //Index of the individual to be processed
      int S,              //Precomputed total number of stages
      int preN,           //Number of "pre latent" stages before stages with visible phenological units
      int visN,           //Number of "visible" stages
      int postN,          //Number of "post latent" stages after stages with visible phenological units
      vector N,           //Precomputed total count of visible phenological units for each individual
      vector N_miss,      //Precomputed count of "latent" phenological units for each individual based on expected number of units and observed number of units. "Count" is a non-negative integer representing an expected value
      matrix counts,      //Matrix of observed counts for visible stages for all individuals
      array[] int stages,      //Precomputed stage numbers if in "latent" pre or post stage or 0 if visible stages present
      vector t_std,       //Times of observation of individuals
      vector O_std,       //Vector of parameters representing mean onset times of all stages (stage 1 is -infinity and not used for probability calculations)
      real sigma_std) {

    real z_lo;
    real z_hi;
    real log_p = 0.0;
    real log_p_miss = 0.0;

    //Check if in "latent" stage and return full density on that stage alone
    //  Stage is set to 0 to indicate there are visible units, as units can be spread over multiple stages and the individual is not "in" a particular stage, per se
    if(stages[ind] == 1) {       //individual in "latent" stage 1
      z_hi = (t_std[ind] - O_std[2]) / sigma_std; //sigma index 1 is stage 2, S-1 is S, stage-1 is stage
      return normal_lccdf(z_hi | 0, 1);     //1 * log(P)
    }
    else if(stages[ind] == S) { //individual in "latent" stage S
      z_lo = (t_std[ind] - O_std[S]) / sigma_std;
      return normal_lcdf(z_lo | 0, 1);      //1 * log(P)
    }
    else if(stages[ind]>0) {    //individual in "latent" stage greater than 1 and less than S
      z_lo = (t_std[ind] - O_std[stages[ind]]) / sigma_std;
      z_hi = (t_std[ind] - O_std[stages[ind]+1]) / sigma_std;
      return log_diff_exp(      //1 * log(P)
          normal_lcdf(z_lo | 0, 1),
          normal_lcdf(z_hi | 0, 1)
          );
    }
    else {                      //individual has phenological units, so not in "latent" stage, so must have 1 or more visible units
      for(i in (preN+1):(preN+visN)) {
        real prop = counts[ind,i]/(N[ind]+N_miss[ind]);
        //real prop = counts[ind,i];
        if(counts[ind,i]>0) {             //If 0 counts in stage, the stage contributes 0 weight to likelihood
          if(i==1) {                      //Case where there are 0 "pre latent" stages. Stage 1 has -infinity onset time
            z_hi = (t_std[ind] - O_std[2]) / sigma_std;     
            log_p += prop * normal_lccdf(z_hi | 0, 1);
          }
          else if(i==S) {                 //Case where there are 0 "post latent" stages. Stage S has infinity cessation time
            z_lo = (t_std[ind] - O_std[S]) / sigma_std;
            log_p += prop * normal_lcdf(z_lo | 0, 1);
          }
          else {                          //Stage has visible units and is between the first and last stages
            z_lo = (t_std[ind] - O_std[i]) / sigma_std;
            z_hi = (t_std[ind] - O_std[i+1]) / sigma_std;

            log_p += prop * log_diff_exp(
                normal_lcdf(z_lo | 0, 1),
                normal_lcdf(z_hi | 0, 1)
                );
          }
        }
      }

      //Calculate probability of phenological unit being in "latent" or "missing" stage where units are not visible and estimate "missing" probability contribution
      if(N_miss[ind]>0) {
        real prop = N_miss[ind]/(N[ind]+N_miss[ind]);
        //real prop = N_miss[ind];
        if(postN==0 && preN>0) {      //only vis and pre stages present
          real z_pre = (t_std[ind] - O_std[preN+1]) / sigma_std;
          log_p_miss = normal_lccdf(z_pre | 0, 1);
        }
        else if(preN==0 && postN>0) {  //only vis and post stages present
          real z_post = (t_std[ind] - O_std[preN+visN+1]) / sigma_std;
          log_p_miss = normal_lcdf(z_post | 0, 1);
        }
        else if(preN>0 && postN>0) {    //mix of vis, pre and post stages
          real z_pre = (t_std[ind] - O_std[preN+1]) / sigma_std;
          real z_post = (t_std[ind] - O_std[preN+visN+1]) / sigma_std;
          log_p_miss = log_sum_exp(normal_lccdf(z_pre | 0, 1),normal_lcdf(z_post | 0, 1));
        }
        else if(preN==0 && postN==0) {
          //pure visible stages
        }
        else {
          reject("Be sure the number of pre-visible stages, number of post-visible stages, and number of stages with visible units are correct.");
        }
        log_p += prop*log_p_miss;
      }
    }
    return(log_p);
  }
}

//In this coding, stage 1 is the stage before the first full stage in the time period
//Stage 1 in this coding doesn't have an onset, but it does have a duration
//Stage 2 onset is the onset of the first full stage
data {
  int<lower=1> N; //number of individuals
  int<lower=2> S; //number of stages total
  int<lower=0> preN;
  int<lower=0> visN;
  int<lower=0> postN;

  int<lower=1> K;                 //Number of covariates

  vector[N] t_raw;
  matrix[N, K] X_raw;

  matrix[N, S] stage_counts;   //Counts of units in each stage. If no visible units, then a 1 for the stage the individual is in and 0 for other stages. For some cases, this is best treated as real rather than int.

  vector<lower=0>[N] minimum_counts;   //minimum total counts of units per individual, including latent pre-developed or post-abscised units, used as hyperparameters if not fixed.

  int<lower=0, upper=1>fixed;     //treat expected_counts as fixed parameters known a priori

  real<lower=1.0> scale;                     //scale to use to transform data

  real<lower=0> T_min;                     //Minimum observation time possible (e.g., 0 for day of year)
  real<lower=T_min> T_max;                 //Maximum observation time possible (e.g., 365 for day of year)

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

  //int<lower=2> S = preN + visN + postN;
  int onesCnt = 0;

  array[N] int stages = rep_array(0,N);     //stages are set to 0 if visible counts are present, otherwise the pre or post stage is recorded for the individual to reduce likelihood computations.
  vector[N] totN = rep_vector(0,N);
  vector[N] missN_obs = rep_vector(0,N);    //mean prior for latent units

  if(T_min !=0) {
    reject("Currently, the minimum response time must be set to 0. Please translate your data accordingly.");
  }

  mean_t = T_min;
  sd_t = T_max / scale;

  real T_min_std = (T_min - mean_t) / sd_t;	//times start at m_raw = 0 for stage 1

  // Standardize time
  for (n in 1:N) {
    t_std[n] = (t_raw[n] - mean_t) / sd_t;
  }

  //Set the stage if in a pre-unit stage so that unnecessary calculations are not performed
  for (n in 1:N) {
    onesCnt=0;
    if(preN>0) {
      for(i in 1:preN) {
        totN[n] += stage_counts[n,i];
        if(stage_counts[n,i] == 1) {
          onesCnt = onesCnt+1;
          if(onesCnt > 1) {
            reject("Only 1 or fewer pre or post stages can be assigned a 1 for an individual in the pre or post stage. All other stages must be assigned a 0.");
          }
          stages[n] = i;
        }
      }
      if(postN>0) {
      }
      //Set the stage if in a post-unit stage so that unnecessary calculations are not performed
      for(i in (preN+visN+1):S) {
        totN[n] += stage_counts[n,i];
        if(stage_counts[n,i] == 1) {
          onesCnt = onesCnt+1;
          if(onesCnt > 1) {
            reject("Only 1 or fewer pre or post stages can be assigned a 1 for an individual in the pre or post stage. All other stages must be assigned a 0.");
          }
          stages[n] = i;
        }
      }
    }
    if(visN>0) {
      for(i in (preN+1):(preN+visN)) {
        totN[n] += stage_counts[n,i];
        if(onesCnt==1 && stage_counts[n,i]>0) {
          reject("If an individual is assigned to a pre or post stage, there can be no visible units developed on the individual, so stage_counts in the stages where units are visible should all be 0.");
        }
      }
    }

  }
  //Check if a pre or post stage is set and visible unit counts are also present. Should not occur
  for(i in 1:N) {
    //missN_obs[i] = fmax(minimum_counts[i] - totN[i], 1e-6); //make sure missing count is positive.
    missN_obs[i] = fmax(minimum_counts[i] - totN[i], 0.5); //make sure missing count is positive.
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
  //vector<lower=0.0>[N] missN_sampled;
  //real<lower=0> phi;
}

//transformed parameters {
  //vector<lower=0.0>[N] missN;
  //real<lower=0> varG = 1.0;

  //if(fixed==1) {
    //missN = missN_obs;
  //}
  //else {
    //missN = missN_sampled;
    //varG = 10.0;
  //}
//}

model {

  //epsilon
  real epsilon = 1e-12;

  // Priors tuned for standardized space
  alpha_d_std ~ normal(anchorMeans_std, anchorSDs_std);
  to_vector(beta_d_std) ~ normal(to_vector(betaMeans_std), to_vector(betaSDs_std));
  sigma_std ~ normal(sigmaMean_std, sigmaSD_std);	

  //phi ~ exponential(1);

  //missN_sampled ~ gamma(rep_vector(phi,N), rep_vector(phi,N) / missN_obs); //positive number of unobserved units yet to develop or abscised.
  //missN_sampled ~ gamma(phi, phi ./ missN_obs); //positive number of unobserved units yet to develop or abscised.
  //missN_sampled ~ gamma(square(missN_obs) / varG, missN_obs / varG); //positive number of unobserved units yet to develop or already abscised.

  //Likelihoods
  vector[S] O_mean_std;
  vector[S-1] D_mean_std;

  for (n in 1:N) {
    for(s in 1:S-1) {

      D_mean_std[s] = alpha_d_std[s] + dot_product(to_vector(beta_d_std[s,]), X_std[n,]);
      if(D_mean_std[s] < 0.0)
        D_mean_std[s] = epsilon;	//Hack to make non-0 probabilities during initiaton and warmup - stops after parameter values get closer to true ones for well-separated stages
                                  //Better: scale times so that sampler initialization is feasible, but softplus is in linear range

      if(s == 1) 
        O_mean_std[1] = T_min_std;
      else 
        O_mean_std[s] = O_mean_std[s-1] + D_mean_std[s-1];

    }
    O_mean_std[S] = O_mean_std[S-1] + D_mean_std[S-1];

    target += log_p_stage_counts_given_t(
        n,
        S,
        preN,
        visN,
        postN,
        totN,
        //missN,
        missN_obs,
        stage_counts,
        stages,
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
