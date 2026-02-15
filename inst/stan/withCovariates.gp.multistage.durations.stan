functions {
	// Log probability of stage at time t
	//The times, mus (mean stage onsets O), sigmas in min max transformed scale
	real log_p_y_given_t(int S, int stage, real t, vector O, vector sigma, int debug) {
	
		real log_p;

		if (debug) {
			print("t=", t, " stage=", stage,
			" O_stage=", O[stage], " sigma_stage=", sigma[stage]);
		}

	//Three cases
I	//	Case I: Stage S: spans the "start" boundary
	//		P(t>m, O_1>t) or P(t<M, t>=O_S) = 1 * (1 - CDF(O_1)) + 1 * CDF(O_S)
	//		P(A U B) = 1 - P(A' ∩ B') = 1 - P(O_S > t, O_1 <= t)
		if (stage == S) {
			real log_p_comp = normal_lccdf(t | O[S], sigma[S]) + normal_lcdf(t | O[1], sigma[1]);
			log_p = log1m_exp(log_p_comp);
			return log_p;
		}

	//	Case II: Internal stages
	//		P(t>=O_{stage}, t<O_{stage+1})
		if (stage >= 1 && stage < S) {
			log_p = normal_lcdf(t | O[stage], sigma[stage]) + normal_lccdf(t | O[stage+1], sigma[stage+1]);
			return log_p;
		}

	//	NA!!!: Case III: Last stage: wraparound to O_1
	//		P(t>=O_{S},t<O_1)
		//if (stage == S) {
			//log_p = normal_lcdf(t | O[S], sigma[S]) + normal_lccdf(t | O[1], sigma[1]);
			//return log_p;
		//}

	reject("Invalid stage: ", stage, " for S=", S);
	return negative_infinity();
	}

  // Functions for E_Kth_shift / E_Kth_approx remain unchanged
	real E_Kth_shift(int n, real sigma, int k) {
		real p = (k - pi() / 8) / (n - pi() / 4 + 1);
		if (p <= 0 || p >= 1) reject("Probability argument to inv_Phi out of bounds: ", p);
		return sigma * inv_Phi(p);
	}

	real E_Kth_approx(int n, real mu, real sigma, int k) {
		return mu + E_Kth_shift(n, sigma, k);
	}
}

data {
	
	//Checks
	int<lower=0, upper=10> debug;
	int<lower=0, upper=1> drop_ll;		//useful if testing prior predictives

	//Handle Extremes
	int<lower=0, upper=1> process_extremes;	//process first and last of the onset events
	int<lower=0> n;				//population size to process extremes

	//Sample size
	int<lower=0> N;				//number of observations

	//Response data
	vector[N] T_raw;			//observed times in the min max transformed scale (between 0 and 1)
	real T_min;				//minimum time on the original scale
	real T_max;				//maximum time on the original scale

	//Stage data
	int<lower=2> S;				//number of stages
	array[N] int<lower=1, upper=S> stage;	//stage designation for each individual

	//Covariate data
	int<lower=1> C;				//number of covariates
	matrix[N,C] X_raw;			//covariate data in the min max transformed scale (between 0 and 1)
	vector[C] mean_X_raw;			//covariate means in the min max transformed scale (bewteen 0 and 1)
	vector[C] min_X;			//min of onset covariates in original scale
	vector<lower=min_X>[C] max_X;		//max of onset covariates in original scale

	int<lower=0> priors;			//0 = flat prior, 1 or above = normal prior

	//Hyperparameters: means and SDs for each parameter in the min max transformed scale

	//	Sigma parameters for onsets
	vector<lower=0>[S] sigmaMean;
	vector<lower=0>[S] sigmaSD;
	
	//	Stage two linear model
	real anchor1Mean;
	real<lower=0> anchor1SD;
	row_vector[C] beta1Mean;
	row_vector<lower=0>[C] beta1SD;

	//	Linear model for stage durations (except stage 1, which is defined in terms of when stage 2 starts, so these start with stage two duration up to stage S duration
	vector[S-1] anchor_dMean;
	vector<lower=0>[S-1] anchor_dSD;
	matrix[S-1,C] beta_dMean;
	matrix<lower=0>[S-1,C] beta_dSD;

	}

transformed data {

	real<lower=0> T_range = T_max - T_min;
	vector<lower=0>[C] range_X;
	vector[C] mean_X;

	//Transform covariate means and ranges to original scale
	for (k in 1:C) {
		range_X[k] = max_X[k] - min_X[k];
		mean_X[k] = mean_X_raw[k] * range_X[k] + min_X[k];
	}
}

parameters {

	//Stage duration increments (covariate-dependent)
	vector<lower=0,upper=1>[S-1] anchor_d_raw;		// anchors for increments : indexed starting at stage 1, 2, 3, ..., S-1 so S-1 elements
	matrix[S-1,C] beta_d_raw;		// covariate slopes for increments

	//Covariate linear model for stage 1 (first stage after "start" of cycle)
	real<lower=0,upper=1> anchor1_raw;	//anchor for stage 1 onset (mean stage 1 onset)
	vector[C] beta1_raw;			//covariate slopes for stage 1 onset

	//Stage-specific SDs
	vector<lower=0.001>[S] sigma_raw;
}

transformed parameters {
	matrix[N,S] O_raw;		//individual-specific stage boundary means based on linear models, indexed from stage 1
	vector[S-1] alpha_d_raw;	//intercepts for durations
	real alpha1_raw;			//intercept for stage 1 onset

	//Compute intercepts
	alpha1_raw = anchor1_raw - dot_product(mean_X_raw, beta1_raw);

	//index 1 is stage 1 duration
	for(i in 1:S-1) {
		alpha_d_raw[i] = anchor_d_raw[i] - dot_product(mean_X_raw, row(beta_d_raw, i));
	}

	//Compute individual-specific stage boundaries
	for (j in 1:N) {
		// Stage 1 anchored to predicted mean
		O_raw[j,1] = alpha1 + dot_product(row(X_raw,j), beta1_raw);

		// Forward stages (1..S-1)
		for (i in 1:S-1) {
			real d_linear = alpha_d_raw[i] + dot_product(row(X_raw,j), row(beta_d_raw,i));
			O_raw[j,i+1] = O_raw[j,i+1] + softplus(d_linear);	//Onset from previous stage plus duration from previous stage determines onset of following stage
		}

		//// NA!!! Wraparound stage 1: last stage + increment
		//real d_linear_S = alpha_d_raw[S-1] + dot_product(row(X_raw,j), row(beta_d_raw,S-1));
		//O_raw[j,1] = O_raw[j,S] + softplus(d_linear_S);
	}
}

model {
	// Priors
	if (priors >= 1) {
		anchor1_raw ~ normal(anchor1Mean, anchor1SD);
		beta1_raw ~ normal(beta1Mean, beta1SD);
		anchor_d_raw ~ normal(anchor_dMean, anchor_dSD);
		to_vector(beta_d_raw) ~ normal(to_vector(beta_dMean), to_vector(beta_dSD));
		sigma_raw ~ normal(sigmaMean, sigmaSD);
	}

  	// Likelihood
	if (!drop_ll) {
		for (j in 1:N) {
			target += log_p_y_given_t(S, stage[j], T_raw[j], to_vector(row(O_raw, j)), sigma_raw, debug);
		}
	}
}

generated quantities {
	real anchor1;
	real alpha1;
	row_vector[C] beta1;

	real alpha1k1;
	real alpha1kN;
	row_vector[C] beta1k1;
	row_vector[C] beta1kN;

	vector[S-1] anchor_d;
	vector[S-1] alpha_d;
	matrix[S-1,C] beta_d;

	real[S-1] alpha_dk1;
	real[S-1] alpha_dkN;
	matrix[S-1, C] beta_dk1;
	matrix[S-1, C] beta_dkN;

	vector[S] sigma;
	
	// Transform parameters for duration increments and sigma to original scale
	for (i in 1:S-1) {
		for (c in 1:C) {
			beta_d[i,c] = (T_range / range_X[c]) * beta_d_raw[i,c];

			//Process extremes if flagged
			if(process_extremes) {
			beta_dk1[i,c] = beta_d[i,c];
			beta_dkN[i,c] = beta_d[i,c];
			}

		//Transform anchor and intercept to original scale
		anchor_d[i] = T_min + T_range * anchor_d_raw[i];
		alpha_d[i] = anchor_d[i] - dot_product(row(beta_d,i), mean_X);


		//Transform sigma to original scale
		sigma[i] = T_range * sigma_raw[i];
		}
	}
	
	//Get the last sigma!
	sigma[S] = T_range * sigma_raw[S];

	if(process_extremes) {
		for (i in 1:S-1) {
			alpha_dk1[i] = alpha_d[i] + E_Kth_shift( n, sigma[i], 1);   //intercept on first
			alpha_dkN[i] = alpha_d[i] + E_Kth_shift( n, sigma[i], n);   //intercept on last
			anchor_dk1[i] = E_Kth_approx(n, anchor_d[i], sigma[i], 1);  //anchor on first
			anchor_dkN[i] = E_Kth_approx(n, anchor_d[i], sigma[i], n);  //anchor on last
		}
	}

	// Transform parameters for stage 1 to original scale
	for (c in 1:C) {
		beta1[c] = (T_range / range_X[c]) * beta1_raw[c];
		if(process_extremes) {
			beta1k1[c] = beta1[c];
			beta1kn[c] = beta1[c];
		}
	}

	anchor1 = T_min + T_range * anchor1_raw;
	alpha1 = anchor1 - dot_product(beta1, mean_X);

	if(process_extremes) {
		alpha1dk1 = alpha2 + E_Kth_shift( n, sigma[1], 1);   //intercept on first
		alpha1dkN = alpha2 + E_Kth_shift( n, sigma[1], n);   //intercept on last
		anchor1dk1 = E_Kth_approx(n, anchor1, sigma[1], 1);  //anchor on first
		anchor1dkN = E_Kth_approx(n, anchor1, sigma[1], n);  //anchor on last
	}
}
