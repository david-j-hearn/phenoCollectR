functions {

// Log probabilities for y in {before, during, after} at time t
// given onset ~ Normal(mu, sigma) and duration D > 0
	real log_p_y_given_t(int stage, real t, real mu, real sigma, real D, int debug) {

		if(debug) {
			print("time, onset, sigma, duration: ", [t, mu, sigma, D]);
		}


		// log P(before) = log P(S > t)
		if(stage == 1) {
			return normal_lccdf(t | mu, sigma);
		}

		// log P(after) = log P(S + D < t) = log P(S < t - D)
		if(stage == 3) {
			return normal_lcdf(t - D | mu, sigma);
		}

		// log P(during) = log( Phi((t-mu)/sigma) - Phi((t-D-mu)/sigma) )
		real log_cdf0 = normal_lcdf(t | mu, sigma);
		real log_cdf1 = normal_lcdf(t - D | mu, sigma);
		return log_diff_exp(log_cdf0, log_cdf1);
	}

	real E_Kth_shift( int n, real sigma, int k) {
		real p = (k - pi() / 8) / (n - pi() / 4 + 1);
		if (p <= 0 || p >= 1) {
			reject("Probability argument to inv_Phi out of bounds: ", p);
		}
		real shift = sigma * inv_Phi(p);
		return(shift);
	}

	real E_Kth_approx( int n, real mu, real sigma, int k) {
		real approx;
		real shift = E_Kth_shift(n, sigma, k);
		approx = mu + shift;
		return(approx);
	}
}

data {

	int<lower=0, upper=10> debug;
	int<lower=0, upper=1> drop_ll;	//useful if testing prior predictives

	int<lower=0, upper=1> process_extremes;	//process first onset / last cessation
	int<lower=0> n;	//population size to process extremes

	//Sample size
	int<lower=0> N; //number of observations

	//Number of covariates for onset (K_O) and for duration (K_D)
	int<lower=1> K_O;
	int<lower=1> K_D;
	int<lower=1> K_union;	//used for cessation and observed (T) models

	int<lower=0, upper=2> priors; //0 = flat prior, 1=population parameter prior, 2=individual parameter prior

	array[K_O] int<lower=1, upper=K_union> idx_O;
	array[K_D] int<lower=1, upper=K_union> idx_D;

	vector[N] T_raw;  //observed times in the min max transformed scale (between 0 and 1)
	real T_min; //minimum time on the original scale
	real T_max; //maximum time on the original scale

	//Stage data
	array[N] int<lower=1, upper=3> stage;

	//Covariate data
	matrix[N,K_O] X_O_raw;          //covariate data in the min max transformed scale (between 0 and 1)
	matrix[N,K_D] X_D_raw;          //covariate data in the min max transformed scale (between 0 and 1)

	vector[K_O] mean_X_O_raw;           //covariate means in the min max transformed scale (bewteen 0 and 1)
	vector[K_D] mean_X_D_raw;           //covariate means in the min max transformed scale (bewteen 0 and 1)

	vector[K_O] min_X_O;            //min of onset covariates in original scale
	vector<lower=min_X_O>[K_O] max_X_O;         //max of onset covariates in original scale

	vector[K_D] min_X_D;            //min of duration covariates in original scale
	vector<lower=min_X_D>[K_D] max_X_D;         //max of duration covariates in original scale

	//Hyperparameters: means and SDs for each parameter in the min max transformed scale
	vector[K_O] betaOnsetMeans;
	vector<lower=0>[K_O] betaOnsetSDs;

	real anchorOnsetMean;
	real<lower=0> anchorOnsetSD;

	vector[K_D] betaDurationMeans;
	vector<lower=0>[K_D] betaDurationSDs;

	real<lower=0> anchorDurationMean;
	real<lower=0> anchorDurationSD;

	//Hyperparameters mean and sd of sigma
	real<lower=0> sigmaMean;
	real<lower=0> sigmaSD;

}

transformed data {

	real epsilon = 1e-6;
	vector[N] mu_D_raw_prior;
	real alpha_D_raw_prior;

	vector[N] mu_O_raw_prior;
	real alpha_O_raw_prior;

	print("Mus: ", [anchorDurationMean*365.0, anchorOnsetMean*365]);
	print("SDs: ", [anchorDurationSD*365.0, anchorOnsetSD*365]);
	print("Mus raw: ", [anchorDurationMean, anchorOnsetMean]);
	print("SDs raw: ", [anchorDurationSD, anchorOnsetSD]);

	real<lower=0> T_range = T_max - T_min;

	vector<lower=0>[K_O] range_X_O; //covariate ranges in the original scale
	vector<lower=0>[K_D] range_X_D; //covariate ranges in the original scale

	vector[K_O] mean_X_O; //mean covariate values in the original scale
	vector[K_D] mean_X_D; //mean covariate values in the original scale

	for (k in 1:K_O) {
		range_X_O[k] = max_X_O[k] - min_X_O[k];
		mean_X_O[k] = mean_X_O_raw[k] * range_X_O[k] + min_X_O[k];
	}

	for (k in 1:K_D) {
		range_X_D[k] = max_X_D[k] - min_X_D[k];
		mean_X_D[k] = mean_X_D_raw[k] * range_X_D[k] + min_X_D[k];
	}
}

parameters {
	//slope parameters, beta, and mean response, (anchor) at mean covariate value, and onset and cessation distribution standard deviation, sigma

	vector[K_O] beta_O_raw;
	real<lower=0,upper=1> anchor_O_raw;

	vector[K_D] beta_D_raw;
	real<lower=0,upper=1-anchor_O_raw> anchor_D_raw;

	real<lower=0.001> sigma_raw;
}

transformed parameters {
	real alpha_O_raw; //intercept of the linear model for Onset
	real alpha_D_raw; //intercept of the linear model for Duration

	vector[N] mu_O_raw;
	vector[N] mu_D_raw;
	vector[N] mu_C_raw;

	alpha_O_raw = anchor_O_raw - dot_product(mean_X_O_raw , beta_O_raw);
	alpha_D_raw = anchor_D_raw - dot_product(mean_X_D_raw , beta_D_raw);

	//Calculate mean onset in transformed scale
	mu_O_raw =  alpha_O_raw + X_O_raw * beta_O_raw;

	//Calculate mean duration in transformed scale, using softplus to assure positivity
	mu_D_raw =  alpha_D_raw + X_D_raw * beta_D_raw;

	//Calculate mean cessation based on onset and duration
	mu_C_raw = mu_O_raw + mu_D_raw;
}

model {
	if(debug) {
	}

	//Calculate the likelihood
	if(!drop_ll) {
		for(i in 1:N) {
			target += log_p_y_given_t(stage[i], T_raw[i], mu_O_raw[i], sigma_raw, mu_D_raw[i], debug);
		}
	}

	if(priors >= 1) {
		beta_O_raw ~ normal( betaOnsetMeans, betaOnsetSDs);
		beta_D_raw ~ normal( betaDurationMeans, betaDurationSDs);

		anchor_O_raw ~ normal( anchorOnsetMean, anchorOnsetSD);
		anchor_D_raw ~ normal( anchorDurationMean, anchorDurationSD);

		sigma_raw ~ normal( sigmaMean, sigmaSD);
	}
}

generated quantities {
	//calculate parameter values in the original scale
	real anchor_O;
	real<lower=0> anchor_D;
	real anchor_C;
	real anchor_T;
	real anchorOk1;		//remove underscore to avoid matching anchor_O
	real anchorCkN;
	real alpha_O;
	real alpha_D;
	real alpha_C;
	real alpha_T;
	real alphaOk1;
	real alphaCkN;
	vector[K_O] beta_O;
	vector[K_D] beta_D;
	vector[K_union] beta_C;
	vector[K_union] beta_T;
	vector[K_O] betaOk1;
	vector[K_union] betaCkN;

	real<lower=0> sigma;

	beta_C = rep_vector(0, K_union);
	beta_T = rep_vector(0, K_union);
	betaCkN = rep_vector(0, K_union);

	// Transform beta to original scale
	for (k in 1:K_O) {
		beta_O[k] = (T_range / range_X_O[k]) * beta_O_raw[k];
		beta_C[idx_O[k]] += beta_O[k];
		beta_T[idx_O[k]] += beta_O[k];
		if(process_extremes) {
			betaOk1[k] = beta_O[k];
			betaCkN[idx_O[k]] += beta_O[k];
		}
	}

	for (k in 1:K_D) {
		beta_D[k] = (T_range / range_X_D[k]) * beta_D_raw[k];
		beta_C[idx_D[k]] += beta_D[k];
		beta_T[idx_D[k]] += beta_D[k]/2;
		if(process_extremes) {
			betaCkN[idx_D[k]] += beta_D[k];
		}
	}


	// Transform anchor to original scale
	anchor_O = T_min + T_range * anchor_O_raw;
	anchor_D = T_range * anchor_D_raw;      //duration anchor is in (0, T_max - T_min)
	anchor_C = anchor_O + anchor_D;
	anchor_T = anchor_O + anchor_D/2;

	// Compute intercept at original scale
	alpha_O = anchor_O - dot_product(beta_O, mean_X_O);
	alpha_D = anchor_D - dot_product(beta_D, mean_X_D);
	alpha_C = alpha_O + alpha_D;
	alpha_T = alpha_O + alpha_D/2;

	// Transform sigma to original scale
	sigma = T_range * sigma_raw;

	if(process_extremes) {
		alphaOk1 = alpha_O + E_Kth_shift( n, sigma, 1);	//intercept on first onset
		alphaCkN = alpha_C + E_Kth_shift( n, sigma, n);	//intercept on last cessation
		anchorOk1 = E_Kth_approx(n, anchor_O, sigma, 1);
		anchorCkN = E_Kth_approx(n, anchor_C, sigma, n);
	}
}
