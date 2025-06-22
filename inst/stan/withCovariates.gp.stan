functions {
	// Implemented in C++ in gp.hpp
	real log_nc(real mu_O, real mu_C, real sigma, int debug);
	real log_nc(real mu_O, real mu_C, real sigma);

	// The likelihood exluding the normalisation constant (log_nc)
	real obs_gp_lpdf(real t, real mu_O, real mu_C, real sigma, int debug) {
		if(t > mu_O) {
			if(debug == 2) {
				print("T: ", t, ", ", [normal_lcdf(mu_C | t, sigma), normal_lcdf(mu_O | t, sigma)]);
			}
			// observed time after onset, CDF is going to be small (and thus precise).
			// straightforward computation
			return(log_diff_exp(normal_lcdf(mu_C | t, sigma), normal_lcdf(mu_O | t, sigma)));
		} else {
			if(debug == 2) {
				print("T - flip: ", t, ", ", [normal_lcdf(2 * t - mu_O | t, sigma), normal_lcdf(2 * t - mu_C | t, sigma),
						log_diff_exp(normal_lcdf(2 * t - mu_O | t, sigma), normal_lcdf(2 * t - mu_C | t, sigma))]);
			}
			// observed time before onset, CDF might be imprecise.
			// flipping around t to get small values again
			return(log_diff_exp(normal_lcdf(2 * t - mu_O | t, sigma), normal_lcdf(2 * t - mu_C | t, sigma)));
		}
	}

	real get_alpha_inv_gamma( real mu_i, real sigma ) {
		real mu = mu_i;
		real epsilon = 1e-6;
		if(mu< sigma) { 
			mu = sigma + epsilon;
		}
		return(2.0 + (mu/sigma)*(mu/sigma));
	}

	real get_beta_inv_gamma( real mu_i, real sigma ) {
		real mu = mu_i;
		real epsilon = 1e-6;
		if(mu< sigma) { 
			mu = sigma + epsilon;
		}
		return((1 + (mu/sigma)*(mu/sigma))*mu);
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
	int<lower=0, upper=1> drop_nc;
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

	real anchorCessationMean;
	real<lower=0> anchorCessationSD;

	//Hyperparameters mean and sd of sigma
	real<lower=0> sigmaMean;
	real<lower=0> sigmaSD;

}

transformed data {

	real epsilon = 1e-6;
	vector[N] mu_D_raw_prior;
	vector[N] invG_alpha_D_prior;
	vector[N] invG_beta_D_prior;
	real alpha_D_raw_prior;

	vector[N] mu_O_raw_prior;
	vector[N] invG_alpha_O_prior;
	vector[N] invG_beta_O_prior;
	real alpha_O_raw_prior;

	real<lower=2> alpha_invGammaO;
	real beta_invGammaO;

	real<lower=2> alpha_invGammaD;
	real beta_invGammaD;

	real<lower=2> alpha_invGammaC;
	real beta_invGammaC;

	real<lower=2> alpha_invGammaS;
	real beta_invGammaS;

	alpha_invGammaO = get_alpha_inv_gamma(anchorOnsetMean , anchorOnsetSD) ;
	beta_invGammaO = get_beta_inv_gamma(anchorOnsetMean , anchorOnsetSD) ;

	//print("Mus: ", [anchorDurationMean*365.0, anchorOnsetMean*365], ", sigma:", sigma*365);
	//print("SDs: ", [anchorDurationSD*365.0, anchorOnsetSD*365], ", sigma:", sigma*365);
	print("Mus: ", [anchorDurationMean*365.0, anchorOnsetMean*365]);
	print("SDs: ", [anchorDurationSD*365.0, anchorOnsetSD*365]);
	print("Mus raw: ", [anchorDurationMean, anchorOnsetMean]);
	print("SDs raw: ", [anchorDurationSD, anchorOnsetSD]);

	alpha_invGammaD = get_alpha_inv_gamma(anchorDurationMean , anchorDurationSD) ;
	beta_invGammaD = get_beta_inv_gamma(anchorDurationMean , anchorDurationSD) ;
	print("alphas: ", [alpha_invGammaO, alpha_invGammaD]);
	print("betas: ", [beta_invGammaO,beta_invGammaD]);

	alpha_invGammaC = get_alpha_inv_gamma(anchorCessationMean , anchorCessationSD) ;
	beta_invGammaC = get_beta_inv_gamma(anchorCessationMean , anchorCessationSD) ;

	alpha_invGammaS = get_alpha_inv_gamma(sigmaMean, sigmaSD) ;
	beta_invGammaS = get_beta_inv_gamma(sigmaMean, sigmaSD) ;

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

	//pre-compute the prior mean, sd, inv. gamma alpha and inv. gamma beta values for each of the covariate values
	alpha_D_raw_prior = anchorDurationMean - dot_product(mean_X_D_raw , betaDurationMeans);
	for(i in 1:N) {
		mu_D_raw_prior[i] = alpha_D_raw_prior + dot_product(X_D_raw[i], betaDurationMeans);
		if(mu_D_raw_prior[i] < anchorDurationSD) {				//needed for inv. gamma, too many rejects outwise
			mu_D_raw_prior[i] = anchorDurationSD + epsilon;
		}
		invG_alpha_D_prior[i] = get_alpha_inv_gamma(mu_D_raw_prior[i], anchorDurationSD);
		invG_beta_D_prior[i] = get_beta_inv_gamma(mu_D_raw_prior[i], anchorDurationSD);
	}
	alpha_O_raw_prior = anchorOnsetMean - dot_product(mean_X_O_raw , betaOnsetMeans);
	for(i in 1:N) {
		mu_O_raw_prior[i] = alpha_O_raw_prior + dot_product(X_O_raw[i], betaOnsetMeans);
		if(mu_O_raw_prior[i] < anchorOnsetSD) {				//needed for inv. gamma, too many rejects outwise
			mu_O_raw_prior[i] = anchorOnsetSD + epsilon;
		}
		invG_alpha_O_prior[i] = get_alpha_inv_gamma(mu_O_raw_prior[i], anchorOnsetSD);
		invG_beta_O_prior[i] = get_beta_inv_gamma(mu_O_raw_prior[i], anchorOnsetSD);
	}
}

parameters {
	//slope parameters, beta, and mean response, (anchor) at mean covariate value, and onset and cessation distribution standard deviation, sigma

	vector[K_O] beta_O_raw;
	real anchor_O_raw;

	vector[K_D] beta_D_raw;
	real<lower=0> anchor_D_raw;

	real<lower=0> sigma_raw;
}



transformed parameters {
	real alpha_O_raw; //intercept of the linear model for Onset
	real alpha_D_raw; //intercept of the linear model for Duration

	vector[N] mu_O_raw;
	vector[N] mu_C_raw;
	vector[N] mu_D_raw;

	alpha_O_raw = anchor_O_raw - dot_product(mean_X_O_raw , beta_O_raw);
	alpha_D_raw = anchor_D_raw - dot_product(mean_X_D_raw , beta_D_raw);

	mu_O_raw =  alpha_O_raw + X_O_raw * beta_O_raw;
	mu_D_raw =  alpha_D_raw + X_D_raw * beta_D_raw;

	mu_C_raw = mu_O_raw + mu_D_raw;
}



model {
	//if(debug) {
	////print("Mus: ", [mu_O, mu_C], ", sigma:", sigma);
	////print("intercept, elev, date, lat: ", [intercept, betaElev, betaDate, betaLat]);
	//}

	//Calculate the likelihood
	if(!drop_ll) {
		for(i in 1:N) {
			target += obs_gp_lpdf(T_raw[i] | mu_O_raw[i], mu_C_raw[i], sigma_raw, debug);
			if(!drop_nc) {
				target += -log_nc(mu_O_raw[i], mu_C_raw[i], sigma_raw, debug); //p(state = 1 | all but t)
			}
		}
	}

	//Define priors

	//Putting priors on the individual mu_C and mu_D reduces divergences rate tremendously!
	//	inverse gamma assures the prior for these values is positive
	//mu_C_raw ~ inv_gamma(alpha_invGammaC,beta_invGammaC);
	//mu_D_raw ~ normal( anchorDurationMean, anchorDurationSD);
	//mu_D_raw ~ inv_gamma(alpha_invGammaD,beta_invGammaD);
	//mu_O_raw ~ inv_gamma(alpha_invGammaO,beta_invGammaO);

	//mu_D_raw =  alpha_D_raw + X_D_raw * beta_D_raw; //DONE ABOVE

	//These hierarchical priors on the individual mean durations and onsets help to regularize, improve HMC diagnostics by reducing divergences, and still maintain accuracy. 
	if(priors >= 2) {
		for(i in 1:N) {
			//target += inv_gamma_lpdf(mu_D_raw[i] | invG_alpha_D_prior[i], invG_beta_D_prior[i]) * 2 / N;
			//target += inv_gamma_lpdf(mu_O_raw[i] | invG_alpha_O_prior[i], invG_beta_O_prior[i]) * 2 / N;
			target += inv_gamma_lpdf(mu_D_raw[i] | invG_alpha_D_prior[i], invG_beta_D_prior[i]);
			target += inv_gamma_lpdf(mu_O_raw[i] | invG_alpha_O_prior[i], invG_beta_O_prior[i]);
		}
	}

	if(priors >= 1) {
		beta_O_raw ~ normal( betaOnsetMeans, betaOnsetSDs);
		beta_D_raw ~ normal( betaDurationMeans, betaDurationSDs);

		anchor_O_raw ~ normal( anchorOnsetMean, anchorOnsetSD);
		anchor_D_raw ~ normal( anchorDurationMean, anchorDurationSD);

		sigma_raw ~ normal( sigmaMean, sigmaSD);
	}

	//anchor_O_raw ~ inv_gamma(alpha_invGammaO, beta_invGammaO);
	//anchor_D_raw ~ inv_gamma(alpha_invGammaD, beta_invGammaD);

	//beta_O_raw ~ normal(betaOnsetMeans, betaOnsetSDs);
	//beta_D_raw ~ normal(betaDurationMeans, betaDurationSDs);

	//sigma_raw ~ inv_gamma(alpha_invGammaS, beta_invGammaS);
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
	//vector[N] mu_Ok1;
	//vector[N] mu_CkN;

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
