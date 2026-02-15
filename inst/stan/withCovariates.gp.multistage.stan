functions {
//Log probabilities for stage at time t
//The times, mus, sigmas in min max transformed scale
//1 <= stage <= S
	real log_p_y_given_t(int S, int stage, real t, vector mu, vector sigma, int debug) {

		if(debug) {
			print("t=", t,
			" stage=", stage,
			" mu_stage=", mu[stage],
			" sigma_stage=", sigma[stage]);
		}

	//Three cases

	//	Case I: stage 1, spans the "start" boundary 
	//		P(t>m, O_2>t) or P(t<M, t>=O_1) = 1 * (1 - CDF(O_2)) + 1 * CDF(O_1)
	//		P(A U B) = 1 - P(A' ∩ B') = 1 - P(O_1 > t, O_2 <= t)
		if(stage == 1) {
			real log_p_comp = normal_lccdf(t | mu[1], sigma[1]) + normal_lcdf(t | mu[2], sigma[2]);   // log P(O1>t, O2<=t)
			return log1m_exp(log_p_comp);         // log(1 - exp(log_p_comp))
		}

	//	Case II: internal stage, stage > 1 and stage < S
	//		P(t>=O_{stage}, t<O_{stage+1})
		if(stage > 1 && stage < S) {
			return normal_lcdf(t | mu[stage], sigma[stage]) + normal_lccdf(t | mu[stage+1], sigma[stage+1]); 
		}

	//	Case III: stage S, last stage
	//		P(t>=O_{S},t<O_1)
		if(stage == S) {
			return normal_lcdf(t | mu[S], sigma[S]) + normal_lccdf(t | mu[1], sigma[1]);   
		}

	reject("Invalid stage: ", stage, " for S=", S);
	return negative_infinity();
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

	//Checks
	int<lower=0, upper=10> debug;
	int<lower=0, upper=1> drop_ll;	//useful if testing prior predictives

	//Handle Extremes
	int<lower=0, upper=1> process_extremes;	//process first onset / last cessation
	int<lower=0> n;	//population size to process extremes

	//Sample size
	int<lower=0> N; //number of observations

	//Response data
	vector[N] T_raw;  //observed times in the min max transformed scale (between 0 and 1)
	real T_min; //minimum time on the original scale
	real T_max; //maximum time on the original scale

	//Stage data
	int<lower=2> S; //number of stages
	array[N] int<lower=1, upper=S> stage;

	//Covariate data
	int<lower=1> C; 		//number of covariates
	matrix[N,C] X_raw;		//covariate data in the min max transformed scale (between 0 and 1)
	vector[C] mean_X_raw;		//covariate means in the min max transformed scale (bewteen 0 and 1)
	vector[C] min_X;            	//min of onset covariates in original scale
	vector<lower=min_X>[C] max_X;	//max of onset covariates in original scale

	int<lower=0> priors; //0 = flat prior, 1 or above = normal prior

	//Hyperparameters: means and SDs for each parameter in the min max transformed scale

	//	Slope coefficients
	matrix[S,C] betaMean;
	matrix<lower=0>[S,C] betaSD;

	//	Anchors (intercepts) 
	vector[S] anchorMean;
	vector<lower=0>[S] anchorSD;

	//	Standard deviations
	vector<lower=0>[S] sigmaMean;
	vector<lower=0>[S] sigmaSD;
}

transformed data {

	real<lower=0> T_range = T_max - T_min;	//range of responses
	vector<lower=0>[C] range_X; 		//covariate ranges in the original scale
	vector[C] mean_X; 			//mean covariate values in the original scale

	//Transform covariate means and ranges to original scale
	for (k in 1:C) {
		range_X[k] = max_X[k] - min_X[k];
		mean_X[k] = mean_X_raw[k] * range_X[k] + min_X[k];
	}
}

parameters {

	//Slope parameters, beta, and anchor at mean covariate value
	matrix[S,C] beta_raw;
	vector<lower=0,upper=1>[S] anchor_raw;
	vector<lower=0.001>[S] sigma_raw;
}

transformed parameters {

	matrix[S,N] O;

	for(j in 1:N) {
		// Stage 2 onset: linear model
		O[j,2] = anchor2_raw[1] + dot_product(row(X_raw,j), beta2_raw);
	
		// Forward stages (3..S)
		for (i in 3:S) {
			real d_linear = alpha_d_raw[i-1] + dot_product(row(X_raw,j), row(beta_d_raw,i-1));
			O[j,i] = O[j,i-1] + softplus(d_linear);
		}

		
		// Wraparound Stage 1: last stage + increment
		real d_linear_S = alpha_d_raw[S] + dot_product(row(X_raw,j), row(beta_d_raw,S));
		O[j,1] = O[j,S] + softplus(d_linear_S);
	}


}

model {

	if(debug) {
	}

	//Calculate the log likelihood
	if(!drop_ll) {
		for(k in 1:S) {
			for(i in 1:N) {
				target += log_p_y_given_t(S, stage[i], T_raw[i], mu_raw[k,i], sigma_raw[k], debug);
			}
		}
	}

	//Set informative priors, if requested
	if(priors >= 1) {
		to_vector(beta_raw) ~ normal(to_vector(betaMean), to_vector(betaSD))
		anchor_raw ~ normal(anchorMean, anchorSD);
		sigma_raw ~ normal(sigmaMeans sigmaSD);
	}
}

generated quantities {

	vector[S] anchor;
	vector[S] anchork1;		
	vector[S] anchorkN;

	vector[S] alpha;
	vector[S] alphak1;
	vector[S] alphakN;

	vector[S,C] beta;
	vector[S,C] betak1;
	vector[S,C] betakN;

	vector<lower=0>[S] sigma;

	// Transform parameters to original scale
	for (k in 1:S) {
		for(i in 1:C) {

			beta[k,i] = (T_range / range_X[i]) * beta_raw[k,i];

			//Calculate slopes for models of extremes (parallel to onset models)
			if(process_extremes) {
				betak1[k,i] = beta[k];
				betakN[k,i] = beta[k];
			}
		}


		// Transform anchor to original scale
		anchor[k] = T_min + T_range * anchor_raw[k];

		// Compute intercept at original scale
		alpha[k] = anchor[k] - dot_product(row(beta,k), mean_X);

		// Transform sigma to original scale
		sigma[k] = T_range * sigma_raw[k];

		//Calculate anchors and intercepts for extremes at original scale
		if(process_extremes) {
			alphak1[k] = alpha[k] + E_Kth_shift( n, sigma[k], 1);	//intercept on first 
			alphakN[k] = alpha[k] + E_Kth_shift( n, sigma[k], n);	//intercept on last 
			anchork1[k] = E_Kth_approx(n, anchor[k], sigma[k], 1);	//anchor on first
			anchorkN[k] = E_Kth_approx(n, anchor[k], sigma[k], n);	//anchor on last
		}
	}
}
