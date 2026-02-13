functions {

	real Phi_stable(real z) {
		real p;
		if (z < 0)
			p = 1 - exp(normal_lccdf(z | 0, 1));
		else
			p = exp(normal_lcdf(z | 0, 1));
		return p;
	}

	real H_norm_cdf_int(real t, real mu, real sigma) {
		// H(t) = (t-mu)*Phi((t-mu)/sigma) + sigma*phi((t-mu)/sigma)
		real z = (t - mu) / sigma;
		real Phi_z = Phi_stable(z);
		real phi_z = exp(normal_lpdf(z | 0, 1));
		return (t - mu) * Phi_z + sigma * phi_z;
		//return (t - mu) * Phi(z) + sigma * exp(normal_lpdf(z | 0, 1));
	}

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

	vector log_p_y_marginal(real mu, real sigma, real D) {

		vector[3] lp;

		// A = ∫_0^1 Phi((t-mu)/sigma) dt
		real A = H_norm_cdf_int(1, mu, sigma) - H_norm_cdf_int(0, mu, sigma);

		// B = ∫_0^1 Phi((t-(mu+D))/sigma) dt
		real B = H_norm_cdf_int(1, mu + D, sigma) - H_norm_cdf_int(0, mu + D, sigma);

		real p_before = 1 - A;
		real p_after  = B;
		real p_during = A - B;

		// clamp for numerical safety
		// numerical safety (also avoids log(0))
		p_before = fmin(1 - 1e-12, fmax(1e-12, p_before));
		p_during = fmin(1 - 1e-12, fmax(1e-12, p_during));
		p_after  = fmin(1 - 1e-12, fmax(1e-12, p_after));

		lp[1] = log(p_before);
		lp[2] = log(p_during);
		lp[3] = log(p_after);

		return lp;

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
	//sample size
	int<lower=0> N;	//sample size

	int<lower=0> n;	//population size

	int<lower=0, upper=1> processExtremes;	//process first onset and last cessation using an asymptotic approximation (numerical integration fails)

	//observations
	vector<lower=0, upper=1>[N] t;	//observations
	array[N] int<lower=1, upper=3> stage; //before, during, after stage

	//scale / translation
	real<lower=0> minResponse;	//minimum possible response (0 for year)
	real maxResponse;	//maximum possible response (365 or 366 for year)

	//hyperparameters
	real<lower=0,upper=1> mean_mean_onset;
	real<lower=0,upper=0.3> sd_mean_onset;
	real<lower=0,upper=1> mean_mean_duration;
	real<lower=0,upper=0.3> sd_mean_duration;
	real<lower=0,upper=0.3> mean_sd;
	real<lower=0,upper=0.3> sd_sd;

	//diagnostics
	int<lower=0, upper=10> debug;
	int<lower=0, upper=1> drop_nc;
	int<lower=0, upper=1> drop_ll;
}

parameters {
	real<lower=0, upper=1> mu_O_raw;
	real uc_sigma_raw;
	real uc_mu_D_raw;
}

transformed parameters {
	real mu_C_raw;
	real<lower=0> mu_D_raw = log1p_exp(uc_mu_D_raw);
	real<lower=0>  sigma_raw = log1p_exp(uc_sigma_raw);
	mu_C_raw = mu_O_raw + mu_D_raw;
}

model {
	if(debug) {
		print("Mus: ", [mu_O_raw, mu_C_raw], ", sigma_raw:", sigma_raw);
	}

	if(!drop_ll) {
		//This needs to be in the for loop when there are covariates when not dropping the normalization constant (nc)
		//Normalization not needed - just probability in the stage...
		//vector[3] log_py_marg = log_p_y_marginal(mu_O_raw, sigma_raw, mu_D_raw);

		for(i in 1:N) {
			target += log_p_y_given_t(stage[i], t[i], mu_O_raw, sigma_raw, mu_D_raw,debug);

		//Normalization not needed - just probability in the stage...
			//if(!drop_nc) {
			//target +=  -log_py_marg[stage[i]];
			//}
		}
	}
//Works with VERY weak, biased priors
	mu_O_raw ~ normal(mean_mean_onset, sd_mean_onset);
	mu_D_raw ~ normal(mean_mean_duration, sd_mean_duration);
	sigma_raw ~ normal(mean_sd, sd_sd);
}

generated quantities {
	real mu_O;
	real mu_D;
	real mu_C;
	real mu_T;
	real<lower=0> sigma;

	real mu_Ok1;
	real mu_CkN;

	mu_O = minResponse + mu_O_raw * (maxResponse - minResponse);
	mu_D = mu_D_raw * (maxResponse - minResponse);
	mu_C = minResponse + mu_C_raw * (maxResponse - minResponse);
	sigma = sigma_raw * (maxResponse - minResponse);
	mu_T = mu_O + mu_D/2;
		if(processExtremes) {
			mu_Ok1 = E_Kth_approx(n, mu_O, sigma, 1) ;
			mu_CkN = E_Kth_approx(n, mu_C, sigma, n) ;
		}
}
