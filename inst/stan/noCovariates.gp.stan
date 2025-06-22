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

	int<lower = 0> n;	//population size

	int<lower=0, upper=1> processExtremes;	//process first onset and last cessation using an asymptotic approximation (numerical integration fails)

	//observations
	vector<lower=0, upper=1>[N] t;	//observations

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
	real<lower=0> mu_D_raw;
	real<lower=0> sigma_raw;
}

transformed parameters {
	real<lower=mu_O_raw, upper=1> mu_C_raw;
	mu_C_raw = mu_O_raw + mu_D_raw;
}

model {
	if(debug) {
		print("Mus: ", [mu_O_raw, mu_C_raw], ", sigma_raw:", sigma_raw);
	}
	if(!drop_nc) {
		target += -N * log_nc(mu_O_raw, mu_C_raw, sigma_raw, debug); //p(state = 1 | all but t)
	}
	if(!drop_ll) {
		for(i in 1:N) {
			target += obs_gp_lpdf(t[i] | mu_O_raw, mu_C_raw, sigma_raw, debug);
		}
	}
	mu_O_raw ~ normal(mean_mean_onset, sd_mean_onset);
	mu_D_raw ~ normal(mean_mean_duration, sd_mean_duration);
	sigma_raw ~ normal(mean_sd, sd_sd);
}

generated quantities {
	real<lower=minResponse,upper=maxResponse> mu_O;
	real<lower=0> mu_D;
	real<lower=minResponse,upper=maxResponse> mu_C;
	real<lower=minResponse,upper=maxResponse> mu_T;
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
