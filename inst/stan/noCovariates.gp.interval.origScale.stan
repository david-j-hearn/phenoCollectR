functions {
        real log_cat_probs(int stage, real t, real mu, real sigma, real D, int debug) {
                if(debug) {
                        print("time, onset, sigma, duration: ", [t, mu, sigma, D]);
                }

                // log P(before) = log P(S > t)
		if(stage==1) {
                	return normal_lccdf(t | mu, sigma);
		}

                // log P(after) = log P(S + D < t) = log P(S < t - D)
		if(stage==3) {
                	return normal_lcdf(t - D | mu, sigma);
		}
                
                // log P(during) = log( Phi((t-mu)/sigma) - Phi((t-D-mu)/sigma) )
                real log_cdf1 = normal_lcdf(t | mu, sigma);
                real log_cdf0 = normal_lcdf(t - D | mu, sigma);
                return log_diff_exp(log_cdf1, log_cdf0);
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

	//scale / translation
	real<lower=0> minResponse;	//minimum possible response (0 for year)
	real maxResponse;	//maximum possible response (365 or 366 for year)

	//observations
	vector<lower=minResponse, upper=maxResponse>[N] t;	//observations
	array[N] int<lower=1, upper=3> stage; //before, during, after stage

	//hyperparameters
	real<lower=minResponse,upper=maxResponse> mean_mean_onset;
	real<lower=minResponse,upper=0.3*maxResponse> sd_mean_onset;
	real<lower=minResponse,upper=maxResponse> mean_mean_duration;
	real<lower=minResponse,upper=0.3*maxResponse> sd_mean_duration;
	real<lower=minResponse,upper=0.3*maxResponse> mean_sd;
	real<lower=minResponse,upper=0.3*maxResponse> sd_sd;

	//diagnostics
	int<lower=0, upper=10> debug;
	int<lower=0, upper=1> drop_nc;
	int<lower=0, upper=1> drop_ll;
}

parameters {
	real<lower=minResponse, upper=maxResponse> mu_O;
	real uc_sigma;
	real uc_mu_D;
}

transformed parameters {
	real<lower=0> sigma = log1p_exp(uc_sigma);
	real<lower=mu_O, upper=maxResponse> mu_C;
	real<lower=0, upper=maxResponse> mu_D = log1p_exp(uc_mu_D);
	mu_C = mu_O + mu_D;
}

model {

	if(debug) {
		print("Mus: ", [mu_O, mu_C], ", sigma:", sigma);
	}

	if(!drop_ll) {
		for(i in 1:N) {
                        target += log_cat_probs(stage[i], t[i], mu_O, sigma, mu_D,debug);
		}
	}

	mu_O ~ normal(mean_mean_onset, sd_mean_onset);
	mu_D ~ normal(mean_mean_duration, sd_mean_duration);
	sigma ~ normal(mean_sd, sd_sd);

}

generated quantities {
	real<lower=minResponse,upper=maxResponse> mu_T;

	real mu_Ok1;
	real mu_CkN;

	mu_T = mu_O + mu_D/2;
		if(processExtremes) {
			mu_Ok1 = E_Kth_approx(n, mu_O, sigma, 1) ;
			mu_CkN = E_Kth_approx(n, mu_C, sigma, n) ;
		}
}
