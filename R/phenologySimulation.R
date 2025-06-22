#' @export
simulateCovariate = function(n, slopeO, interceptO, sigma, slopeD, interceptD, minCovariate, maxCovariate) {

x = runif(n, minCovariate, maxCovariate)
O = interceptO + slopeO * x + rnorm(n,0,sigma)
D = interceptD + slopeD * x	#no intrinsic variance for duration under GP
C = O + D
Ts = runif(n, O, C)

out = list(
	X = x,
	O = O,
	D = D,
	C = C,
	Ts = Ts
	)

return(out)

}

#' @export
#' @importFrom MASS mvrnorm
simulateCorrelatedCovariateData = function(n, beta, cov_names, response_name = "Y",
		Sigma = NULL, anchor = 0, noise_sd = 1) {
# Load required package
	if (!requireNamespace("MASS", quietly = TRUE)) {
		stop("Please install the 'MASS' package with install.packages('MASS')")
	}

# Check input
		p  =  length(beta)
		if (length(cov_names) != p) {
			stop("Length of cov_names must match length of beta.")
		}

# Default covariance matrix (identity)
	if (is.null(Sigma)) {
		Sigma  =  diag(p)
	} else {
		if (!all(dim(Sigma) == c(p, p))) {
			stop("Sigma must be a square matrix with dimensions matching length of beta.")
		}
	}

# Simulate covariates
	mu  =  rep(0, p)
		X  =  mvrnorm(n = n, mu = mu, Sigma = Sigma)
		colnames(X)  =  cov_names

# Linear predictor
		linear_part  =  X %*% beta

# Adjust intercept for anchor
		intercept  =  anchor - mean(linear_part)

# Simulate noise
		epsilon  =  rnorm(n, mean = 0, sd = noise_sd)

# Simulate response
		Y  =  as.vector(intercept + linear_part + epsilon)

# Combine into dataframe
		df  =  data.frame(Y = Y, X)
		names(df)[1]  =  response_name

		return(df)
}

simulatePopulation.BB  =  function(N, minResponse=0, maxResponse=365, mu_O, sigma_O, mu_D, sigma_D, mins=1, maxs=3000, betaDuration=NA) {

#maintain original copies
	onset_mean = mu_O
		duration_mean = mu_D
		onset_sd = sigma_O
		duration_sd = sigma_D

		if(is.na(betaDuration) || !betaDuration) {
			stop("The duration mean and SD should be provided in terms of the beta distribution that samples the raw durations. This will be scaled by (1 - mu_O) after mu_O has been scaled between 0 and 1 to assure that duration + onset falls between minResponse and maxResponse times. Take your mu_D that you want to simulate and divide it by (1 - (mu_O-minResponse)/(maxResponse-minResponse)) to get the duration to input here. Note that sigma_D will also be scaled. Please input betaDistribution=T when you call this function to acknowledge understanding of this convention.")
		}

#scale to beta distribution support
	onset_mean = (onset_mean-minResponse)/(maxResponse-minResponse)
		onset_sd = onset_sd / (maxResponse - minResponse)
		duration_mean = duration_mean / (maxResponse - minResponse)
		duration_sd = duration_sd / (maxResponse - minResponse)

# Convert mean/sd to alpha/beta for onset
		alpha_s = beta_alpha(onset_mean, onset_sd)
		beta_s = beta_beta(onset_mean, onset_sd)

# Convert mean/sd to alpha/beta for duration
		alpha_d = beta_alpha(duration_mean, duration_sd)
		beta_d = beta_beta(duration_mean, duration_sd)


		if(alpha_s<mins || alpha_s>maxs || alpha_d<mins || alpha_d>maxs || beta_s<mins || beta_s>maxs || beta_d<mins || beta_d>maxs) {
			print(c(alpha_s,beta_s,alpha_d,beta_d))
				return(list(error_m = "Infeasible parameter value provided as input under beta onset, beta duration model. Try different parameter values, paying close attention to the scale of sigma. Sigmas that are too large cannot be accommodated by the beta distribution.", error=T))

		}

# Simulate onset times
	t_start  =  rbeta(N, alpha_s, beta_s)

# Simulate durations scaled to (0, 1 - onset)
		raw_duration  =  rbeta(N, alpha_d, beta_d)
		duration  =  raw_duration * (1 - t_start)

		t_end  =  t_start + duration

		t_start = minResponse + t_start * (maxResponse-minResponse)
		t_end = minResponse + t_end * (maxResponse-minResponse)

# Uniform sample between start and end for each individual
		#observed  =  runif(N, minResponse = t_start, maxResponse = t_end)
        observed  =  runif(N, min = t_start, max = t_end) # fixed

		cessation_sd = sd(t_end)

		return(list(
					error = F,
					error_m = "No errors detected during simulation under beta onset, beta duration model.",
					N = N,
					minResponse = minResponse,
					maxResponse = maxResponse,
					O = t_start,
					Ts = observed,
					C = t_end,
					mu_O = mu_O,
					sigma_O = sigma_O,
					mu_D_raw = mu_D, #for raw duration
					sigma_D_raw = sigma_D, #for raw duration
					mu_D = duration_mean * (1 - onset_mean) * (maxResponse - minResponse), #infers the beta parameters for the raw distribution, here is the scaled result for the duration
					mu_C = minResponse + (onset_mean + duration_mean * ( 1 - onset_mean)) * (maxResponse - minResponse),
					sigma_C = cessation_sd,
					CkN = max(t_end),
					Ok1 = min(t_start),
					alpha_s = alpha_s,
					beta_s = beta_s,
					alpha_d = alpha_d,
					beta_d = beta_d
					) )
}

simulatePopulation.GP = function(N, mu_O, mu_C, sigma, minResponse=0, maxResponse=365) {
	n = N
		sd = sigma
		mu_O = (mu_O-minResponse) / (maxResponse - minResponse)
		mu_C = (mu_C-minResponse) / (maxResponse - minResponse)
		sd = sd / (maxResponse - minResponse)
		if(mu_O - 2*sigma < minResponse) {
			warning(paste("The onset is close to the beginning of the range of possible times. Multiple samples are likely to be rejected to maintain all observations within the range between ", minResponse, " and ", maxResponse, ". Your results are therefore likely going to be skewed."))
		}
	if(mu_C + 2*sigma > maxResponse) {
		warning(paste("The cessation is close to the upper end of the range of possible times. Multiple samples are likely to be rejected to maintain all observations within the range between ", minResponse, " and ", maxResponse, ". Your results are therefore likely going to be skewed."))
	}
	if(mu_O>=mu_C) {
		stop("The mean onset must be before the mean cessation. Quitting.")
	}
	d = mu_C - mu_O
		start = rnorm(n,mu_O,sd)
		start = start[start>0]
		start = start[start+d<1]
		tot = length(start)
		nRej = 0
		while(tot<n) {
			nRej = nRej + (n-tot)
				temp = rnorm(n-tot,mu_O,sd)
				temp = temp[temp>0]
				temp = temp[temp+d<1]
				start = c(start,temp)
				tot = tot + length(temp)
		}
	if(nRej>0.1*n) {
		warning(paste("Warning: rejected start times below 0 ", nRej, " times. This can bias inferences. If ", nRej, "is less than, say, 0.1% of the sample size (", (n*0.01), " in your case), this should be ok, otherwise, the model inference will be highly biased"))
			return(list(error=T))
	}
	end = start + d
		obs = runif(n,start,end)

		output = list(
				error = F,
				error_m = "No errors detected during simulation under GP model.",
				N = n,
				minResponse = minResponse,
				maxResponse = maxResponse,
				Ts = minResponse + obs*(maxResponse-minResponse),
				O = minResponse + start * (maxResponse-minResponse),
				C = minResponse + end * (maxResponse-minResponse),
				mu_O = minResponse + mu_O * (maxResponse - minResponse),
				sigma = sigma,
				mu_D = d * (maxResponse - minResponse),
				mu_C = minResponse + mu_C * (maxResponse - minResponse),
				CkN = minResponse + max(end) * (maxResponse - minResponse),
				Ok1 = minResponse + min(start) * (maxResponse - minResponse),
				rejected = nRej
			     )
		return(output)
}

#' @export
simulatePopulation =  function(N, mu_O, sigma_O, mu_D_raw, sigma_D=NA, minResponse=0, maxResponse=365, mins=1.5, maxs=3000, type="GP") {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D_raw, sigma_D=sigma_D, N=NA, n=N, minResponse=minResponse, maxResponse=maxResponse) # some redundancy with below...
	if(N <= 0) {
		stop("Population size must be positive during simulations.")
	}
	if(sigma_O <=0 || sigma_D <=0 ) {
		stop("Standard deviations (sigmas) must be positive.")
	}
	if(mu_O<=minResponse || mu_O>=maxResponse) {
		stop("Mean onset must be in the interval (minResponse, maxResponse).");
	}
	if(mu_D_raw <= 0 || mu_D_raw >= maxResponse) {
		stop("Raw mean duration must be in the interval (0, maxResponse).");
	}
	if(mins>=maxs) {
		stop("The smallest beta distribution shape parameter must be less than the largest beta distribution shape parameter.")
	}
	if(mins<0 || maxs<0) {
		stop("Beta shape parameters must have positive values.")
	}
	if(type == "BB") {
		warning("The input mean duration, ", mu_D_raw, ", will be scaled so that all individual phenophases fit between ", minResponse, " and ", maxResponse, ".")
			if(is.na(sigma_D)) {
				stop("Missing the standard deviation of the phenophase duration distribution.")
			}
		simulatePopulation.BB(N=N, minResponse=minResponse, maxResponse=maxResponse, mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D_raw, sigma_D=sigma_D, mins=mins, maxs=maxs, betaDuration=T)
	}
	else if(type == "GP") {
		mu_C = mu_O + mu_D_raw
			simulatePopulation.GP(N=N, mu_O=mu_O, mu_C=mu_C, sigma=sigma_O, minResponse=minResponse, maxResponse=maxResponse)
	}
	else {
		stop(paste("Unsupported type: ", type))
	}
}

#' @importFrom copula getSigma
#' @importFrom MASS mvrnorm
sample_conditional_covariates  =  function(x_target, x_column = 1, covars, copula_fit, n_samples = 1000) {
  Sigma  =  getSigma(copula_fit@copula)
  K  =  ncol(Sigma)
  all_names  =  colnames(covars)

  # Empirical CDF of the target covariate
  ecdf_target  =  ecdf(covars[[x_column]])
  u_x  =  ecdf_target(x_target)
  z_x  =  qnorm(u_x)

  # Partition correlation matrix
  idx_y  =  setdiff(1:K, x_column)
  Sigma_xx  =  Sigma[x_column, x_column]
  Sigma_yy  =  Sigma[idx_y, idx_y]
  Sigma_yx  =  Sigma[idx_y, x_column]

  # Conditional normal parameters
  mu_y  =  Sigma_yx / Sigma_xx * z_x
  Sigma_cond  =  Sigma_yy - tcrossprod(Sigma_yx) / Sigma_xx

  # Simulate latent Z_{2:K}
  z_y  =  mvrnorm(n_samples, mu = mu_y, Sigma = Sigma_cond)

  # Inverse probit and quantile transform
  u_y  =  pnorm(z_y)
  x_y  =  matrix(NA, nrow = n_samples, ncol = length(idx_y))
  for (j in seq_along(idx_y)) {
    x_y[, j]  =  quantile(covars[[idx_y[j]]], u_y[, j], type = 8)
  }

  df  =  as.data.frame(x_y)
  colnames(df)  =  all_names[idx_y]
  df[[all_names[x_column]]]  =  x_target
  df  =  df[all_names]  # order columns
  return(df)
}
