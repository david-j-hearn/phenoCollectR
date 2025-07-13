#' Simulate a response variable based on a linear model with one covariate
#'
#' @description Simulate the response variable for a single linear model of a covariate based on a degenerate Gaussian process with deterministic duration (GP). See Hearn et al. (XXXX) for details. 
#'
#' @param n Sample size 
#' @param slopeO Slope of the onset model
#' @param interceptO Intercept of the onset model
#' @param sigma Standard deviation of the onset distribution and of the cessation distribution
#' @param slopeD Slope of the duration model
#' @param interceptD Intercept of the duration model
#' @param minCovariate Minimum possible value of the response
#' @param maxCovariate Maximum possible value of the response
#'
#' @return A list with the values of the covariate (X), onset times (O), phenophase durations (D), cessation times (C), and observed times (Ts)
#' @export
#'
#' @examples
#' #Set the parameters
#' n=1000 #sample size
#' slopeO = 1 #set the slope of the onset model
#' slopeD = 0.5 #set the slope of the duration model
#' interceptO = 100 #set the intercept of the onset model
#' interceptD = 20 #set the intercept of the duration mdoel
#' sigma = 7 #set the standard deviation of the onset distribution and of the cessation distribution
#' minCovariate = -2 #set the minimum value of the covariate
#' maxCovariate = 30 #set the maximum value of the covariate
#' data = simulateCovariate(n=n, slopeO=slopeO, interceptO=interceptO, sigma=sigma, slopeD=slopeD, interceptD=interceptD, minCovariate=minCovariate, maxCovariate=maxCovariate)
#' #plot the simulated observed collection times
#' plot(data$X, data$Ts, xlab="Mean spring temperature", ylab="Day of year", col="purple", main=NULL, pch=16)
#' points(data$X, data$O, col="red", pch=16, cex=0.3)
#' points(data$X, data$C, col="blue", pch=16, cex=0.3)
#' #Plot the line that passes through the mean observed collection times, onset and cessation
#' abline(a = interceptO + interceptD/2, b = slopeO + slopeD/2, col = "purple", lwd = 2)
#' abline(a = interceptO, b = slopeO, col = "red", lwd = 2)
#' abline(a = interceptO + interceptD, b = slopeO + slopeD, col = "blue", lwd = 2)
#' #plot phenophases for each individual in gray
#' segments(x0 = data$X, y0 = data$O, x1 = data$X, y1 = data$C, col = "gray", lwd=0.5)
simulateCovariate = function(n, slopeO, interceptO, sigma, slopeD, interceptD, minCovariate, maxCovariate) {

x = runif(n, minCovariate, maxCovariate)
O = interceptO + slopeO * x + rnorm(n,0,sigma)
D = interceptD + slopeD * x	#no intrinsic variance for duration under GP
C = O + D
Ts = runif(n, O, C) #works because duration is a constant size interval and time period bounds are not enforced.

out = list(
	X = x,
	O = O,
	D = D,
	C = C,
	Ts = Ts
	)

return(out)
}

#' Simulate a response variable based on multiple correlated covariates
#'
#' @description Simulate the values of a response variable based on simulated values of multiple correlated covariates. Assumes a linear model with normally distributed variation.
#'
#' @param n Sample size
#' @param beta Vector of slope coefficients of the covariates
#' @param cov_names Vector of the names of the covariates
#' @param mu Vector of the means of the covariates; must be of same length as beta. (default: all 0 mean)
#' @param response_name The name of the response variable; must be of same length as beta. (default: "Y")
#' @param Sigma The covariance matrix. Must be of dimension length(beta) X length(beta) (default: identity matrix)
#' @param anchor Marginal mean value of the response variable
#' @param noise_sd Standard deviation of the noise in the response variable
#'
#' @return A data frame with columns labeled with the variable names and rows representing simulation replicates.
#' @export
#' @importFrom MASS mvrnorm
#'
#' @examples
#' #Set the model parameters
#' slopes = c(1,2,3)
#' means = c(10,20,30)
#' covariance_names = c("x1", "x2", "x3") 
#' response_name = "y"
#' covariance_matrix = matrix(c( 1.0, 0.5, 0.3, 0.5, 2.0, 0.4, 0.3, 0.4, 1.5), nrow = 3, byrow = TRUE)
#' mean_response = 100
#' noise = 3
#' n=1000
#' #Simulate the data
#' simulated_data = simulateCorrelatedCovariateData(n=n, beta=slopes, cov_names=covariance_names, mu = means, Sigma=covariance_matrix, anchor=mean_response, response_name = response_name, noise_sd = noise)
#' #Make a scatter plot of the simulated data
#' plot(simulated_data$x1, simulated_data$y, main=NULL, xlab="X1", ylab="Y")
simulateCorrelatedCovariateData = function(n, beta, cov_names, mu = NULL, response_name = "Y",
		Sigma = NULL, anchor = 0, noise_sd = 1) {
# Load required package
	if (!requireNamespace("MASS", quietly = TRUE)) {
		stop("Please install the 'MASS' package with install.packages('MASS')")
	}

# Check input
		p = length(beta)
		if (length(cov_names) != p) {
			stop("Length of cov_names must match length of beta.")
		}

# Default covariance matrix (identity)
	if (is.null(Sigma)) {
		Sigma = diag(p)
	} else {
		if (!all(dim(Sigma) == c(p, p))) {
			stop("Sigma must be a square matrix with dimensions matching length of beta.")
		}
	}

# Simulate covariates
	if(is.null(mu)) {
	mu = rep(0, p)
	}
	if(length(mu) != p) {
		stop("The length of the mu vector must match the length of the beta vector.")
	}
		X = mvrnorm(n = n, mu = mu, Sigma = Sigma)
		colnames(X) = cov_names

# Linear predictor
		linear_part = X %*% beta

# Adjust intercept for anchor
		intercept = anchor - mean(linear_part)

# Simulate noise
		epsilon = rnorm(n, mean = 0, sd = noise_sd)

# Simulate response
		Y = as.vector(intercept + linear_part + epsilon)

# Combine into dataframe
		df = data.frame(Y = Y, X)
		names(df)[1] = response_name

		return(df)
}

simulatePopulation.BB = function(N, minResponse=0, maxResponse=365, mu_O, sigma_O, mu_D, sigma_D, mins=1, maxs=3000, betaDuration=NA) {
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

		cnt = 1
		t_start = numeric(N)
		t_end = numeric(N)
		Ts = numeric(N)
		duration_raw = numeric(N)
		duration = numeric(N)
		observed = numeric(N)

		t_start = rbeta(N, alpha_s, beta_s)
		duration_raw = rbeta(N, alpha_d, beta_d)
		duration = duration_raw * (1 - t_start)
		t_end =  t_start + duration
		observed = runif(N) 
		#sampled times of individuals in the phenophase is a subset of all the randomly observed times of the individuals (not all individuals are in the phenophase at a randomly observed time)
		Ts = observed[observed>t_start & observed<t_end]

		#There are more efficient ways to do this
		#while(cnt <= N) {
			##sample the start time
			#t_start[cnt] = rbeta(1, alpha_s, beta_s)
			##sample the unscaled duration
			#duration_raw[cnt] = rbeta(1, alpha_d, beta_d)
			##scale the duration
			#duration[cnt] = duration_raw[cnt] * (1 - t_start[cnt])
			##calculate the end time
			#t_end[cnt] = t_start[cnt] + duration[cnt]
			##sample the observed time
			#observed = runif(1) #CORRECT, but cessation time distribution is biased
			##check if observed time is in phenophase
			#if(observed>t_start[cnt] && observed<t_end[cnt]) { #CORRECT
				#Ts[cnt] = observed #CORRECT
				##Ts[cnt] = runif(1, t_start[cnt], t_end[cnt]) #BIASED for BB model
				#cnt = cnt+1
			#}
		#}

		t_start = minResponse + t_start * (maxResponse-minResponse)
		t_end = minResponse + t_end * (maxResponse-minResponse)
		Ts = minResponse + Ts * (maxResponse-minResponse)
		duration = duration * (maxResponse-minResponse)
		duration_raw = duration_raw * (maxResponse-minResponse)

		#MC estimate
		cessation_sd = sd(t_end)
		duration_sd = sd(duration)

		return(list(
					error = F,
					error_m = "No errors detected during simulation under beta onset, beta duration model.",
					N = N,
					minResponse = minResponse,
					maxResponse = maxResponse,
					O = t_start,
					Ts = Ts,
					C = t_end,
					D = duration,
					D_raw = duration_raw,
					CkN = max(t_end),
					Ok1 = min(t_start),
					R = max(t_end) - min(t_start),
					mu_O = mu_O,
					sigma_O = sigma_O,
					mu_D_raw = mu_D, #for raw duration
					sigma_D_raw = sigma_D, #for raw duration
					mu_D = duration_mean * (1 - onset_mean) * (maxResponse - minResponse), #infers the beta parameters for the raw distribution, here is the scaled result for the duration
					sigma_D = duration_sd,
					mu_C = minResponse + (onset_mean + duration_mean * ( 1 - onset_mean)) * (maxResponse - minResponse),
					sigma_C = cessation_sd,
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
		#if(mu_O - 2*sigma < minResponse) {
			#warning(paste("The onset is close to the beginning of the range of possible times. Multiple samples are likely to be rejected to maintain all observations within the range between ", minResponse, " and ", maxResponse, ". Your results are therefore likely going to be skewed."))
		#}
	#if(mu_C + 2*sigma > maxResponse) {
		#warning(paste("The cessation is close to the upper end of the range of possible times. Multiple samples are likely to be rejected to maintain all observations within the range between ", minResponse, " and ", maxResponse, ". Your results are therefore likely going to be skewed."))
	#}
	if(mu_O>=mu_C) {
		stop("The mean onset must be before the mean cessation. Quitting.")
	}
	d = mu_C - mu_O
		start = rnorm(n,mu_O,sd)
		#start = start[start>0]
		#start = start[start+d<1]
		#tot = length(start)
		#nRej = 0
		#while(tot<n) {
			#nRej = nRej + (n-tot)
				#temp = rnorm(n-tot,mu_O,sd)
				#temp = temp[temp>0]
				#temp = temp[temp+d<1]
				#start = c(start,temp)
				#tot = tot + length(temp)
		#}
	##if(nRej>0.1*n) {
		#warning(paste("Warning: rejected start times below 0 ", nRej, " times. This can bias inferences. If ", nRej, "is less than, say, 0.1% of the sample size (", (n*0.01), " in your case), this should be ok, otherwise, the model inference will be highly biased"))
			#return(list(error=T))
	#}
	end = start + d
	obs = runif(n,start,end)


	k = (maxResponse-minResponse)

				Ts = minResponse + obs*k
				O = minResponse + start * k
				C = minResponse + end * k
				CkN = minResponse + max(end) * k
				Ok1 = minResponse + min(start) * k

	Ok1 = Ok1%%k
	CkN = CkN%%k
	C = C%%k
	O = O%%k
	Ts = Ts%%k

	C[C < 0] = C[C < 0] + k
	O[O < 0] = O[O < 0] + k
	Ts[Ts < 0] = Ts[Ts < 0] + k
	Ok1[Ok1 < 0] = Ok1[Ok1 < 0] + k
	CkN[CkN < 0] = CkN[CkN < 0] + k

		output = list(
				error = F,
				error_m = "No errors detected during simulation under GP model.",
				N = n,
				minResponse = minResponse,
				maxResponse = maxResponse,
				Ts = Ts,
				O = O,
				C = C,
				Ok1 = Ok1,
				CkN = CkN,
				mu_O = minResponse + mu_O * (maxResponse - minResponse),
				sigma = sigma,
				mu_D = d * (maxResponse - minResponse),
				mu_C = minResponse + mu_C * (maxResponse - minResponse)
			     )
		return(output)
}

#' Simulate phenological states for individuals of a population
#'
#' @description Simulate the phenophase onset, duration, and cessation times for individuals of a population. Phenological extremes (first onset, last cessation) are also provided. 
#'
#' Data can be simulated under either the beta onset, beta duration model (BB) or under the Gaussian process model (GP). 
#' 
#' Note that for the BB model, the number of observed times will be less than the number of individuals in the population because individuals are viewed at random times, and not all individuals are in the phenophase at the randomly observed time. Individuals with longer phenophases are more likely to be recorded. Only simulated specimens in the phenophase are recorded. This is not the case for the GP model, because under the GP model, every individual is equally likely to be recorded, so the observed time of each individual can be randomly sampled between the onset time and the cessation time. 
#' 
#' Additionally, durations are scaled under the BB model so that all events are guaranteed to occur in the time period, whereas durations are not scaled under the GP model, and times that fall outside of the time period are "wrapped" around into the next time period.
#'
#' @param N Population size
#' @param mu_O Mean of the onset times of individuals in the population
#' @param sigma_O Standard deviation of the onset times of individuals in the population
#' @param mu_D_raw Unscaled mean duration. This is scaled under the BB model to assure onset times and cessation times fall between the minimum and maximum response times
#' @param sigma_D Unscaled standard deviation of the durations
#' @param minResponse The minimum possible response time (default: 0)
#' @param maxResponse The maximum possible response time (default: 365)
#' @param mins The minimum allowable value for the beta distribution shape parameters. Must be positive and less than maxs. (default: 1.5)
#' @param maxs The maximum allowable value for the beta distribution shape parameters. Must be positive and greater than mins. (default: 3000)
#' @param type The model type: BB or GP (default: GP)
#'
#' @return A list with the input parameter values, the simulated onsets (O), durations (D), cessations (C), observed collection times (Ts). Specifics of the output depend on whether the BB or the GP model was applied.
#' @export
#'
#' @examples
#' #Set the parameters
#' N = 1000000
#' mu_O = 30
#' sigma_O = 10
#' mu_D_raw = 40
#' sigma_D = 14
#' #Simulate the data under the beta onset, beta duration (BB) model
#' data = simulatePopulation(N=N, mu_O=mu_O, sigma_O=sigma_O, mu_D_raw=mu_D_raw, sigma_D=sigma_D, type="BB")
#' #Plot histograms of the phenological values
#' xlim = c(min(data$O, data$C), max(data$O, data$C))
#' breaks = seq(xlim[1],xlim[2], length.out=100)
#' hist(data$O, col=rgb(1,0,0,0.3), xlab="Day of year", probability=TRUE, breaks=breaks, xlim=xlim, main=NULL) #Onset
#' hist(data$Ts, col=rgb(1,0,1,0.3), probability=TRUE, breaks=breaks, add=TRUE) #Observed collection times
#' hist(data$C, col=rgb(0,0,1,0.3), probability=TRUE, breaks=breaks, add=TRUE) #Cessation
#' abline(v=data$Ok1, col="yellow") #First onset for the population as yellow vertical line
#' abline(v=data$CkN, col="cyan") #Last cessation for the population as cyan vertical line
#' curve(dO(x,mu_O=mu_O, sigma_O=sigma_O, type="BB"), add=TRUE, col="red") #overlay theoretical density curve for onset
#' curve(dT(x,mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D_raw, sigma_D=sigma_D, type="BB"), add=TRUE, col="purple") #overlay theoretical density curve for observed times
#' curve(dC(x,mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D_raw, sigma_D=sigma_D, type="BB"), add=TRUE, col="blue") #overlay theoretical density curve for cessation
#' #Use the above parameter values, but simulate data under the GP model (default)
#' #!!Note the "wrapping around" effect for which values below 0 get wrapped to the end of the "previous" time period 
#' data = simulatePopulation(N=N, mu_O=mu_O, sigma_O=sigma_O, mu_D_raw=mu_D_raw)
#' #Plot histograms of the phenological values
#' dev.new()
#' xlim = c(min(data$O, data$C), max(data$O, data$C))
#' breaks = seq(xlim[1],xlim[2], length.out=100)
#' hist(data$O, col=rgb(1,0,0,0.3), xlab="Day of year", probability=TRUE, breaks=breaks, xlim=xlim, main=NULL) #Onset
#' hist(data$Ts, col=rgb(1,0,1,0.3), probability=TRUE, breaks=breaks, add=TRUE) #Observed collection times
#' hist(data$C, col=rgb(0,0,1,0.3), probability=TRUE, breaks=breaks, add=TRUE) #Cessation
#' abline(v=data$Ok1, col="yellow") #First onset for the population as yellow vertical line
#' abline(v=data$CkN, col="cyan") #Last cessation for the population as cyan vertical line
#' curve(dO(x,mu_O=mu_O, sigma_O=sigma_O), add=TRUE, col="red") #overlay theoretical density curve for onset
#' curve(dT(x,mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D_raw, sigma_D=sigma_D), add=TRUE, col="purple") #overlay theoretical density curve for observed times
#' curve(dC(x,mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D_raw, sigma_D=sigma_D), add=TRUE, col="blue") #overlay theoretical density curve for cessation
simulatePopulation =  function(N, mu_O, sigma_O, mu_D_raw, sigma_D=NA, minResponse=0, maxResponse=365, mins=1.5, maxs=3000, type="GP") {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D_raw, sigma_D=sigma_D, N=NA, n=N, minResponse=minResponse, maxResponse=maxResponse) # some redundancy with below...
	if(N <= 0) {
		stop("Population size must be positive during simulations.")
	}
	if(sigma_O <=0 ) {
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
    if(sigma_D <=0 ) {
        stop("Standard deviations (sigmas) must be positive.")
    }
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
sample_conditional_covariates = function(x_target, x_column = 1, covars, copula_fit, n_samples = 1000) {
  Sigma = getSigma(copula_fit@copula)
  K = ncol(Sigma)
  all_names = colnames(covars)

  # Empirical CDF of the target covariate
  ecdf_target = ecdf(covars[[x_column]])
  u_x = ecdf_target(x_target)
  z_x = qnorm(u_x)

  # Partition correlation matrix
  idx_y = setdiff(1:K, x_column)
  Sigma_xx = Sigma[x_column, x_column]
  Sigma_yy = Sigma[idx_y, idx_y]
  Sigma_yx = Sigma[idx_y, x_column]

  # Conditional normal parameters
  mu_y = Sigma_yx / Sigma_xx * z_x
  Sigma_cond = Sigma_yy - tcrossprod(Sigma_yx) / Sigma_xx

  # Simulate latent Z_{2:K}
  z_y = mvrnorm(n_samples, mu = mu_y, Sigma = Sigma_cond)

  # Inverse probit and quantile transform
  u_y = pnorm(z_y)
  x_y = matrix(NA, nrow = n_samples, ncol = length(idx_y))
  for (j in seq_along(idx_y)) {
    x_y[, j] = quantile(covars[[idx_y[j]]], u_y[, j], type = 8)
  }

  df = as.data.frame(x_y)
  colnames(df) = all_names[idx_y]
  df[[all_names[x_column]]] = x_target
  df = df[all_names]  # order columns
  return(df)
}
