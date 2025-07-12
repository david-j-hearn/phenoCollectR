#overlap
getProportionOverlap.OC.GP = function(mu_O, mu_C, sigma, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma, mu_D=NA, sigma_D=NA, mu_C=mu_C, sigma_C=NA, N=NA, n=NA, minResponse=minResponse, maxResponse=maxResponse)
	f <- function(x) dnorm(x, mu_O, sigma)
	g <- function(x) dnorm(x, mu_C, sigma)
	overlap <- integrate(function(x) pmin(f(x), g(x)), minResponse,maxResponse)$value
	return(overlap)
}
#peak
getPeak.T.GP = function(mu_O, mu_C) {
	parameter_checks(mu_O=mu_O, sigma_O=NA, mu_D=NA, sigma_D=NA, mu_C=mu_C, sigma_C=NA, N=NA, n=NA, minResponse=NA, maxResponse=NA)
	peak = mu_O + (mu_C - mu_O) / 2
	return(peak)
}
#Expected:
#first start
E.Ok1.GP = Vectorize(function(N, mu_O, sigma, minResponse=0, maxResponse=365, threshApprox=NA) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma, mu_D=NA, sigma_D=NA, mu_C=NA, sigma_C=NA, N=N, n=NA, minResponse=minResponse, maxResponse=maxResponse)
						 sigma_O = sigma
						 integrand <- function(x) { x * dOk1.GP(x, N=N, mu_O=mu_O, sigma=sigma) }
						 vals = integrate(integrand, minResponse, maxResponse, rel.tol = 1e-12, abs.tol = 1e-12, subdivisions = 10000L)$value

						 if(!is.na(threshApprox) && threshApprox>0) {
							 vals[is.na(vals)] = Inf
							 vals[is.nan(vals)] = Inf
							 vals.approx = E.Kth.approx(N=N, mu=mu_O, sigma=sigma_O, k=1)
							 vals[abs(vals - vals.approx) > threshApprox] = vals.approx[abs(vals - vals.approx) > threshApprox]
							 nRep = sum(abs(vals - vals.approx) > threshApprox)
							 if(nRep>0) {
								 warning(paste(nRep, " instances failed numerical integration and were replaced by an asymptotic approximation."))
							 }
						 }
						 return(vals)
})
#start
E.O.GP = function(mu_O) {
	parameter_checks(mu_O=mu_O, sigma_O=NA, mu_D=NA, sigma_D=NA, mu_C=NA, sigma_C=NA, N=NA, n=NA, minResponse=NA, maxResponse=NA)
	return(mu_O)
}
#observed
E.T.GP = function(mu_O, mu_D) {
	parameter_checks(mu_O=mu_O, sigma_O=NA, mu_D=mu_D, sigma_D=NA, mu_C=NA, sigma_C=NA, N=NA, n=NA, minResponse=NA, maxResponse=NA)
	return(mu_O + mu_D/2)
}

#cessation
E.C.GP = function(mu_C) {
	parameter_checks(mu_O=NA, sigma_O=NA, mu_D=NA, sigma_D=NA, mu_C=mu_C, sigma_C=NA, N=NA, n=NA, minResponse=NA, maxResponse=NA)
	return(mu_C)
}

#last cessation
E.CkN.GP = Vectorize(function(N, mu_C, sigma, minResponse=0, maxResponse=365, threshApprox=NA) {
	parameter_checks(mu_O=NA, sigma_O=NA, mu_D=NA, sigma_D=NA, mu_C=mu_C, sigma_C=sigma, N=N, n=NA, minResponse=minResponse, maxResponse=maxResponse)
						 sigma_C = sigma
						 integrand <- function(x) { x * dCkN.GP(x, N=N, mu_C=mu_C, sigma=sigma_C) }
						 vals = integrate(integrand, minResponse, maxResponse, rel.tol = 1e-12, abs.tol = 1e-12, subdivisions = 10000L)$value
						 #vals = integrate(integrand, minResponse, maxResponse)$value

						 if(!is.na(threshApprox) && threshApprox>0) {
							 vals[is.na(vals)] = Inf
							 vals[is.nan(vals)] = Inf
							 #warning("Using CkN approximation")
							 vals.approx = E.Kth.approx(N=N, mu=mu_C, sigma=sigma_C, k=N)
							 vals[abs(vals - vals.approx) > threshApprox] = vals.approx[abs(vals - vals.approx) > threshApprox]
							 nRep = sum(abs(vals - vals.approx) > threshApprox)
							 if(nRep>0) {
								 warning(paste(nRep, " instances failed numerical integration and were replaced by an asymptotic approximation."))
							 }
						 }

						 return(vals)
})

E.D.GP = function(mu_O, mu_C) {
	parameter_checks(mu_O=mu_O, sigma_O=NA, mu_D=NA, sigma_D=NA, mu_C=mu_C, sigma_C=NA, N=NA, n=NA, minResponse=NA, maxResponse=NA)
	return(mu_C - mu_O)
}

SD.Ok1.GP = function(N, mu_O, sigma, minResponse=0, maxResponse=365, intFailLow=NA, intFailHigh=NA) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma, mu_D=NA, sigma_D=NA, mu_C=NA, sigma_C=NA, N=N, n=NA, minResponse=minResponse, maxResponse=maxResponse)
	E = E.Ok1.GP(N, mu_O, sigma, minResponse, maxResponse)
	integrand <- function(x) { x * x * dOk1.GP(x, N, mu_O, sigma) }
	E2 = integrate(integrand, minResponse, maxResponse)$value
	retVal = E2
	if(!is.na(intFailLow)) {
		retVal[retVal < intFailLow] <- NA
	}
	if(!is.na(intFailHigh)) {
		retVal[retVal > intFailHigh] <- NA
	}
	E2 = retVal
	var = E2 - E^2
	if(var < 0) { stop("Error calculating variance for first onset times") }
	return(sqrt(var))
}

SD.O.GP = function(sigma) {
	parameter_checks(mu_O=NA, sigma_O=sigma, mu_D=NA, sigma_D=NA, mu_C=NA, sigma_C=NA, N=NA, n=NA, minResponse=NA, maxResponse=NA)
	return(sigma)
}

SD.T.GP = function(mu_O, mu_C, sigma, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma, mu_D=NA, sigma_D=NA, mu_C=mu_C, sigma_C=NA, N=NA, n=NA, minResponse=minResponse, maxResponse=maxResponse)
	E = E.T.GP(mu_O, mu_C)
	integrand <- function(x) { x * x * dT.GP(x, mu_O, mu_C, sigma, minResponse, maxResponse) }
	E2 = integrate(integrand, minResponse, maxResponse)$value
	var = E2 - E^2
	if(var < 0) {
		var = var(rT.GP(10000,mu_O, mu_C, sigma, minResponse, maxResponse))
		#stop("Error calculating variance for observed times")
	}
	return(sqrt(var))
}

SD.C.GP = function(sigma) {
	parameter_checks(mu_O=NA, sigma_O=NA, mu_D=NA, sigma_D=NA, mu_C=NA, sigma_C=sigma, N=NA, n=NA, minResponse=NA, maxResponse=NA)
	return(sigma)
}

SD.CkN.GP = function(N, mu_C, sigma, minResponse=0, maxResponse=365, intFailLow=NA, intFailHigh=NA) {
	parameter_checks(mu_O=NA, sigma_O=NA, mu_D=NA, sigma_D=NA, mu_C=mu_C, sigma_C=sigma, N=N, n=NA, minResponse=minResponse, maxResponse=maxResponse)
	E = E.CkN.GP(N, mu_C, sigma, minResponse, maxResponse)
	integrand <- function(x) { x * x * dCkN.GP(x, N, mu_C, sigma) }
	E2 = integrate(integrand, minResponse, maxResponse)$value
	retVal = E2
	if(!is.na(intFailLow)) {
		retVal[retVal < intFailLow] <- NA
	}
	if(!is.na(intFailHigh)) {
		retVal[retVal > intFailHigh] <- NA
	}
	E2 = retVal
	var = E2 - E^2
	if(var < 0) { stop("Error calculating variance for first onset times") }
	return(sqrt(var))
}

SD.D.GP = function() {
	warning("Duration has a Dirac delta distribution, so the variance is 0.")
	return(0)
}

PI.Ok1.GP = function(N, mu_O, sigma, alpha=0.05) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma, mu_D=NA, sigma_D=NA, mu_C=NA, sigma_C=NA, N=N, n=NA, minResponse=NA, maxResponse=NA)
	if(alpha<=0 || alpha>=1) {
		stop("Alpha must be between 0 and 1 exclusive.")
	}
	lower = alpha/2
	upper = 1 - alpha/2
	return(qOk1.GP(c(lower, upper), N, mu_O, sigma))
}

PI.O.GP = function(mu_O, sigma, alpha=0.05) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma, mu_D=NA, sigma_D=NA, mu_C=NA, sigma_C=NA, N=NA, n=NA, minResponse=NA, maxResponse=NA)
	if(alpha<=0 || alpha>=1) {
		stop("Alpha must be between 0 and 1 exclusive.")
	}
	lower = alpha/2
	upper = 1 - alpha/2
	return(qO.GP(c(lower, upper), mu_O, sigma))
}

PI.T.GP = function(mu_O, mu_C, sigma, alpha=0.05) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma, mu_D=NA, sigma_D=NA, mu_C=NA, sigma_C=NA, N=NA, n=NA, minResponse=NA, maxResponse=NA)
	if(alpha<=0 || alpha>=1) {
		stop("Alpha must be between 0 and 1 exclusive.")
	}
	lower = alpha/2
	upper = 1 - alpha/2
	return(qT.GP(c(lower, upper), mu_O, mu_C, sigma))
}

PI.C.GP = function(mu_C, sigma, alpha=0.05) {
	parameter_checks(mu_O=NA, sigma_O=NA, mu_D=NA, sigma_D=NA, mu_C=mu_C, sigma_C=sigma, N=NA, n=NA, minResponse=NA, maxResponse=NA)
	if(alpha<=0 || alpha>=1) {
		stop("Alpha must be between 0 and 1 exclusive.")
	}
	lower = alpha/2
	upper = 1 - alpha/2
	return(qC.GP(c(lower, upper), mu_C, sigma))
}

PI.CkN.GP = function(N, mu_C, sigma, alpha=0.05) {
	parameter_checks(mu_O=NA, sigma_O=NA, mu_D=NA, sigma_D=NA, mu_C=mu_C, sigma_C=sigma, N=N, n=NA, minResponse=NA, maxResponse=NA)
	if(alpha<=0 || alpha>=1) {
		stop("Alpha must be between 0 and 1 exclusive.")
	}
	lower = alpha/2
	upper = 1 - alpha/2
	return(qCkN.GP(c(lower, upper), N, mu_C, sigma))
}

PI.D.GP = function(mu_O, mu_C) {
	parameter_checks(mu_O=mu_O, sigma_O=NA, mu_D=NA, sigma_D=NA, mu_C=mu_C, sigma_C=NA, N=NN, n=NA, minResponse=NA, maxResponse=NA)
	d = mu_C - mu_O
	warning("Duration has a Dirac delta distribution, so 100% of the samples fall at the duration.")
	return(c(d,d))
}

#' @importFrom posterior as_draws_df
#' @importFrom dplyr filter %>%
runStan.NoCovariates.T.GP = function(fileOrData, minResponse=0, maxResponse=365, hyperparameters = c(100,7,50,7,10,7), dataProvided=FALSE, runMAP=TRUE, setStringent=FALSE, processExtremes=TRUE, N=500, maxDiv=0, threshApprox=NA, ...) {
	options(mc.cores = 4)

	cat("Starting Stan run.\n")

	#get and scale observed times
	if(dataProvided) {
		observed = (fileOrData - minResponse ) / (maxResponse - minResponse)
	}
	else {
		observed = getObservations(fileOrData,minResponse=minResponse,maxResponse=maxResponse)
	}

	if(processExtremes && !is.na(N)) {
		processExtremes = 1
	}
	else {
		N = 0
		processExtremes = 0
	}

	#set the hyperparameters
	mean_mean_onset = (hyperparameters[1] - minResponse) / (maxResponse-minResponse)
	sd_mean_onset = hyperparameters[2] / (maxResponse-minResponse)
	mean_mean_duration = hyperparameters[3] / (maxResponse-minResponse)
	sd_mean_duration = hyperparameters[4] / (maxResponse-minResponse)
	mean_sd = hyperparameters[5] / (maxResponse-minResponse) 		#taken just from the onset distribution sigma, as, in the GP model, duration is deterministic
	sd_sd = hyperparameters[6] / (maxResponse-minResponse) 			#taken just from the onset distribution sigma, as, in the GP model,  duration is deterministic

	#prepare data for Stan
	data <- list(
				 N = length(observed),	#sample size
				 n = N,
				 processExtremes = processExtremes,
				 t = observed,
				 minResponse = minResponse,
				 maxResponse = maxResponse,
				 debug = 0,
				 drop_nc = 0,
				 drop_ll = 0,
				 mean_mean_onset = mean_mean_onset,
				 sd_mean_onset = sd_mean_onset,
				 mean_mean_duration = mean_mean_duration,
				 sd_mean_duration = sd_mean_duration,
				 mean_sd = mean_sd,
				 sd_sd = sd_sd
	)

	cat("Attempting to compile Stan model.\n")
	m = tryCatch({
	  noCovariates.gp_file <- system.file("stan", "noCovariates.gp.stan"
	                                      , package = "phenoCollectR")
	  gp_file <- system.file("stan", "gp.hpp", package = "phenoCollectR")
		cmdstanr::cmdstan_model(stan_file = noCovariates.gp_file
		                        , stanc_options = list("allow-undefined")
		                        ,  user_header = gp_file)
	}, error = function(e) {
		ret = list(
				   error = TRUE,
				   errorM = e$message
		)
		return(ret)
	})


	#sample the posterior based on the provided observations and hyperparameter values
	cat("Attempting to sample Stan model.\n")
	res = tryCatch({
		if(setStringent) {
			m$sample(data = data, adapt_delta = 0.99, max_treedepth = 15)
		}
		else {
			m$sample(data = data)
		}
	}, error = function(e) {
		ret = list(
				   error = TRUE,
				   errorM = e$message
		)
		return(ret)
	})

	summary_ok <- tryCatch({
		res$summary()
		TRUE
		}, error = function(e) {
		FALSE
	})

	print("Starting summary")
	if(!summary_ok) {
		ret = list(
				   error = TRUE,
				   errorM = "Apparently Stan finished but with failed chains."
		)
		return(ret)
	}
	print("Finishing summary")

	cat("Checking divergences.\n")
	nDiv = sum(res$diagnostic_summary()$num_divergent)
	if(nDiv > maxDiv) {
		ret = list(
				   error = TRUE,
				   errorM = paste("Stan run failed. ", nDiv, " divergences encountered.")
		)
		return(ret)
	}


	if(runMAP) {
		#get GP model MAP
		cat("Attempting to find MAP parameter values.\n")
		resMAP = tryCatch({
			chain <- 1
			iter <- 999
			init_f_samp <- function() {
				init_df <- res$draws() %>% as_draws_df() %>% filter(.chain == chain, .iteration == iter)
				as.list(init_df)
			}
			m$optimize(data, init = init_f_samp)
		}, error = function(e) {
			ret = list(
					   error = TRUE,
					   errorM = e$message
			)
			return(ret)
		})
	}
	else { resMAP = NA }

	result = list(
				  minResponse = minResponse,
				  maxResponse = maxResponse,
				  runMAP = runMAP,
				  setStringent = setStringent,
				  N = N,
				  maxDiv = maxDiv,
				  threshApprox = threshApprox,
				  responseData = minResponse + observed * (maxResponse - minResponse),			#prescaled
				  hyperparameters = hyperparameters,			#prescaled
				  model = m, 				#Stan model
				  sample = res,				#HMC sample
				  MAP = resMAP,				#NA if runMap is false
				  data = data,				#data provided to Stan that includes the scaled hyperparameters and scaled observations
				  error = FALSE,
				  errorM = "Results were obtained, but check diagnostics of the sample."
	)
	return(result)
}

runStan.WithCovariates.T.GP = function(responseData, minResponse=0, maxResponse=365, onsetCovariateData, durationCovariateData, onsetHyperBeta, onsetHyperAnchor, durationHyperBeta, durationHyperAnchor, sigmaHyper=c(10,50), dataProvided=FALSE, setStringent=TRUE, maxDiv=0, processExtremes=TRUE, N=500, priorLevel=2, ...) {
	options(mc.cores = 4)

	if(processExtremes && !is.na(N)) {
		processExtremes = 1
	}
	else {
		N = 0
		processExtremes = 0
	}

	range = maxResponse-minResponse

	cat("Processing data.\n")
	if(dataProvided) {
		observed = (responseData - minResponse ) / (maxResponse - minResponse)
		covariatesOnset = processCovariates(onsetCovariateData)
		covariatesDuration = processCovariates(durationCovariateData)
	}
	else {
		#get scaled data
		observed = getObservations(responseData,minResponse=minResponse,maxResponse=maxResponse)
		covariatesOnset = getCovariates(onsetCovariateData) #min max scaled
		covariatesDuration = getCovariates(durationCovariateData) #min max scaled

		#get unscaled hyperparameters
		onsetHyperBeta = getHyperparameters(onsetHyperBeta)
		onsetHyperAnchor = getHyperparameters(onsetHyperAnchor)
		durationHyperBeta = getHyperparameters(durationHyperBeta)
		durationHyperAnchor = getHyperparameters(durationHyperAnchor)
		#cessationHyperAnchor = getHyperparameters(cessationHyperAnchor)
	}

	K_O = covariatesOnset$K
	K_D = covariatesDuration$K

	cov_O = colnames(covariatesOnset$covariates)
	cov_D = colnames(covariatesDuration$covariates)

	cov_union = union(cov_O, cov_D)

	idx_O = match(cov_O, cov_union)
	idx_D = match(cov_O, cov_union)

	#get the number of covariates for onset and for duration
	n = N	#switch variable naming
	N = length(observed)
	N_O = nrow(covariatesOnset$scaledCovariates)
	N_D = nrow(covariatesDuration$scaledCovariates)


	cat("Checking data compatibility.\n")
	#check data dimensions
	if(N != N_O || N != N_D) {
		return(list(error=TRUE, error_m="Be sure the sample size is the the same for the response (observed) values, for the covariate values for the onset, and for the covariate values for the duration. Rows should be parallel. Row 1 in the response data vector should correspond to the same individual in row 1 of each covariate data frame. If files names are input, corresponding files should be text formatted with tab-separated columns which should have headers. If data frames are input, these should have column names."))
	}

	#check hyperparameter specifications
	if(nrow(onsetHyperBeta) != K_O || ncol(onsetHyperBeta) != 2 || nrow(durationHyperBeta) != K_D || ncol(durationHyperBeta) != 2 || length(onsetHyperAnchor) != 2 || length(durationHyperAnchor) != 2 || length(sigmaHyper) != 2) {
		return(list(error=TRUE, error_m="The input hyperparameter information should be a data frame with one row for each covariate, or provide a file name with corresponding file having the data frame as a text, tab-separated table. Each row should have two columns. The first column provides the mean value of the parameter's prior distribution, and the second column provides the SD of the parameter's prior distribution. The covariates in rows, top to bottom, should match the covariates in columns, left to right, in the input covariate data frames. Values should be in the original scale (e.g., units of days, or for slopes, total days changed over range of the covariate) of the observations. If files names are input, corresponding files should have headers. If data frames are input, these should have column labels of 'mean
					and 'sd'."))

	}
	cat("Setting hyperparameters.\n")
	covariatesOnsetRanges = covariatesOnset$maxs - covariatesOnset$mins
	covariatesDurationRanges = covariatesDuration$maxs - covariatesDuration$mins

	meanSlopeO = (onsetHyperBeta[[1]]/range) * covariatesOnsetRanges
	sdSlopeO = (onsetHyperBeta[[2]]/range) * covariatesOnsetRanges
	meanAnchorO = ((onsetHyperAnchor[1]-minResponse)/range)
	sdAnchorO = (onsetHyperAnchor[2]/range)

	meanSlopeD = (durationHyperBeta[[1]]/range) * covariatesDurationRanges
	sdSlopeD = (durationHyperBeta[[2]]/range) * covariatesDurationRanges
	meanAnchorD = ((durationHyperAnchor[1])/range)
	sdAnchorD = (durationHyperAnchor[2]/range)

	#meanAnchorC = ((cessationHyperAnchor[1]-minResponse)/range)
	#sdAnchorC = (cessationHyperAnchor[2]/range)

	meanSigma = (sigmaHyper[1]/range)
	sdSigma = (sigmaHyper[2]/range)

	cat("Preparing data for Stan.\n")
	stanData <- list(
					 N = N,
					 K_O = K_O,
					 K_D = K_D,
					 K_union = length(cov_union),

					 idx_O = idx_O,
					 idx_D = idx_D,

					 priors = priorLevel,

					 T_raw = observed,
					 T_min = minResponse,
					 T_max = maxResponse,

					 process_extremes = processExtremes,
					 n = n,

					 min = minResponse,
					 maxResponse = maxResponse,

					 X_O_raw = as.matrix(covariatesOnset$scaledCovariates, ncol=K_O),
					 X_D_raw = as.matrix(covariatesDuration$scaledCovariates, ncol=K_D),

					 mean_X_O_raw = covariatesOnset$scaledMeans,
					 mean_X_D_raw = covariatesDuration$scaledMeans,

					 min_X_O = covariatesOnset$mins,
					 max_X_O = covariatesOnset$maxs,

					 min_X_D = covariatesDuration$mins,
					 max_X_D = covariatesDuration$maxs,

					 betaOnsetMeans = meanSlopeO,
					 betaOnsetSDs = sdSlopeO,
					 anchorOnsetMean = meanAnchorO,
					 anchorOnsetSD = sdAnchorO,

					 betaDurationMeans = meanSlopeD,
					 betaDurationSDs = sdSlopeD,
					 anchorDurationMean = meanAnchorD,
					 anchorDurationSD = sdAnchorD,

					 #anchorCessationMean = meanAnchorC,
					 #anchorCessationSD = sdAnchorC,

					 sigmaMean = meanSigma,
					 sigmaSD = sdSigma,

					 debug = FALSE,
					 drop_nc = FALSE,
					 drop_ll = FALSE
	)

	cat("Attempting to compile Stan model.\n")
	withCovariates.gp_file <- system.file("stan", "withCovariates.gp.stan"
	                                      , package = "phenoCollectR")
	gp_file <- system.file("stan", "gp.hpp", package = "phenoCollectR")
	m = tryCatch({
		cmdstanr::cmdstan_model(stan_file = withCovariates.gp_file
		                        , stanc_options = list("allow-undefined")
		                        , user_header = gp_file)
	}, error = function(e) {
		message <- conditionMessage(e)
		cat(paste("Error compiling Stan model: ", message, ".\n"))
		ret = list(
				   error = TRUE,
				   errorM = paste(message, ": Could not compile Stan model.")
		)
		return(ret)
	})

	cat("Attempting to run Stan.\n")
	res = tryCatch({
		fit = if(setStringent) {
			mres = m$sample(data = stanData, adapt_delta = 0.99, max_treedepth = 15,...)
		}
		else {
			mres = m$sample(data = stanData,...)
		}

	}, error = function(e) {
		cat(paste("ERROR running Stan: ", e$message, ".\n"))
		ret = list(
				   error = TRUE,
				   errorM = paste(e$message, ": Could not sample Stan model"),
				   data = stanData
		)
		return(ret)
	})

	nDiv = 50000
	chain_ok <- tryCatch({
		s = res$diagnostic_summary()
		print(s)
		nDiv = sum(s$num_divergent)
		TRUE
	}, error = function(e) {
		FALSE
	})

	if(!chain_ok) {
		cat("Stan finished but with no summary available.\n")
		ret = list(
				   error = TRUE,
				   errorM = "Apparently Stan finished but with failed chains."
		)
		return(ret)
	}

	cat("Checking divergences. (with covariates)\n")
	if(nDiv > maxDiv) {
		ret = list(
				   error = TRUE,
				   errorM = paste("Stan run failed. ", nDiv, " divergences encountered.")
		)
		return(ret)
	}

	output = list(
				  data = stanData,
				  responseData = responseData,
				  onsetCovariateData = onsetCovariateData,
				  durationCovariateData = durationCovariateData,
				  onsetHyperBeta = onsetHyperBeta,
				  onsetHyperAnchor = onsetHyperAnchor,
				  durationHyperBeta = durationHyperBeta,
				  durationHyperAnchor = durationHyperAnchor,
				  #cessationHyperAnchor = cessationHyperAnchor,
				  sigmaHyper = sigmaHyper,
				  result=res,
				  model=m,
				  maxDiv = maxDiv,
				  minResponse = minResponse,
				  maxResponse = maxResponse,
				  setStringent = setStringent,
				  error = FALSE,
				  error_m = "No error was detected during the Stan run"
	)
	return(output)
}

