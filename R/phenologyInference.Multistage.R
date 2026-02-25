runStan.WithCovariates.Multistage.durations.GP = function(responseData=NULL, stage=NULL, minResponse=0, maxResponse=365, covariateData=NULL, onsetHyperBeta=NULL, onsetHyperAnchor=NULL, durationHyperBetaMean=NULL, durationHyperBetaSD=NULL, durationHyperAnchor=NULL, sigmaHyper=NULL, setStringent=TRUE, priorLevel=1, maxDiv=0, processExtremes=TRUE, N=500, debug=FALSE, ...) {

	if(minResponse != 0) {
		stop("minResponse must be 0.")
	}

	range = maxResponse - minResponse

	options(mc.cores = 4)

	print("Running runStan.WithCovariates.multistage.durations.GP")

	if(processExtremes && !is.null(N)) {
		if(N>0) {
		processExtremes = 1
		}
	}
	else {
		N = 0
		processExtremes = 0
	}


	cat("Processing data.\n")
	observed = responseData
	covariates = processCovariates(covariateData)

	if(nrow(durationHyperBetaSD) != nrow(durationHyperBetaMean)) {
		stop("The number of stages (rows) in the durationHyperBetaSD and durationHyperBetaMean matrices need to be the same.")
	}

	S = nrow(durationHyperBetaSD)+1; #before first onset is another stage to avoid wraparound
	K = covariates$K

	cov_names = colnames(covariates$scaledCovariates)

	#get the number of covariates for onset and for duration
	n = N	#switch variable naming, this one for population size, not sample size
	N = length(observed)
	N_C = nrow(covariates$scaledCovariates)
	N_S = length(stage) 

	cat("Checking data compatibility.\n")
	#check that stage values match
	if(N != N_S) {
		return(list(error=TRUE, error_m="When running the multistage model, please include a vector of stages of the individuals that were observed. The vector of stages should be the same length as the vector of response data. For each individual sampled, there should be a response time and there should be a stage associated with the time. If the time is before the fist phenophase of the time period, the stage is the last stage, if during the first phenophase, the stage is 1, and so forth."))
	}
	#check data dimensions
	if(N != N_C) {
		return(list(error=TRUE, error_m="Be sure the sample size is the the same for the response (observed) values and for the covariate values. Rows should be parallel. Row 1 in the response data vector should correspond to the same individual in row 1 of the covariate data frame. The input data frames should have column names."))
	}

	#check hyperparameter specifications
	if(nrow(onsetHyperBeta) != K || ncol(onsetHyperBeta) != 2 || nrow(durationHyperBetaMean) != (S-1) || ncol(durationHyperBetaMean) != K || nrow(durationHyperBetaSD) != (S-1) || ncol(durationHyperBetaSD) != K || length(onsetHyperAnchor) != 2 || nrow(durationHyperAnchor) != S-1 || ncol(durationHyperAnchor) != 2 || ncol(sigmaHyper) != 2 || nrow(sigmaHyper)!=S) {
		return(list(error=TRUE, error_m="Hyperparameters should be as follows: onsetHyperAnchor: a 2-element vector with mean and sd of the onset time; onsetHyperBeta: a data frame with two columns (mean and sd of each covariate slope of the onset model), one covariate per row; durationHyperAnchor: a data frame with two columns (mean and sd of the duration), one pair for each stage(rows), excluding the last stage (S-1 rows); durationHyperBetaMean: a matrix with S-1 rows and K columns with each stage X covariate mean slope in cells; durationHyperBetaSD: a matrix with S-1 rows and K columns with each stage X covariate slope sd in cells; sigmaHyper: a data frame with two colums (mean and sd of each model sigma) and S rows, one for each stage."))

	}

	cat("Setting hyperparameters.\n")
	#covariateRanges = covariates$maxs - covariates$mins
	covariateRanges = 1.0
	range = 1.0

	meanSlopeO = (onsetHyperBeta[[1]]/range) * covariateRanges
	sdSlopeO = (onsetHyperBeta[[2]]/range) * covariateRanges
	meanAnchorO = ((onsetHyperAnchor[1]-minResponse)/range)
	sdAnchorO = (onsetHyperAnchor[2]/range)

	meanAnchorD = durationHyperAnchor[[1]]/range
	sdAnchorD = durationHyperAnchor[[2]]/range

	meanSigma = sigmaHyper[[1]]/range
	sdSigma = sigmaHyper[[2]]/range

	meanSlopeD = durationHyperBetaMean/range
	meanSlopeD = sweep(meanSlopeD, 2, covariateRanges, `*`)
	sdSlopeD = durationHyperBetaSD/range
	sdSlopeD = sweep(sdSlopeD, 2, covariateRanges, `*`)

	cat("Preparing data for Stan.\n")
	stanData <- list( 
	 #debug = debug,				#Enable debugging features in Stan (scalar)
	 #drop_ll = FALSE,			#Do not include likelihoods, for prior predictives (scalar)
	 #priors = priorLevel,			#Type of prior (0 = flat, other = normal)
	 #process_extremes = processExtremes,	#Process extremes (scalar)
	 #n = n,					#Population size to estimate extremes (scalar)
	 N = N,					#Sample size (scalar)
	 t_raw = observed,			#Observed collection times, scaled (N vector)
	 #T_min = minResponse,			#Minimum collection time, original scale (scalar) ASSUMED TO BE 0
	 T_max = maxResponse,			#Maximum collection time, original scale (scalar)
	 #S = S,					#Number of stages (scalar)
	 S = S+1,					#Number of stages (scalar) - no wraparound -> S+1
	 stage = stage,				#Observed stages (N vector)
	 K = K,					#Number of covariates (scalar)
	 #X_raw = as.matrix(covariates$scaledCovariates, ncol=K),	#The scaled covariate data (N X K matrix)
	 #mean_X_raw = covariates$scaledMeans,	#The scaled covariate means (C vector)
	 #min_X = covariates$mins,		#The minimum values of covariates, original scale (C vector)
	 #max_X = covariates$maxs,		#The maximum values of covariates, original scale (C vector)
	 #sigmaMean = meanSigma,			#Hyperparameter mean sigma (S vector)
	 #sigmaSD = sdSigma,			#Hyperparameter sd sigma (S vector)
	 #anchor1Mean = meanAnchorO,		#Hyperparameter mean onset, stage 1 (scalar)
	 #anchor1SD = sdAnchorO,			#Hyperparameter sd onset, stage 1 (scalar)
	 #beta1Mean = meanSlopeO,		#Hyperparameter mean slopes, stage 1 (C vector)
	 #beta1SD = sdSlopeO,			#Hyperparameter sd slopes, stage 1 (C vector)
	 #anchorDurationMean = meanAnchorD,	#Hyperparameter mean duration, without last stage (C vector)
	 #anchorDurationSD = sdAnchorD,		#Hyperparameter sd duration, without last stage (C vector)
	 #betaDurationMean = meanSlopeD,		#Hyperparameter mean duration, without last stage (S-1 X C matrix)
	 #betaDurationSD = sdSlopeD,		#Hyperparameter sd duration, without last stage (S-1 X C matrix)
	 X_raw = as.matrix(covariates$covariates, ncol=K)	#The scaled covariate data (N X K matrix)
	)

	cat("Attempting to compile Stan model in file withCovariates.gp.multistage.durations.stan.\n")
	withCovariates.multistage.durations.gp_file <- system.file("stan", "withCovariates.gp.multistage.durations.stan"
	                                      , package = "phenoCollectR")
	m = tryCatch({
		cmdstanr::cmdstan_model(stan_file = withCovariates.multistage.durations.gp_file)
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
	      data = stanData,				#Data passed to Stan
	      responseData = responseData,		#Original response data
	      minResponse = minResponse,		#Minimum response
	      maxResponse = maxResponse,		#Maximum response
	      covariateData = covariateData,		#Original covariateData
	      onsetHyperBeta = onsetHyperBeta,		#Prior for stage 1 slopes
	      onsetHyperAnchor = onsetHyperAnchor,	#Prior for stage 1 anchor
	      durationHyperBetaMean = durationHyperBetaMean,	#Prior for mean duration slopes
	      durationHyperBetaSD = durationHyperBetaSD,	#Prior for sd duration slopes
	      durationHyperAnchor = durationHyperAnchor,	#Prior for duration
	      sigmaHyper = sigmaHyper,			#Prior for onset and duration SDs
	      result=res,				#The output from the Stan run
	      model=m,					#The model used in Stan
	      priorLevel=priorLevel,			#Which prior to use (flat / normal)
	      processExtremes=processExtremes,		#Whether extremes were processed
	      N=n,					#Population size used for extremes
	      maxDiv = maxDiv,				#Max tolerated divergences
	      setStringent = setStringent,		#Stringent run?
	      error = FALSE,				#Errorr?
	      error_m = "No error was detected during the Stan run"	#Error message
	)
	return(output)
}
