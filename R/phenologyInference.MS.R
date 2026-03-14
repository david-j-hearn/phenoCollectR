#' @importFrom copula pobs normalCopula fitCopula 
runStan.WithCovariates.Multistage.durations.GP = function(responseData=NULL, stage=NULL, nStages, nCovariates, minResponse=0, maxResponse=365, covariateData=NULL, onsetHyperBeta=NULL, onsetHyperAnchor=NULL, durationHyperBetaMean=NULL, durationHyperBetaSD=NULL, durationHyperAnchor=NULL, sigmaHyper=NULL, setStringent=TRUE, priorLevel=1, maxDiv=0, calculatePPD=FALSE, nXs=100, nReps=10, targetCovariateName="cov1", debug=FALSE, ...) {

	#checkInput and checkPriors should have already been called from runStanPhenology


	if(minResponse!=0) {
		stop("Minimum response must be 0")
	}

	xPPD <- as.data.frame(matrix(NA_real_, nrow = nReps*nXs , ncol = nCovariates))

	if(calculatePPD && nXs>0 && nReps>0 && !is.null(targetCovariateName)) {
		if(nCovariates<=1) {
			stop("Currently, PPD can only be estimated for models with more than 1 covariate.")
		}

		cat("Calculating PPD. This will lengthen the Stan runtime.\n")
		cat("Setting up copula for sampling covariates.\n")

		target_idx = which(colnames(covariateData) == targetCovariateName)

		U  =  pobs(as.matrix(covariateData))  # Pseudo-observations (uniform marginals)
		cop  =  normalCopula(dim = ncol(U), dispstr = "un")
		fitC  =  fitCopula(cop, U, method = "ml")
		cov_names <- paste0("X", 1:nCovariates)
		colnames(xPPD) <- cov_names

		minX = min(covariateData[targetCovariateName])
		maxX = max(covariateData[targetCovariateName])
		Xs = seq(from=minX, to=maxX, length.out=nXs)
		for(i in 1:nXs) {
			startInd = (i-1)*nReps + 1
			endInd = i*nReps
			xPPD[startInd:endInd,] = sample_conditional_covariates(
				       x_target = Xs[i],
				       x_column = target_idx,
				       covars = covariateData,
				       copula_fit = fitC,
				       n_samples = nReps
				       )[1:nReps,]
		}
	}
	else if(calculatePPD) {
		stop("Please make sure the number of increments along the x-axis (nXs) and the number of replicates (nReps) are both 1 or greater. Also, provide the name of the covariate that you want to appear on the x-axis (targetCovariateName)")
	}
	else if(!calculatePPD){
		nReps=1
		nXs=1
		xPPD <- as.data.frame(matrix(rep(0,nReps*nXs*nCovariates), nrow = nReps*nXs , ncol = nCovariates))
	}

	options(mc.cores = 4)
	cat("Running runStan.WithCovariates.multistage.durations.GP\n")

	cat("Processing data.\n")
	observed = responseData
	covariates = processCovariates(covariateData)

	S = nStages
	K = nCovariates

	#get the number of covariates for onset and for duration
	N = length(observed)

	cat("Setting hyperparameters. To be done in Stan\n")

	betaMeans = rbind(onsetHyperBeta[[1]], durationHyperBetaMean)
	betaSDs = rbind(onsetHyperBeta[[2]], durationHyperBetaSD)

	anchorMeans = c(onsetHyperAnchor[1],durationHyperAnchor[[1]])
	anchorSDs = c(onsetHyperAnchor[2],durationHyperAnchor[[2]])

	sigmaMean = sigmaHyper[1]
	sigmaSD = sigmaHyper[2]

	cat("Preparing data for Stan.\n")
	stanData <- list( 
			 #debug = debug,				#Enable debugging features in Stan (scalar)
			 #drop_ll = FALSE,			#Do not include likelihoods, for prior predictives (scalar)
			 #priors = priorLevel,			#Type of prior (0 = flat, other = normal)
			 #process_extremes = processExtremes,	#Process extremes (scalar)
			 #n = n,					#Population size to estimate extremes (scalar)
			 N = N,					#Sample size (scalar)
			 t_raw = observed,			#Observed collection times, scaled (N vector)
			 T_max = maxResponse,			#Maximum collection time, original scale (scalar)
			 stage = stage,				#Observed stages (N vector)
			 S = S+1,					#Number of stages (scalar) - no wraparound -> S+1
			 X_raw = as.matrix(covariates$covariates, ncol=K),	#The scaled covariate data (N X K matrix)
			 K = K,					#Number of covariates (scalar)
			 ppd = calculatePPD,			#Boolean to calculate posterior predictive data
			 nXs = nXs,				#number of x iterations of target covariate for PPD
			 nReps = nReps,				#number of replicates per Stan sample per stage per nXs iterations
			 xPPD = as.matrix(xPPD, ncol=K),	#the data over which to calculate PPD
			 betaMeans = betaMeans,			#joined onset and duration model slope coefficients
			 betaSDs = betaSDs, 			#joined onset and duration model sd on slope coefficients
			 anchorMeans = anchorMeans, 		#joined onset and duration model anchors (with standardized data these are the marginal mean responses)
			 anchorSDs = anchorSDs, 		#joined onset and duration model sd on anchors 
			 sigmaMean = sigmaMean,
			 sigmaSD = sigmaSD
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
		      #priorLevel=priorLevel,			#Which prior to use (flat / normal)
		      #processExtremes=processExtremes,		#Whether extremes were processed
		      #N=n,					#Population size used for extremes
		      maxDiv = maxDiv,				#Max tolerated divergences
		      setStringent = setStringent,		#Stringent run?
		      error = FALSE,				#Errorr?
		      error_m = "No error was detected during the Stan run"	#Error message
	)
	return(output)
}

#' @importFrom copula pobs normalCopula fitCopula 
runStan.WithCovariates.Multistage.overlap.GP = function(responseData=NULL, stageCounts=NULL, nStages, nCovariates, minResponse=0, maxResponse=365, covariateData=NULL, onsetHyperBeta=NULL, onsetHyperAnchor=NULL, durationHyperBetaMean=NULL, durationHyperBetaSD=NULL, durationHyperAnchor=NULL, sigmaHyper=NULL, setStringent=TRUE, priorLevel=1, maxDiv=0, calculatePPD=FALSE, nXs=100, nReps=10, targetCovariateName="cov1", debug=FALSE, ...) {

	#checkInput and checkPriors should have already been called from runStanPhenology

  #Stages are tricky - the onset of the first stage is not discoverable from this approach; instead, the start of the period through to the second stage is the first stage, so it may overlap with pre-reproductive and the first reproductive stage, and the last stage extends from the last reproductive period through to the period after reproduction.

	if(minResponse!=0) {
		stop("Minimum response must be 0")
	}

	#xPPD <- as.data.frame(matrix(NA_real_, nrow = nReps*nXs , ncol = nCovariates))

	#if(calculatePPD && nXs>0 && nReps>0 && !is.null(targetCovariateName)) {
		#if(nCovariates<=1) {
			#stop("Currently, PPD can only be estimated for models with more than 1 covariate.")
		#}

		#cat("Calculating PPD. This will lengthen the Stan runtime.\n")
		#cat("Setting up copula for sampling covariates.\n")
#
		#target_idx = which(colnames(covariateData) == targetCovariateName)

		#U  =  pobs(as.matrix(covariateData))  # Pseudo-observations (uniform marginals)
		#cop  =  normalCopula(dim = ncol(U), dispstr = "un")
		#fitC  =  fitCopula(cop, U, method = "ml")
		#cov_names <- paste0("X", 1:nCovariates)
		#colnames(xPPD) <- cov_names
#
		#minX = min(covariateData[targetCovariateName])
		#maxX = max(covariateData[targetCovariateName])
		#Xs = seq(from=minX, to=maxX, length.out=nXs)
		#for(i in 1:nXs) {
			#startInd = (i-1)*nReps + 1
			#endInd = i*nReps
			#xPPD[startInd:endInd,] = sample_conditional_covariates(
				       #x_target = Xs[i],
				       #x_column = target_idx,
				       #covars = covariateData,
				       #copula_fit = fitC,
				       #n_samples = nReps
				       #)[1:nReps,]
		#}
	#}
	#else if(calculatePPD) {
		#stop("Please make sure the number of increments along the x-axis (nXs) and the number of replicates (nReps) are both 1 or greater. Also, provide the name of the covariate that you want to appear on the x-axis (targetCovariateName)")
	#}
  if(calculatePPD) {
    stop("Posterior predictive plotting is not yet implemented for within-individual overlapping multistage analyses.")
  }
	else if(!calculatePPD){
		nReps=1
		nXs=1
		xPPD <- as.data.frame(matrix(rep(0,nReps*nXs*nCovariates), nrow = nReps*nXs , ncol = nCovariates))
	}

	options(mc.cores = 4)
	cat("Running runStan.WithCovariates.multistage.overlap.GP\n")

	cat("Processing data.\n")
	observed = responseData
	covariates = processCovariates(covariateData)

	S = nStages
	K = nCovariates

  #stageCounts = cbind(rep(0,nrow(stageCounts)),stageCounts)

	#get the number of covariates for onset and for duration
	N = length(observed)

	cat("Setting hyperparameters. Scaling to be done in Stan\n")



	betaMeans = rbind(onsetHyperBeta[[1]], durationHyperBetaMean)
	betaSDs = rbind(onsetHyperBeta[[2]], durationHyperBetaSD)

	anchorMeans = c(onsetHyperAnchor[1],durationHyperAnchor[[1]])
	anchorSDs = c(onsetHyperAnchor[2],durationHyperAnchor[[2]])

	sigmaMean = sigmaHyper[1]
	sigmaSD = sigmaHyper[2]

print("betaMeans")
print(betaMeans)
print("betaSDs")
print(betaSDs)
print("anchorMeans")
print(anchorMeans)
print("anchorSDs")
print(anchorSDs)


	cat("Preparing data for Stan.\n")
	stanData <- list( 
			 #debug = debug,				#Enable debugging features in Stan (scalar)
			 #drop_ll = FALSE,			#Do not include likelihoods, for prior predictives (scalar)
			 #priors = priorLevel,			#Type of prior (0 = flat, other = normal)
			 #process_extremes = processExtremes,	#Process extremes (scalar)
			 #n = n,					#Population size to estimate extremes (scalar)
			 N = N,					#Sample size (scalar)
			 t_raw = observed,			#Observed collection times, scaled (N vector)
			 T_max = maxResponse,			#Maximum collection time, original scale (scalar)
			 stage_counts = stageCounts,				#Observed stage counts of developmental unit (N X S matrix)
			 S = S,					#Number of stages (scalar) - no wraparound
			 X_raw = as.matrix(covariates$covariates, ncol=K),	#The scaled covariate data (N X K matrix)
			 K = K,					#Number of covariates (scalar)
			 #ppd = calculatePPD,			#Boolean to calculate posterior predictive data
			 #nXs = nXs,				#number of x iterations of target covariate for PPD
			 #nReps = nReps,				#number of replicates per Stan sample per stage per nXs iterations
			 #xPPD = as.matrix(xPPD, ncol=K),	#the data over which to calculate PPD
			 betaMeans = betaMeans,			#joined onset and duration model slope coefficients
			 betaSDs = betaSDs, 			#joined onset and duration model sd on slope coefficients
			 anchorMeans = anchorMeans, 		#joined onset and duration model anchors (with standardized data these are the marginal mean responses)
			 anchorSDs = anchorSDs, 		#joined onset and duration model sd on anchors 
			 sigmaMean = sigmaMean,
			 sigmaSD = sigmaSD
	)

	cat("Attempting to compile Stan model in file withCovariates.gp.multistage.overlap.stan.\n")
	withCovariates.multistage.overlap.gp_file <- system.file("stan", "withCovariates.gp.multistage.overlap.stan"
								   , package = "phenoCollectR")
	m = tryCatch({
		cmdstanr::cmdstan_model(stan_file = withCovariates.multistage.overlap.gp_file)
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
		      #priorLevel=priorLevel,			#Which prior to use (flat / normal)
		      #processExtremes=processExtremes,		#Whether extremes were processed
		      #N=n,					#Population size used for extremes
		      maxDiv = maxDiv,				#Max tolerated divergences
		      setStringent = setStringent,		#Stringent run?
		      error = FALSE,				#Errorr?
		      error_m = "No error was detected during the Stan run"	#Error message
	)
	return(output)
}
