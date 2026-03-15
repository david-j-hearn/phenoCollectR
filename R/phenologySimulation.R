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
#' @keywords internal
#'
#' @examples
#' \donttest{
#' #Set the parameters
#' n=1000 #sample size
#' slopeO = 1 #set the slope of the onset model
#' slopeD = 0.5 #set the slope of the duration model
#' interceptO = 100 #set the intercept of the onset model
#' interceptD = 20 #set the intercept of the duration mdoel
#' sigma = 7 #set the standard deviation of the onset distribution and of the cessation distribution
#' minCovariate = -2 #set the minimum value of the covariate
#' maxCovariate = 30 #set the maximum value of the covariate
#' data = simulateCovariate(n=n, slopeO=slopeO, interceptO=interceptO, sigma=sigma, slopeD=slopeD
#'                          , interceptD=interceptD, minCovariate=minCovariate
#'                          , maxCovariate=maxCovariate)
#' #plot the simulated observed collection times
#' plot(data$X, data$Ts, xlab="Mean spring temperature", ylab="Day of year", col="purple"
#'      , main=NULL, pch=16)
#' points(data$X, data$O, col="red", pch=16, cex=0.3)
#' points(data$X, data$C, col="blue", pch=16, cex=0.3)
#' #Plot the line that passes through the mean observed collection times, onset and cessation
#' abline(a = interceptO + interceptD/2, b = slopeO + slopeD/2, col = "purple", lwd = 2)
#' abline(a = interceptO, b = slopeO, col = "red", lwd = 2)
#' abline(a = interceptO + interceptD, b = slopeO + slopeD, col = "blue", lwd = 2)
#' #plot phenophases for each individual in gray
#' segments(x0 = data$X, y0 = data$O, x1 = data$X, y1 = data$C, col = "gray", lwd=0.5)
#' }
#' @noRd
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

#' Simulate phenological times in a population with unimodal phenophase
#'
#' @description Simulate a population of individuals with simulated onset, duration, cessation, and observed times based on multiple correlated covariates for onset and for duration. Assumes a linear model with normally distributed variation. Defaults produced simulated values with no covariates.
#'
#' @param n Sample size
#' @param betaOnset Named vector of slope coefficients of the covariates for onset. Names are the names of the covariates. (default: NULL)
#' @param betaDuration Named vector of slope coefficients of the covariates for duration. Names are the names of the covariates. (default: NULL)
#' @param covariateNamesOnset Vector of the names of the covariates for onset; must be same length as betaOnset  (default: NULL)
#' @param covariateNamesDuration Vector of the names of the covariates for duration; must be same length as betaDuration  (default: NULL)
#' @param covariateMeans Named vector of the means of the covariates for onset. Names are the names of all the covariate variables. (default: NULL)
#' @param covarianceMatrix The covariance matrix for all covariates with named rows and columns. Must be real, positive semidefinite. Either the covariance matrix or the correlation matrix is needed. Names are the names of all the covariate variables. (default: NULL)
#' @param correlationMatrix The correlation matrix for all covariates with named rows and columns. Either the covariance matrix or the correlation matrix is needed. If the correlation matrix is used, the standard deviations of the covariates need to be supplied in the named covariateStandardDeviations vector parameter. Names are the names of all the covariate variables.  (default: NULL)
#' @param covariateStandardDeviations Named vector of the standard deviations of all the covariate variables. Names are the names of all the covariate variables. (default: NULL)
#' @param meanOnset Marginal mean value of the onset times. (default: 0)
#' @param meanDuration Marginal mean value of the duration times. (default: 0)
#' @param sigmaOnset Standard deviation of the onset times. Duration is assumed fixed in this model. (default: 1)
#' @param minResponse The minimum time of observations. default: 0, which is the only fully supported value.
#' @param maxResponse The maximum time of observations. default: 365, representing the number of days in a year. 
#'
#' @return Dataframe with columns labeled with the variable names, including all covariates, the simulated onset times, durations, cessations, sampled times, stage (before, during, after phenophase) and rows representing simulated individuals
#' @importFrom MASS mvrnorm
#'
#' @examples
#' \donttest{
#' #Set the model parameters
#' #Variable names
#' covariate_namesOnset = c("o1", "o2", "c1") 
#' covariate_namesDuration = c("d1", "c1") 
#' allCovariates = union(covariate_namesOnset, covariate_namesDuration) #onset first, then duration
#' #Set up parameters to simulate covariate values
#' correlationMatrix = matrix(c( 1.0, 0.5, 0.3, 0.1, 
#'			              0.5, 1.0, 0.4, 0.2,
#'				      0.3, 0.4, 1.0, 0.0,
#'				      0.1, 0.2, 0.0, 1.0), nrow = 4
#'                            , byrow = TRUE)
#' colnames(correlationMatrix) = allCovariates
#' rownames(correlationMatrix) = allCovariates
#' covariateStandardDeviations = c(1,2,2,5)
#' names(covariateStandardDeviations) = allCovariates
#' covariateMeans = c(10,20,30,10)
#' names(covariateMeans) = allCovariates
#' #Set up the onset parameters
#' slopesOnset = c(1,2,3)
#' names(slopesOnset) = covariate_namesOnset
#' mean_responseOnset = 150
#' noiseOnset = 3
#' #Set up the duration parameters
#' slopesDuration = c(1,3)
#' names(slopesDuration) = covariate_namesDuration
#' mean_responseDuration = 30
#' #Sample size
#' n=1000
#' #Simulate the data
#' simulated_data = simulatePopulationLatentIntervalStates(n=n,
#' 			betaOnset=slopesOnset, betaDuration=slopesDuration,
#'			covariateNamesOnset=covariate_namesOnset, covariateNamesDuration=covariate_namesDuration,
#'			covariateMeans = covariateMeans,
#'			correlationMatrix = correlationMatrix,
#'			covariateStandardDeviations = covariateStandardDeviations,
#'			meanOnset = mean_responseOnset, meanDuration = mean_responseDuration,
#'			sigmaOnset = noiseOnset)
#' #Make a scatter plot of the simulated data
#' plot(simulated_data$o1[simulated_data$stage==2], simulated_data$observedTime[simulated_data$stage==2], main=NULL, xlab="O1", ylab="Onset")
#' #Check simulated values against actual
#' #Should be o1: 10, o2: 20, c1: 30, d1: 10, onsets: 150, durations: 30, cessations: 180, observedTimes: 366/2, stages: around 2
#' colMeans(simulated_data)
#' attach(simulated_data)
#' #The correlation matrix should be similar to correlationMatrix defined above
#' cor(data.frame(o1,o2,c1,d1))
#' #Coefficients for the onset model should be slopesOnset: 1, 2, 3
#' summary(lm(onset ~ o1 + o2 + c1))
#' #Coefficients for the duration model should be slopesDuration: 1, 3
#' summary(lm(duration ~ d1 + c1))
#' #The observed time coefficients should be betaOnset + betaDuration/2
#' #With the above numbers that's o1: 1, o2: 2, c1: (3 + 3/2), d1: 1/2
#' summary(lm(observedTime[stage==2] ~ o1[stage==2] + o2[stage==2] + c1[stage==2] + d1[stage==2]))
#' #The cessation coefficients should be betaOnset + betaDuration
#' #With the above numbers that's o1: 1, o2: 2, c1: (3 + 3), d1: 1
#' summary(lm(cessation ~ o1 + o2 + c1 + d1))
#' }
#' @noRd
simulatePopulationLatentIntervalStates = function(n, 
					covariateNamesOnset=NULL, covariateNamesDuration=NULL, 
					covariateMeans = NULL, 
					covarianceMatrix = NULL,
					correlationMatrix = NULL, 
					covariateStandardDeviations = NULL,
					betaOnset=NULL, betaDuration=NULL, 
					meanOnset = 0, meanDuration = 0, 
					sigmaOnset = 1, 
					minResponse=0, maxResponse=365,
					seed = NULL) {

	#Basic checks
	if(n<=0) {
		stop("n should be a positive integer.")
	}
	if(meanOnset <= minResponse || meanDuration <= minResponse ) {
		stop("Mean onset and mean duration must be greater than 0")
	}
	if(meanOnset >= maxResponse) {
		stop("Mean onset must be before the end of the time period")
	}
	if(meanDuration >= maxResponse-meanOnset) {
		stop("Mean duration must occur within the window between the mean onset time and the end of the time period.")
	}

	#Simulate covariate data
	covs_all = union(covariateNamesOnset,covariateNamesDuration)
	if(!is.null(covs_all)) {
	X = simulateCorrelatedCovariateData(n=n, covariateNames=covs_all, covariateMeans=covariateMeans, Sigma = covarianceMatrix, R = correlationMatrix, covariateSDs = covariateStandardDeviations, seed = seed)  
	data = X
	}
	else {
		data = data.frame(matrix(nrow = n, ncol = 0))
	}

	#Simulate onset data
	if(is.null(covariateNamesOnset)) {
		onset = rnorm(n,meanOnset,sigmaOnset)
		onsetData = data.frame(onset=onset)
	}
	else {
		#Pull out data for onset
		onsetCovariates = X[, covariateNamesOnset, drop=FALSE]; 
		onsetData = simulateResponseData(X=onsetCovariates, beta=betaOnset, response_name = "onset", anchor = meanOnset, noise_sd = sigmaOnset) 
	}

	#Simulate duration data
	if(is.null(covariateNamesDuration)) {
		duration = rep(meanDuration,n)
		durationData = data.frame(duration=duration)
	}
	else {
		#Pull out data for duration
		durationCovariates = X[, covariateNamesDuration, drop=FALSE]; 
		durationData = simulateResponseData(X=durationCovariates, beta=betaDuration, response_name = "duration", anchor = meanDuration, noise_sd = 0) 
	}

	#Simulate observed times based on simulated onset and duration times
	onset = onsetData$onset
	duration = durationData$duration
	cessation = onset+duration
	observedTime = runif(n,minResponse,maxResponse)
	stage = 1L + (observedTime >= onset) + (observedTime >= cessation)

	#Prepare output data frame
	data$onset = onset
	data$duration = duration
	data$cessation = cessation
	data$observedTime = observedTime
	data$stage = stage
	return(data)
}

#' Simulate a response variable based on a linear model of covariates
#'
#' @description Simulate the values of a response variable based on values of multiple covariates. Assumes a linear model with normally distributed variation.
#'
#' @param X A data frame with labeled columns indicating covariate names and rows with covariate data
#' @param beta Vector of slope coefficients of the covariates
#' @param response_name The name of the response variable (default: "Y")
#' @param Sigma The covariance matrix. Must be of dimension length(beta) X length(beta) (default: identity matrix)
#' @param anchor Marginal mean value of the response variable
#' @param noise_sd Standard deviation of the noise in the response variable
#'
#' @return A data frame with columns labeled with the variable names and rows representing simulation replicates.
#' @importFrom MASS mvrnorm
#'
#' @examples
#' \donttest{
#' ##First, simulate covariate data
#' #Set the model parameters for covariates
#' covariate_names = c("x1", "x2", "x3") 
#' means = c(10,20,30)
#' names(means) = covariate_names
#' correlation_matrix = matrix(c( 1.0, 0.5, 0.3, 
#'				 0.5, 1.0, 0.4, 
#'				 0.3, 0.4, 1), nrow = 3
#'                            , byrow = TRUE)
#' rownames(correlation_matrix) = covariate_names
#' colnames(correlation_matrix) = covariate_names
#' covariateSDs = c(1,2,4)
#' names(covariateSDs) = covariate_names
#' n=1000
#' X = simulateCorrelatedCovariateData(n=n, covariateNames=covariate_names, means=means, R=correlation_matrix, covariateSDs=covariateSDs)
#'
#' ##Next, simulate response data based on simulated covariate data
#' #Set the model parameters
#' slopes = c(1,2,3)
#' names(slopes) = covariate_names
#' response_name = "y"
#' mean_response = 100
#' noise = 3
#' #Simulate the data
#' simulated_data = simulateResponseData(X=X, beta=slopes, 
#'                                                  response_name = response_name,
#'                                                  anchor=mean_response,
#'                                                  noise_sd = noise)
#' #Make a scatter plot of the simulated data
#' plot(simulated_data$x1, simulated_data$y, main=NULL, xlab="X1", ylab="Y")
#' #Present some stats to check match with true parameters
#' #Coefficients should match slopes for the covariates
#' summary(lm(simulated_data$y ~ simulated_data$x1 + simulated_data$x2 + simulated_data$x3))
#' #Correlation matrix should be close to the true correlation matrix
#' cor(X)
#' #Means should resemble the provided ones for the response and covariates
#' colMeans(simulated_data)
#' }
#' @noRd
simulateResponseData = function(X, beta, response_name = "Y", anchor = 0, noise_sd = 0) {
# Basic checks
	if (!is.data.frame(X) && !is.matrix(X)) {
	  stop("X must be a data.frame or matrix.")
	}

	X <- as.data.frame(X)

	n = nrow(X)

	xn <- colnames(X)
	bn <- names(beta)

	if (is.null(xn)) stop("X must have column names.")
	if (is.null(bn)) stop("beta must be a named vector.")

  # Check name matching
	missing <- setdiff(xn, bn)
	extra   <- setdiff(bn, xn)

	if (length(missing) > 0) {
	  stop("beta is missing coefficients for: ", paste(missing, collapse = ", "))
	}
	if (length(extra) > 0) {
	  warning("beta has extra coefficients not in X: ", paste(extra, collapse = ", "))
	}

  # Align beta to X columns
	beta <- beta[xn]

  # Linear predictor
	linear_part <- as.vector(as.matrix(X) %*% beta)

  # Adjust intercept for anchor
	intercept <- anchor - mean(linear_part)

  # Simulate noise
	epsilon <- rnorm(n, mean = 0, sd = noise_sd)

  # Simulate response
	Y <- intercept + linear_part + epsilon

  # Combine into dataframe
	df <- data.frame(Y = Y, X, check.names = FALSE)
	names(df)[1] <- response_name

	return(df)
}

#' Simulate a correlation matrix
#'
#' @description Simulates a correlation matrix, coviariate variances, and provides the covariance matrix as well.
#' 
#' @param C The number of covariates (default: 3)
#' @param K The number of hidden factors that influence the correlations of the covariates (default: 2)
#' @param loading_sd Standard deviation of the strength of covariation. Defines the average difference from the mean covariance, which is 0 (default: 0.7)
#' @param var_min The minimum variance of a covariate, which is uniformly sampled between var_min and var_max. Must be positive and less than var_max. (default: 0.3)
#' @param var_max The maximum variance of a covariate, which is uniformly sampled between var_min and var_max. Must be positive and more than var_min. (default: 1)
#'
#' @return A list with the correlation matrix, R, the covariance matrix, Sigma, the loading matrix Lambda defined by the normally distributed loading, the diagonal matrix of variances, Psi, covariateSDs, where Sigma is Lambda %*% t(Lambda) + Psi
#'
#' @examples
#' \donttest{
#' nCovariates = 5
#' nLatentFactorsUnderlyingCorrelations = 2
#' averageDeviationFromMeanCovariance = 0.7
#' minimumVariance = 0.3
#' maximumVariance = 1.0
#' simulateFactorCorrelation(nCovariates, nLatentFactorsUnderlyingCorrelations, averageDeviationFromMeanCovariance, minimumVariance, maximumVariance)
#' }
#' @noRd
simulateFactorCorrelation <- function(C = 3, 
				      K = 2, 
				      loading_sd = 0.7, 
				      var_min = 0.3, 
				      var_max = 1.0) {

	stopifnot(C >= 1)
	stopifnot(K >= 1)
	
	# If C = 1, K must effectively be 1 (otherwise Lambda is 1xK, still OK,
	# but the result is just a scalar variance anyway).
	K_use <- min(K, C)
	
	# Loadings matrix: C x K_use
	Lambda <- matrix(rnorm(C * K_use, 0, loading_sd), nrow = C, ncol = K_use)
	
	# Unique variances (length C)
	psi <- runif(C, var_min, var_max)
	
	# Make Psi safely a C x C matrix
	Psi <- diag(as.numeric(psi), nrow = C, ncol = C)
	
	# Covariance
	Sigma <- Lambda %*% t(Lambda) + Psi
	
	# Correlation: for P=1, cov2cor returns 1
	R <- if (C == 1) matrix(1, 1, 1) else cov2cor(Sigma)
	
	return(list(R = R, Sigma = Sigma, Lambda = Lambda, Psi = Psi, covariateSDs = sqrt(psi)))
}

#' Simulate a slope value for a phenological stage response to covariates
#'
#' @description Simulates a slope value for a phenological stage response to covariates
#' 
#' @param windowBelow The mean duration of the stage before the current one. (default: 10)
#' @param windowAbove The mean duration of the current stage. (default: 10)
#' @param minCovariate The minimum value of the covariate. (default: -10)
#' @param maxCovariate The maximum value of the covariate. (default: 10)
#' @param stageMinimumSeparation A value describing the minimum amount of separation between the previous and the subsequent stage. (default: 10)
#'
#' @return A slope value giving rise to data that fit, on average, within the input constraints
#'
#' @examples
#' \donttest{
#' previousStageMeanDuration = 50
#' currentStageMeanDuration = 60
#' minCovariate = -4
#' maxCovariate = 27
#' stageMinimumSeparation = 5
#' slope = simulateCovariateSlope(previousStageMeanDuration, currentStageMeanDuration, minCovariate, maxCovariate, stageMinimumSeparation=10)
#' n = 100
#' noise = 3
#' x = runif(n, minCovariate,maxCovariate)
#' y = rnorm(n, previousStageMeanDuration + x * slope, noise) 
#' plot(x,y)
#' }
#' @noRd
simulateCovariateSlope = function(windowBelow=10, windowAbove=10, minCovariate=-10, maxCovariate=10, stageMinimumSeparation=10) {
	#print("below")
	#print(windowBelow)
	#print("above")
	#print(windowAbove)
	#print("min sep")
	#print(stageMinimumSeparation)
	maxD = min(windowBelow, windowAbove) - stageMinimumSeparation
	if(maxD < 0 ) { return(0) }
	Bmax = maxD /(maxCovariate - minCovariate)
	return(runif(1,-Bmax,Bmax))
}

#nSDNormal number of standard deviations of the mean to use for normal distribution sampling
#nSDOverlap number of standard deviations of the mean onset for each stage used to determine region of overlap
resampleBiasedData = function(simulatedData=NULL, resample=FALSE, centeredNormal=FALSE, overlapOnly=FALSE, excludeLastStage=FALSE, replace=FALSE, newSampleSize=1, nSDNormal=1, nSDOverlap=3, minResponse=0, maxResponse=365) {
	#simulatedDataOrig = simulatedData
	if(is.null(simulatedData)) {
		stop("Please provide the output from simulateMultistageData and a bias type as input")
	}
	if(minResponse != 0) {
		stop("minResponse must be kept at 0.")
	}
	if(maxResponse <= 0) {
		stop("maxResponse must be greater than 0.")
	}
	if(newSampleSize < 1 && resample) {
		stop("The new sample size must be at least 1.")
	}
	if( (nSDNormal <= 0 && centeredNormal) || (nSDOverlap <= 0 && overlapOnly) ) {
		stop("The number of standard deviation units must be positive.")
	}
	if(!simulatedData$nonCyclical) {
		stop("Please renumber your stages so that the times between 0 and your first stage are labeled 1, stage 1 times are labeled 2, and so forth so that after the start of the last stage up to the end of the time period is nStages + 1. This function will still treat times as cyclical, especially if centeredNormal is set to TRUE or excludeLastStage is set to TRUE. In both cases, the samples from the last stage wrap around the 'start' point.")
	}

	nOrig = length(simulatedData$outputData$sampledTime)
	if(nOrig < newSampleSize && !replace && resample) {
		print(nOrig)
		print(newSampleSize)
		stop("The new sample size must be smaller than or equal to the original sample size.")
	}

	nStages = simulatedData$nStages+1	#Non cyclical, so that times between 0 and the first stage are labeled 1...
	nCovariates = simulatedData$nCovariates

	meanX = colMeans(simulatedData$outputData[,paste0("cov", 1:simulatedData$nCovariates)])		#estimate covariate means from data

	onsets = numeric(nStages)
	sigmas = numeric(nStages)
	slopes = matrix(rep(0, nStages*nCovariates), nrow=nStages, ncol=nCovariates)
	for(i in 1:nStages) {
		if(i==1) {
			onsets[i] = 0
			sigmas[i] = 0
			slopes[i,] = rep(0,nCovariates)
		}
		else if(i==2) {
			onsets[i] = simulatedData$stage1OnsetMean
			sigmas[i] = simulatedData$stage1OnsetSD
			slopes[i,] = simulatedData$stage1OnsetCovariateSlopes
		}
		else {
			onsets[i] = onsets[i-1] + simulatedData$stageDurationMeans[i-2]
			sigmas[i] = simulatedData$stageDurationSDs[i-2]
			slopes[i,] = simulatedData$stageDurationCovariateSlopes[i-2]
		}
		print(onsets[i])
	}

	if(centeredNormal) {	#Normally distributed bias centered over middle stage - resample times, but stages still based on previous onset times for the individual
		print("Centered normal sampling")
		nOrig = length(simulatedData$outputData$sampledTime)		#redundant, but safer if code ever reorganized
		meanInd = nStages/2					#determine where to center the normal distribution
		mean = (onsets[meanInd] + onsets[meanInd+1]) / 2
		sd = (onsets[nStages] - onsets[2]) / (2*nSDNormal) 	#start of "second" stage and start of "last" stage are mean +/- nSD
		print(mean)
		print(sd)
		simulatedData$outputData$sampledTime = rnorm(n=nOrig, mean=mean, sd=sd)

		#Wrap times that looped over the start of the cycle - assumes not more than 1 looping, which should be very rare, since the SDs of the normal distributions are much smaller than the range of response times 
		simulatedData$outputData$sampledTime[simulatedData$outputData$sampledTime<0] = maxResponse + simulatedData$outputData$sampledTime[simulatedData$outputData$sampledTime<0] 
		simulatedData$outputData$sampledTime = simulatedData$outputData$sampledTime %% maxResponse

		for(i in 1:nOrig) {
			for(s in 1:(nStages)) {
				if(s==1) {
					if(simulatedData$outputData$sampledTime[i]>=0 && 
					   simulatedData$outputData$sampledTime[i]<=simulatedData$outputData[i, nCovariates+s]) {
						simulatedData$outputData$sampledStage[i]=s
					}
				}
				else if(s<nStages) {
					if(simulatedData$outputData$sampledTime[i]>=simulatedData$outputData[i,nCovariates+s-1] && 	#output doesn't store 0 as first stage
					   simulatedData$outputData$sampledTime[i]<=simulatedData$outputData[i, nCovariates+s]) {
						simulatedData$outputData$sampledStage[i]=s
					}
				}
				else if(s==nStages) {
					if(simulatedData$outputData$sampledTime[i]>=simulatedData$outputData[i,nCovariates+s-1] && 
					   simulatedData$outputData$sampledTime[i]<=maxResponse) {
						simulatedData$outputData$sampledStage[i]=s
					}
				}
			}
		}
	}

	if(excludeLastStage) {
		print("Excluding last stage")
		simulatedData$outputData = simulatedData$outputData[ 
								    (simulatedData$outputData$sampledStage != nStages &
								     simulatedData$outputData$sampledStage != 1 ),]
		print(length(simulatedData$outputData$sampledTime))
	}


	stageCenter = numeric(nStages)
	stageLowerOverlap = numeric(nStages)
	stageUpperOverlap = numeric(nStages)
	if(overlapOnly) {	#Overlapping onset / cessation regions only

		nOO = nrow(simulatedData$outputData)

		position = rep(0,nOO)
		simulatedData$outputData$position = position

		temp <- as.data.frame(matrix(numeric(0), ncol = ncol(simulatedData$outputData)))
		colnames(temp) = colnames(simulatedData$outputData)

		print("Sampling only overlapping regions")
		print(nOO)
		for(i in 1:nStages) {
			if(i == 1) {
				#shiftI = minResponse	#doesn't matter
			}
			else {
				shiftI = onsets[i] +  (as.matrix(simulatedData$outputData[,1:nCovariates]) - meanX) %*% as.numeric(slopes[i,])
			}
			if(i == nStages) {
				#shiftII = maxResponse	#doesn't matter
			}
			else {
				shiftII = onsets[i+1] + (as.matrix(simulatedData$outputData[,1:nCovariates]) - meanX) %*% as.numeric(slopes[i+1,])
			}
			for(j in 1:nOO) {		##assumes that only two stages are overlapping at an average onset of a stage
				if(i==1) {
					ll = minResponse
					lu = minResponse
					ul = shiftII[j] - nSDOverlap * sigmas[i+1]
					uu = shiftII[j] + nSDOverlap * sigmas[i+1]
				}
				else if(i>1 && i<nStages) {
					ll = shiftI[j] - nSDOverlap * sigmas[i]
					lu = shiftI[j] + nSDOverlap * sigmas[i]
					ul = shiftII[j] - nSDOverlap * sigmas[i+1]
					uu = shiftII[j] + nSDOverlap * sigmas[i+1]
				}
				else if(i==nStages) {
					ll = shiftI[j] - nSDOverlap * sigmas[i]
					lu = shiftI[j] + nSDOverlap * sigmas[i]
					ul = maxResponse
					uu = maxResponse
				}

				t = simulatedData$outputData$sampledTime[j] 
				if(simulatedData$outputData$sampledStage[j] == i)  {
					if( t >= ll && t < lu) {
						simulatedData$outputData$position[j] = 1
						temp[nrow(temp) + 1, ] = simulatedData$outputData[j,]
						stageLowerOverlap[i] = stageLowerOverlap[i] + 1
					}
					else if( t > ul && t <= uu) {
						simulatedData$outputData$position[j] = 2
						temp[nrow(temp) + 1, ] = simulatedData$outputData[j,]
						stageUpperOverlap[i] = stageUpperOverlap[i] + 1
					}
					else if(t >= lu && t <= ul) {	
						simulatedData$outputData$position[j] = 3
						#center - do nothing - NOTHING IS ADDED TO TEMP!!!
					}
				}
			}
		#print(paste0("overlapping count: ", stageLowerOverlap[i] + stageUpperOverlap[i], " stage ", i))
		}

		#print("before unique")
		#print(nrow(temp))
		temp = unique(temp)
		#print(temp)
		simulatedData$outputData = temp
		#print("after assignment")
		#print(simulatedData$outputData)
		#print("after unique")
		#print(length(simulatedData$outputData$sampledTime))
	}


	if(resample) {	#Resampled original data
		nOrig = length(simulatedData$outputData$sampledTime)
		cat(paste("Resampling from ", nOrig, " original individuals.\n"))
		if(nOrig < newSampleSize) {
			fact = 2*newSampleSize/nOrig
			stop(paste("The new sample size must be smaller than or equal to the original sample size. Please increase the original sample size by at least ", fact, " times more."))
		}
		inds = sample(1:nOrig, newSampleSize, replace=FALSE)
		simulatedData$outputData = simulatedData$outputData[inds,]
		print(length(simulatedData$outputData$sampledTime))
	}


	#get the final breakdown of sample sizes in zones of overlap and between zones of overlap - this is inefficient because calculations are redone but conceptually easier
	stageCenter = numeric(nStages)
	stageLowerOverlap = numeric(nStages)
	stageUpperOverlap = numeric(nStages)

	nFinal = nrow(simulatedData$outputData)
	print("Partitioning overlap / central samples")
	print(nFinal)
	for(i in 1:(nStages)) {
		#shiftI = numeric(nFinal)
		#shiftII = numeric(nFinal)
		if(i == 1) {
			#shiftI = minResponse	#doesn't matter
		}
		else {
			shiftI = onsets[i] +  (as.matrix(simulatedData$outputData[,1:nCovariates]) - meanX) %*% as.numeric(slopes[i,])
		}
		if(i == nStages) {
			#shiftII = maxResponse	#doesn't matter
		}
		else {
			shiftII = onsets[i+1] + (as.matrix(simulatedData$outputData[,1:nCovariates]) - meanX) %*% as.numeric(slopes[i+1,])
		}
		for(j in 1:nFinal) {		##assumes that only two stages are overlapping at an average onset of a stage
			if(i==1) {
				ll = minResponse
				lu = minResponse
				ul = shiftII[j] - nSDOverlap * sigmas[i+1]
				uu = shiftII[j] + nSDOverlap * sigmas[i+1]
			}
			else if(i>1 && i<nStages) {
				ll = shiftI[j] - nSDOverlap * sigmas[i]
				lu = shiftI[j] + nSDOverlap * sigmas[i]
				ul = shiftII[j] - nSDOverlap * sigmas[i+1]
				uu = shiftII[j] + nSDOverlap * sigmas[i+1]
			}
			else if(i==nStages) {
				ll = shiftI[j] - nSDOverlap * sigmas[i]
				lu = shiftI[j] + nSDOverlap * sigmas[i]
				ul = maxResponse
				uu = maxResponse
			}
			t = simulatedData$outputData$sampledTime[j] 
			if(simulatedData$outputData$sampledStage[j] == i) {  
				if(!is.null(simulatedData$outputData$position)) {
					if(simulatedData$outputData$position[j]==1) {
						stageLowerOverlap[i] = stageLowerOverlap[i] + 1
					}
					else if(simulatedData$outputData$position[j]==2) {
						stageUpperOverlap[i] = stageUpperOverlap[i] + 1
					}
					else if(simulatedData$outputData$position[j]==2) {
						stageCenter[i] = stageCenter[i] + 1
					}
				}

				else if(t >= ll && t < lu) {
					stageLowerOverlap[i] = stageLowerOverlap[i] + 1
				}
				else if(t > ul && t <= uu) {
					stageUpperOverlap[i] = stageUpperOverlap[i] + 1
				}
				else if(t >= lu && t <= ul) {
					stageCenter[i] = stageCenter[i] + 1
					#if(overlapOnly) {
						#print("in center")
						#print(t)
						#print("boundary lower")
						#print(lu)
						#print("boundary upper")
						#print(ul)
						#stop("should not be any center points.")
					#}
				}
			}
		}
		#print(paste0("overlapping count: ", stageLowerOverlap[i] + stageUpperOverlap[i], " stage ", i))
	}

	simulatedData$stageCenter = stageCenter
	simulatedData$stageLowerOverlap = stageLowerOverlap
	simulatedData$stageUpperOverlap = stageUpperOverlap

	return(simulatedData)
}

#' Simulate covariate data, phenological stage onset times and sampled times and stage data
#'
#' @description Simulate covariate data, phenological stage onset times and sampled times and stage data. Usable with both presence-only and multistage analyses. For presence-only analyses, remove all rows in the output data set (outputData) with sampled times (sampledTime) outside the stage of iterest.
#' 
#' @param n The sample size. (default: 1000)
#' @param nStages The number of stages to simulate. (default: 2)
#' @param nonCyclical Encode the times before the first stage onset as a different stage and the times after the last stage onset as the last stage only (no wraparound). (default: TRUE)
#' @param stageNames A vector of the names of stages. (default: NULL)
#' @param stage1OnsetMean The mean value of the first stage's onset. (default: NULL)
#' @param stage1OnsetSD The standard deviation of the first stage's onset. (default: NULL)
#' @param stage1OnsetCovariateSlopes A vector with the first stage's onset model coefficients. (default: NULL)
#' @param stageDurationMeans A vector with all but the last stage's mean durations. (default: NULL)
#' @param stageDurationSDs A vector with all but the last stage's duration standard deviations. (default: NULL)
#' @param stageDurationCovariateSlopes A matrix with the slopes for each coefficient of the duration model for each stage. Columns are the covariates, and rows are the stages. (default: NULL)
#' @param stageMinimumSeparation The minimum separation, on average, between stage onsets. (default: 10)
#' @param minStageVariance The minimum allowable variance in durations for a stage. (default: 0)
#' @param maxStageVariance The maximum allowable variance in durations for a stage. (default: 3.0)
#' @param meanOnsetSpread A value or a vector of values (one for each stage, including before stage one and after the last stage), indicating the relative duration, on average, of that stage. Positive numbers only. Smaller numbers result in greater variation of durations. (default: 5.0)
#' @param minResponse The minimum response time. (default: 0)
#' @param maxResponse The maximum response time. (default: 365)
#' @param nCovariates The number of covariates. (default: 1)
#' @param X An n X nCovariates matrix with the covariate data. (default: NULL)
#' @param covariateMeans A vector of covariate means. (default: NULL)
#' @param covariateSDs A vector of covariate standard deviations. (default: NULL)
#' @param covariateNames A vector of covariate names. (default: NULL)
#' @param R A correlation matrix for the covariates. Used if no Sigma covariance matrix is provided. (default: NULL)
#' @param Sigma A covariance matrix for the covariates. (default: NULL)
#' @param nHiddenFactors The number of hidden factors that link to the covariates, if their correlation structure is to be simulated. (default: 2)
#' @param hiddenFactorStrength The average deviation from the mean covariance of 0. (default: 0.7)
#' @param minCovariateVariance The minimum allowable variance of a covariate. (default: 1.0)
#' @param maxCovariateVariance The maximum allowable variance of a covariate. (default: 25.0)
#' @param minCovariateMean The minimum mean value of a covariate (default: -10.0)
#' @param maxCovariateMean The maximum mean value of a covariate (default: 10.0)
#' @param seed The random number generator seed. (default: NULL)
#'
#' @return A list with the following: 
#'
#'		Data frame of sampled data for simulated individuals (outputData): covariate data (X), stage onsets (by name), sampled times (sampledTime), sampled stages (sampledStage)
#'
#'		Parameters for simulation: 
#'
#'			Covariate correlation structure: R, Sigma 
#'
#'			Covariate descriptive statistics: covariateSDs, covariateMeans
#'
#'			Covariate data: X
#'
#'			Stage 1 model: stage1OnsetMean, stage1OnsetSD, stage1OnsetCovariateSlopes
#'
#'			Duration models for all but the last stage: stageDurationMeans, stageDurationSDs, stageDurationCovariateSlopes
#'
#' @export
#'
#' @examples
#' \donttest{
#' ##Set basic simulation parameters 
#' numberStages = 9
#' numberCovariates = 3
#' sampleSize = 50
#'
#' ##Simulate data
#' simulatedData = simulateMultistageData(n=sampleSize, nStages=numberStages, nCovariates=numberCovariates)
#'
#' ##Plot simulated data with the first covariate along the x-axis.
#' #	set colors for the stages
#' #stageColors = rainbow(numberStages) 
#' #stageColors = viridisLite::viridis(numberStages)
#' stageColors = c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499") #Paul Tol
#' #	create the plot
#' plotMultistageSimulation(simulatedData=simulatedData, targetCovariateIndex=1, stageColors=stageColors)
#' }
simulateMultistageData = function(n=1000, 
				  nStages=2, 
				  nonCyclical=TRUE,
				  stageNames=NULL,
				  stage1OnsetMean=NULL,
				  stage1OnsetSD=NULL,
				  stage1OnsetCovariateSlopes=NULL,
				  stageDurationMeans=NULL,
				  stageDurationSDs=NULL,		
				  stageDurationCovariateSlopes=NULL,
				  stageMinimumSeparation=10,		#used to determine how much time between stages, at a minimum (on average, without including noise)
				  minStageVariance=16.0,
				  maxStageVariance=100.0,
				  meanOnsetSpread=5.0,		#dirichlet alphas, integer or vector of length nStages, larger values mean more even spread of durations, if not given
				  minResponse=0.0,
				  maxResponse=365.0, 		
				  nCovariates=1, 
				  X=NULL,			#allow users to enter the covariate information already
				  covariateMeans=NULL, 
				  covariateSDs=NULL, 
				  covariateNames=NULL,
				  R=NULL, 
				  Sigma=NULL, 
				  nHiddenFactors=2, 
				  hiddenFactorStrength=0.7, 
				  minCovariateVariance=16.0,
				  maxCovariateVariance=100.0,
				  minCovariateMean=-5.0,
				  maxCovariateMean=5.0,
				  seed=NULL) {

  if(!nonCyclical) {
    stop("Currently, stages must be coded so that they are \"unrolled\" and numbered from 1 to nStages+1. So, nonCyclical must be set to TRUE.")
  }

	#if(nonCyclical) {
		#individuals = as.data.frame( matrix(NA_real_, nrow = n, ncol = nStages+1))
		#colnames(individuals) = paste0("stage", 1:(nStages+1))
	#}
	#else {
		#individuals = as.data.frame( matrix(NA_real_, nrow = n, ncol = nStages))
		#colnames(individuals) = paste0("stage", 1:nStages)
	#}

	#Basic checks
	if(minResponse!=0) {
		stop("Minimum responses other than 0 are not supported.")
	}

	if(nStages<2) {
		stop("There must be at least 2 stages, even if these are just 'active' and 'dormant'.")
	}

	if(nCovariates<1) {
		stop("There must be at least 1 covariate.")
	}

	#Generate covariate names, if needed
	if(is.null(covariateNames) || length(covariateNames) < nCovariates) {
		covariateNames = paste0("cov", 1:nCovariates)
	}

	#Generate stage names, if needed
	if(is.null(stageNames) || length(stageNames) < nStages) {
		stageNames = paste0("stage", 1:nStages)
	}

	#Simulate covariate data
	#	check if simulation is needed
	simulate=FALSE
	if(is.null(X)) { simulate=TRUE }
	if(!simulate) { 
		if(nrow(X) != n) {
			warning("Simulating covariate data since the requested sample size does not match the given covariate data.")
			simulate=TRUE
		}
	}

	if(simulate) {
		#simulate correlation matrix if needed
		if((is.null(R) || is.null(covariateSDs)) && is.null(Sigma)) {
			RInfo=simulateFactorCorrelation(C=nCovariates, K = nHiddenFactors, loading_sd = hiddenFactorStrength, var_min = minCovariateVariance, var_max = maxCovariateVariance) 
			Sigma = RInfo$Sigma
			colnames(Sigma) = covariateNames
			rownames(Sigma) = covariateNames
			covariateSDs = RInfo$covariateSDs
			names(covariateSDs) = covariateNames
		}
	
		#simulate mean covariate values if needed
		if(is.null(covariateMeans)) {
			if(minCovariateMean>=maxCovariateMean) {
				stop("The minimum covariate mean must be smaller than the maximum covariate mean, or provide the means themselves.")
			}
			covariateMeans = runif(nCovariates,minCovariateMean,maxCovariateMean)
			names(covariateMeans) = covariateNames
		}

		#simulate covariate data - is a data frame
		X = simulateCorrelatedCovariateData(n=n, covariateNames=covariateNames, covariateMeans=covariateMeans, Sigma=Sigma, R=R, covariateSDs=covariateSDs, seed=seed) 
		
	}

	#simulate response data
	#	simulate stage 1 onset SD, if needed
	if(is.null(stage1OnsetSD)) {
		if(minStageVariance<0 || minStageVariance >= maxStageVariance) {
			stop("The minimum stage variance must be positive and less than the maximum stage variance.")
		}
		stage1OnsetSD = sqrt(runif(1, minStageVariance, maxStageVariance))
	}
	#	simulate onset means for each stage, and convert these to durations, as needed
	if(is.null(stage1OnsetMean) || is.null(stageDurationMeans) || length(stageDurationMeans)!=nStages-1) {
		if(length(meanOnsetSpread)==1) {
			meanOnsetSpread = rep(meanOnsetSpread,nStages+1)
		}
		if(any(meanOnsetSpread<0)) {
			stop("Each mean onset spread must be positive.")
		}
		if(length(meanOnsetSpread)!=nStages+1) {
			stop("Provide a vector of length nStages+1 representing the weight to give each stage (including before stage one and after the last stage) when randomly sampling stage lengths. Larger values give more weight and (non-intuitively) less variability in possible spread.")
		}
		means = phenoCollectR:::rdirichlet(1,meanOnsetSpread) * (maxResponse - minResponse) + minResponse
		stage1OnsetMean = means[1]	
		stageDurationMeans = means[2:nStages]
	}
	#	simulate duration SDs, if needed
  if(!is.null(stageDurationSDs)) {
    warning("Setting all duration SDs to 0. Only the first full stage has intrinsic noise.")
		stageDurationSDs = rep(0,nStages-1)
  }
	if(is.null(stageDurationSDs) || length(stageDurationSDs)!=nStages-1) {
		#stageDurationSDs = sqrt(runif(nStages-1, minStageVariance, maxStageVariance))
		stageDurationSDs = rep(0,nStages-1)
	}

	#	simulate stage 1 slopes for onset model, if needed
	if(is.null(stage1OnsetCovariateSlopes) || length(stage1OnsetCovariateSlopes)!=nCovariates) {
		stage1OnsetCovariateSlopes = rep(0,nCovariates)
		d1 = stage1OnsetMean
		d2 = ifelse(stageDurationMeans[1]<stageMinimumSeparation*5,stageMinimumSeparation*5,stageDurationMeans[1])
		#print("Stage 1")
		for(i in 1:nCovariates) {
			#the current / next is different here, since this is the onset model, not the duration model
			stage1OnsetCovariateSlopes[i] = simulateCovariateSlope(windowBelow=d1, windowAbove=d2, minCovariate=min(X[,i]), maxCovariate=max(X[,i]), stageMinimumSeparation = stageMinimumSeparation) 
		}
	}
	#	simulate stage duration covariate slopes, if needed
	if(!is.null(stageDurationCovariateSlopes) && (ncol(stageDurationCovariateSlopes)!=nCovariates || nrow(stageDurationCovariateSlopes)!=nStages-1)) {
		stop("Duration Covariate Slopes are provided, but the dimensions do no match: should be nStages-1 X nCovariates")
	}
	if(is.null(stageDurationCovariateSlopes) || ncol(stageDurationCovariateSlopes)!=nCovariates || nrow(stageDurationCovariateSlopes)!=nStages-1) {
		#print(stageDurationCovariateSlopes)
		stageDurationCovariateSlopes = matrix(0,ncol=nCovariates,nrow=nStages-1)
		colnames(stageDurationCovariateSlopes) = covariateNames
		rownames(stageDurationCovariateSlopes) = stageNames[1:nStages-1]
		for(i in 1:nCovariates) {
			#print("covariate")
			#print(i)
			mC = min(X[,i])
			MC = max(X[,i])
			for(j in 1:(nStages-1)) {
				#print("stage")
				#print(j)
				windowBelow = stageDurationMeans[j]
				if(j == nStages-1) { #wrap around
					windowAbove = stage1OnsetMean
				}
				else {
					windowAbove = stageDurationMeans[j+1]
				}
				stageDurationCovariateSlopes[j,i] = simulateCovariateSlope(windowBelow=windowBelow, windowAbove=windowAbove, minCovariate=mC, maxCovariate=MC, stageMinimumSeparation = stageMinimumSeparation)
			}
		}
	}

	#Simulate individuals' stages
	#	Basic checks
	if(nrow(X) != n) {
		stop("The number of rows in the covariate data matrix does not match the requested sample size.")
	}
	if(ncol(X) != nCovariates) {
		stop("The number of columns in the covariate data matrix does not match the requested number of covariates.")
	}

	stages = rep(0,n)
	onsets = matrix(0,nrow=n,ncol=nStages)
	colnames(onsets) = stageNames
	#onsetCumTot = 0
	for(i in 1:n) {
		#covariates = as.vector(X[i,]) #get vector of covariate values
		covariates = as.numeric(X[i,]) #get vector of covariate values
	  #onsetCumTot = 0
		for(j in 1:nStages) {
			if(j==1) {
				slopes = stage1OnsetCovariateSlopes
				#print(slopes)
				#print(covariateMeans)
				intercept = stage1OnsetMean - sum(as.numeric(slopes) * as.numeric(covariateMeans))
				sd = stage1OnsetSD
				#onsets[i,j] = rnorm(1, softplus(intercept + sum(slopes * covariates)),sd)
				onsets[i,j] = rnorm(1, intercept + sum(slopes * covariates),sd) #negative possible for first onset
				#if(times[i]>0 && times[i]<=onsets[i,j]) {
					#stages[i] = nStages
				#}
				#onsetCumTot = onsets[i,j]
			}
			else {
				slopes = as.vector(stageDurationCovariateSlopes[j-1,])
				intercept = stageDurationMeans[j-1] - sum(as.numeric(slopes) * as.numeric(covariateMeans))
				#sd = stageDurationSDs[j-1]    #Should be 0
        sd = 0  #Make it official
				#print(slopes)
				#print(covariateMeans)
				duration = softplus(rnorm(1,intercept + sum(slopes * covariates),sd)) #Softplus ok here. Sd is 0.
				#duration = rnorm(1,intercept + sum(slopes * covariates),sd)
				#onsets[i,j] = onsetCumTot + duration
        onsets[i,j] = onsets[i,j-1] + duration
				#diff = onsets[i,j] - onsets[i,j-1]  #should never happen due to softplus
				#if(diff < 0) {
					#onsets[i,j] = onsets[i,j-1]+1e-5 #stage of 1e-5 duration
				#}
				#if(times[i]>onsets[i,j-1] && times[i]<=onsets[i,j]) {
					#stages[i] = j-1
				#}
				#onsetCumTot = onsetCumTot + stageDurationMeans[j-1]
				#onsetCumTot = onsetCumTot + onsets[i,j]
			}
		}
		#if(times[i]>=onsetCumTot) {
			#stages[i]=nStages
		#}
	}
	#
	#simulate the onset times for each stage and simulate sampling each individual in the population
	times = runif(n,minResponse,maxResponse)
  #minO = min(onsets)
  #maxO = max(onsets)
	#times = runif(n,ifelse(minResponse<minO,minResponse,minO),ifelse(maxResponse>maxO, maxResponse,maxO))

	for(i in 1:n) {
		for(j in 1:nStages) {
			if(nonCyclical) {
				if(j==1) {
					#if(times[i]>0 && times[i]<=onsets[i,j]) {
					if(times[i]<=onsets[i,j]) {   #allow negative
						stages[i] = 1
					}
				}
				else {
					if(times[i]>onsets[i,j-1] && times[i]<=onsets[i,j]) {
						stages[i] = j
					}
				}
				if(j==nStages && times[i]>=onsets[i,j]) {   #allow times past maxResponse without wraparound (nonCyclical=TRUE)
					stages[i] = nStages+1
				}
			}
			else {
        stop("nonCyclical must be set to TRUE.")
				if(j==1) {
					if(times[i]>0 && times[i]<=onsets[i,j]) {
						stages[i] = nStages
					}
				}
				else {
					if(times[i]>onsets[i,j-1] && times[i]<=onsets[i,j]) {
						stages[i] = j-1
					}
				}
				if(j==nStages && times[i]>=onsets[i,j]) {
					stages[i] = nStages
				}
			}
		}
	}

	outputData = cbind(X,as.data.frame(onsets))
	outputData$sampledTime = times
	outputData$sampledStage = stages

	return(list(outputData=outputData, 
		    nStages = nStages,
		    nCovariates = nCovariates,
		    nonCyclical = nonCyclical,
		    R=R,
		    Sigma=Sigma,
		    covariateSDs=covariateSDs,
		    covariateMeans=covariateMeans,
		    X=X,
		    stage1OnsetMean=stage1OnsetMean,
		    stage1OnsetSD=stage1OnsetSD,
		    stage1OnsetCovariateSlopes=stage1OnsetCovariateSlopes,
		    stageDurationMeans=stageDurationMeans,
		    stageDurationSDs=stageDurationSDs,
		    stageDurationCovariateSlopes=stageDurationCovariateSlopes))
}


#' Simulate covariate data and phenological unit stages nested within individuals
#'
#' @description Simulate covariate data, phenological unit stage counts, onset times, and sampled times. An example would be an individual plant with multiple flowers as the phenological unit, with flower units in possibly different phenological stages (bud, flower, immature fruit, mature fruit). The function also optionally models transitions from a pre-(phenological unit) stage, to stages with multiple phenological units, to post-(phenological unit) stages after those units have senesced.
#' 
#' @param n The sample size. (default: 500)
#' @param meanUnits The mean number of phenological units that an individual develops during a phenological time period (e.g., an individual plant produces multiple flower units during a phenological cycle). (default: 20)
#' @param nStages The number of stages to simulate. This optionally includes the 'pre-' and 'post-' stages (see parameter 'includePrePost').  (default: 4; one pre-, one post-, and 2 internal stages with multiple phenological units)
#' @param nCovariates The number of covariates to simulate. (default: 2)
#' @param includePrePost Include pre- and post- stages for which phenological units have not yet developed (pre-) or have senesced (post-). (default: TRUE)
#' @param ... Additional optional parameters passed to simulateMultistageData (called from within this simulateMultistageOverlapData function). 
#'
#' @return A list with the following: 
#'
#'		Data frame named 'simulatedData' with the output from the called simuateMultistageData function
#'    A list of 'simulatedIndividuals' with:
#'      The time when an individual was observed
#'      A vector, 'phenologicalUnitStages', of the stages of the phenological units within the individual. If the individual is pre- or post-, there will be a single number representing that stage.
#'      A vector, 'stageCounts', of the counts of phenological units within each of the nStages stages. If the individual is pre- or post-, there will be a 1 in the index position of that stage and zeroes in the other stages. Otherwise, there will be a count of the number of phenological units in the stage.
#'      A vector, 'pi', of the probabilities of being in each stage at the sampled time.
#'
#' @export
#'
#' @examples
#' \donttest{
#' ##Set basic simulation parameters 
#' numberStages = 5           #Since 'includePrePost' is set to TRUE below, this include pre-(phenological unit), phenological unit stage 1, phenological unit stage 2, and post-(phenological unit) stages. An example would be a plant with the following stages: the plant hasn't produced any flower units (pre-), the plant has flower units are open flowers (phenologial unit stage 1), the plant has flower units that have developed into fruits (phenological unit stage 2), and the plant has fruits that have senesced (post-). 
#' numberCovariates = 3
#' sampleSize = 500
#' meanUnits = 30             #Set the mean number of phenological units (e.g., flowers) that an individual develops during a phenological cycle.
#' includePrePost = TRUE      #The first and last stages will be pre- and post- (reproductive) stages, evaluated at the level of the individual (1 for in the stage, 0 for not in the stage)
#'
#' ##Simulate data
#' sim = simulateMultistageOverlapData(n=sampleSize, meanUnits=meanUnits, nStages=numberStages, nCovariates=numberCovariates, includePrePost=includePrePost)
#'
#' ##Plot simulated data with the first covariate along the x-axis and the predominant stage color coded. Note, this plot does not represent all the stages of all the phenological units, but rather a single 'x' for each indivdiual color-coded by the predominant stage.
#' #	set colors for the stages
#' #stageColors = rainbow(numberStages) 
#' #stageColors = viridisLite::viridis(numberStages)
#' stageColors = c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499") #Paul Tol
#' #	create the plot
#' plotMultistageSimulation(simulatedData=sim$simulatedData, targetCovariateIndex=1, stageColors=stageColors)
#' }
simulateMultistageOverlapData = function(n=500, meanUnits = 20, nStages=4, nCovariates=2, includePrePost=TRUE, ...) {
 
  if(includePrePost && nStages<=3) {
    stop("Provide at least 4 stages when simulating pre- and post- stages, so that there is a minimum of two 'middle' stages.")
  } 

  simulatedData = simulateMultistageData(n=n, nStages=nStages-1, nCovariates=nCovariates, ...)    #breaks one stage into pre-, post-, so nStages-1 results in nStages
  individuals = vector("list", n)
  nUnits = rpois(n, meanUnits)
 
  onsetCols = (nCovariates+1):(nCovariates+nStages-1) 
  for(i in 1:n) {
    onsets = as.numeric(simulatedData$outputData[i,onsetCols])
    #print("onsets from simulation")
    #print(onsets)
    #print("time")
    #print(simulatedData$outputData$sampledTime[i])
    pi = stageProbabilities(t = simulatedData$outputData$sampledTime[i], 
              onsets = onsets, 
              SD = simulatedData$stage1OnsetSD    #just one SD with these data
         )
    #print("pi")
    #print(pi)
    #if(includePrePost) {
      ##Determine if individual is in pre- or post- stage
      #stage = simulatedData$outputData$sampledStage[i]
      #units = stage
      #if(stage == 1) {
        #stageCounts = rep(0,nStages)
        #stageCounts[1] = 1
      #}
      #else if(stage == nStages) {
        #stageCounts = rep(0,nStages)
        #stageCounts[nStages] = 1
      #}
      #else {
        #piPP = pi[2:(nStages-1)] / sum(pi[2:(nStages-1)])
        #units = sample(2:(nStages-1), nUnits[i], replace=TRUE, prob=piPP)
        #stageCounts = tabulate(units, nbins=nStages)
      #}
    #}
    #else {
      units = sample(1:nStages, nUnits[i], replace=TRUE, prob=pi)
      stageCounts = tabulate(units, nbins=nStages)
    #}

    if(includePrePost) {
      if(stageCounts[1] == nUnits[i]) {
        stageCounts = rep(0,nStages)
        stageCounts[1] = 1
      }
      else if(stageCounts[nStages] == nUnits[i]) {
        stageCounts = rep(0,nStages)
        stageCounts[nStages] = 1
      }
      else {
        stageCounts[1] = 0
        stageCounts[nStages] = 0
      }
    }
    individuals[[i]] = list(
      sampledTime = simulatedData$outputData$sampledTime[i],
      phenologicalUnitStages = units,
      stageCounts = stageCounts,
      pi = pi
    )
  }
  return(list(simulatedIndividuals = individuals, simulatedData = simulatedData))
}

stageProbabilities = function(t, onsets, SD) {
#print(t)
#print(onsets)
#print(SDs)
  S = length(onsets)
#print("# stages")
#print(S)
  pi = numeric(S+1) # S is the number of onsets (not including baseline), so number of stages is S+1
  pi[1] = 1 - pnorm((t - onsets[1])/SD)         #before onset 1
  pi[S+1] = pnorm((t - onsets[S])/SD)             #after onset S

  for(i in 1:(S-1)) {                           
    pi[i+1] = pnorm((t - onsets[i])/SD) - pnorm((t - onsets[i+1])/SD)
  }
  pi = pmax(pi, 0)
  return(pi / sum(pi))
}

#' Simulate linearly correlated covariate data
#'
#' @description Simulate linearly correlated covariate data for any positive number of covariates. 
#'
#' @param n Sample size
#' @param covariateNames Vector of the names of the covariates
#' @param means Vector of the means of the covariates; must be of same length as beta. (default: all 0 mean)
#' @param Sigma The covariance matrix. Must be real positive semidefinite and of dimension length(beta) X length(beta) (default: identity matrix if no correlation matrix provided)
#' @param R The correlation matrix. Diagonal must be all ones, and off diagonal entries must vary between -1 and 1. Needs value of covariateSDs parameter to be applied. (default: identity matrix)
#' @param covariateSDs A vector of standard deviations of each covariate. 
#' @param seed An optional seed for the random number generator.
#'
#' @return A data frame with columns labeled with the covariate names and n rows representing simulation replicates.
#' @importFrom MASS mvrnorm
#'
#' @examples
#' \donttest{
#' #Set the model parameters
#' covariate_names = c("x1", "x2", "x3") 
#' means = c(10,20,30)
#' names(means) = covariate_names
#' correlation_matrix = matrix(c( 1.0, 0.5, 0.3, 
#'				 0.5, 1.0, 0.4, 
#'				 0.3, 0.4, 1), nrow = 3
#'                            , byrow = TRUE)
#' rownames(correlation_matrix) = covariate_names
#' colnames(correlation_matrix) = covariate_names
#' covariateSDs = c(1,2,4)
#' names(covariateSDs) = covariate_names
#' n=1000
#' #Simulate the data
#' X = simulateCorrelatedCovariateData(n=n, covariateNames=covariate_names, covariateMeans=means, R=correlation_matrix, covariateSDs=covariateSDs)
#' #Make a scatter plot of the simulated data
#' plot(X$x1, X$x2, main=NULL, xlab="X1", ylab="X2")
#' #Correlation matrix should be close to the true correlation matrix
#' cor(X)
#' }
#' @noRd
simulateCorrelatedCovariateData = function (n, covariateNames=NULL, covariateMeans=NULL, Sigma = NULL, R = NULL, covariateSDs = NULL, seed = NULL)  {

	if(is.null(covariateNames)) {
		return(rep(0,n))
	}

#Rename n
	N = n

#Set seed
	if (!is.null(seed)) { set.seed(seed) }

# Load required package
	if (!requireNamespace("MASS", quietly = TRUE)) {
		stop("Please install the 'MASS' package with install.packages('MASS')")
	}

  # ---- all covariates (union) ----
  all_covs = covariateNames
  p <- length(all_covs)

  # ---- checks ----
  if (is.null(Sigma) && (is.null(R) || is.null(covariateSDs))) {
    warn("No covariance or correlation matrix with scales provided. Using identity matrix for covariance matrix.")
	Sigma = diag(p)
  	rownames(Sigma) = all_covs
	colnames(Sigma) = all_covs
  }

  # ---- build Sigma_full ----
  if (!is.null(Sigma)) {
    if (!is.matrix(Sigma)) stop("Sigma must be a matrix.")

    rn <- rownames(Sigma)
    cn <- colnames(Sigma)
    if (is.null(rn) || is.null(cn)) {
      stop("Sigma must have rownames and colnames.")
    }
    if (!identical(rn, cn)) {
      stop("Sigma rownames and colnames must match and be in the same order.")
    }

    missing <- setdiff(all_covs, rn)
    if (length(missing) > 0) {
      stop("Sigma is missing covariates: ", paste(missing, collapse = ", "))
    }

    Sigma_use <- Sigma[all_covs, all_covs, drop = FALSE]

  } else {
    # R + covariateSDs path
    if (!is.matrix(R)) stop("R must be a matrix.")

    rn <- rownames(R)
    cn <- colnames(R)
    if (is.null(rn) || is.null(cn)) {
      stop("R must have rownames and colnames.")
    }
    if (!identical(rn, cn)) {
      stop("R rownames and colnames must match and be in the same order.")
    }

    missing <- setdiff(all_covs, rn)
    if (length(missing) > 0) {
      stop("R is missing covariates: ", paste(missing, collapse = ", "))
    }

    R_use <- R[all_covs, all_covs, drop = FALSE]

    # covariateSDs can be named or unnamed
    if (is.null(names(covariateSDs))) {
      if (length(covariateSDs) != nrow(R)) {
        stop("If covariateSDs is unnamed, it must have length equal to nrow(R).")
      }
      sds_named <- covariateSDs
      names(sds_named) <- rownames(R)
    } else {
      sds_named <- covariateSDs
    }

    missing_sds <- setdiff(all_covs, names(sds_named))
    if (length(missing_sds) > 0) {
      stop("covariateSDs is missing covariates: ", paste(missing_sds, collapse = ", "))
    }

    sds_use <- sds_named[all_covs]
    p <- length(sds_use)
    if(p==1) {
	D <- diag(as.numeric(sds_use), nrow = p, ncol = p)
    }
    else {
	    D = sds_use
    }
    Sigma_use <- diag(D) %*% R_use %*% diag(D)
    #print(Sigma_use)
    #print(R_use)
    #print(sds_use)
    #print(D)
    #print(diag(D))
  }

  # ---- means ----
  if (is.null(covariateMeans)) {
    means_use <- rep(0, p)
    names(means_use) <- all_covs
  } else {
    if (is.null(names(covariateMeans))) {
      stop("means must be a named vector with names matching covariates.")
    }
    missing_means <- setdiff(all_covs, names(covariateMeans))
    if (length(missing_means) > 0) {
      stop("means is missing covariates: ", paste(missing_means, collapse = ", "))
    }
    means_use <- covariateMeans[all_covs]
  }

  # ---- PD check ----
  ev <- eigen(Sigma_use, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-10) {
    stop("Sigma is not positive definite (min eigenvalue = ", signif(min(ev), 4), ").")
  }

  # ---- simulate ----
  X <- MASS::mvrnorm(n = N, mu = means_use, Sigma = Sigma_use)
  colnames(X) <- all_covs
  X <- as.data.frame(X)

return(X)
}

#' Simulate a response variable based on multiple correlated covariates
#'
#' @description Simulate the values of a response variable based on simulated values of multiple correlated covariates. Assumes a linear model with normally distributed variation.
#'
#' @param n Sample size
#' @param beta Vector of slope coefficients of the covariates
#' @param covariateNames Vector of the names of the covariates
#' @param mu Vector of the means of the covariates; must be of same length as beta. (default: all 0 mean)
#' @param response_name The name of the response variable (default: "Y")
#' @param Sigma The covariance matrix. Must be of dimension length(beta) X length(beta) (default: identity matrix). Either this or the correlation matrix R with covariate standardiviations (covariateSDs) should be provided, and if none, the identity matrix is default. Default: identity matrix
#' @param R The correlation matrix. Either with the standard deviations of each covariate (covariateSDs) or the Sigma covariance matrix should be given, and if none, the identity matrix is default. Default: identity matrix
#' @param anchor Marginal mean value of the response variable
#' @param noise_sd Standard deviation of the noise in the response variable
#'
#' @return A data frame with columns labeled with the variable names and rows representing simulation replicates.
#' @importFrom MASS mvrnorm
#'
#' @examples
#' \donttest{
#' #Set the model parameters
#' covariate_names = c("x1", "x2", "x3") 
#' slopes = c(1,2,3)
#' names(slopes) = covariate_names
#' means = c(10,20,30)
#' names(means) = covariate_names
#' response_name = "y"
#' covariance_matrix = matrix(c( 1.0, 0.5, 0.3, 0.5, 2.0, 0.4, 0.3, 0.4, 1.5), nrow = 3
#'                            , byrow = TRUE)
#' colnames(covariance_matrix) = covariate_names
#' rownames(covariance_matrix) = covariate_names
#' mean_response = 100
#' noise = 3
#' n=50000
#' #Simulate the data
#' simulated_data = simulateCorrelatedCovariateAndResponseData(n=n, beta=slopes, covariateNames=covariate_names
#'                                                  , cov_means = means, Sigma=covariance_matrix
#'                                                  , anchor=mean_response
#'                                                  , response_name = response_name
#'                                                  , noise_sd = noise)
#' #Make a scatter plot of the simulated data
#' plot(simulated_data$x1, simulated_data$y, main=NULL, xlab="X1", ylab="Y")
#' ##Compare inferences to true values
#' #coefficients should be 1, 2, 3
#' summary(lm(simulated_data$y ~ simulated_data$x1 + simulated_data$x2 + simulated_data$x3))
#' #Means should be 100, 10, 20, 30
#' colMeans(simulated_data)
#' }
#' @noRd
simulateCorrelatedCovariateAndResponseData = function(n, beta, covariateNames, cov_means = NULL, response_name = "Y", Sigma = NULL, R = NULL, covariateSDs = NULL, anchor = 0, noise_sd = 1) {

# Simulate covariates
	X = simulateCorrelatedCovariateData(n=n, covariateNames=cov_names, covariateMeans=cov_means, Sigma = Sigma, R = R, covariateSDs = covariateSDs)  

# Simulate response
	Y = simulateResponseData(X=X, beta=beta, response_name=response_name, anchor=anchor, noise_sd=noise_sd) 

#Return output
	return(Y)
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
			#print(c(alpha_s,beta_s,alpha_d,beta_d))
				return(list(error_m = "Infeasible parameter value provided as input under beta onset, beta duration model. Try different parameter values, paying close attention to the scale of sigma. Sigmas that are too large cannot be accommodated by the beta distribution.", error=T))

		}

		cnt = 1
		#t_start = numeric(N)
		#t_end = numeric(N)
		#duration_raw = numeric(N)
		#duration = numeric(N)
		#observed = numeric(N)
		Ts = rep(NA, N)

		t_start = rbeta(N, alpha_s, beta_s)
		duration_raw = rbeta(N, alpha_d, beta_d)
		duration = duration_raw * (1 - t_start)
		t_end =  t_start + duration
		observed = runif(N) 
		#sampled times of individuals in the phenophase is a subset of all the randomly observed times of the individuals (not all individuals are in the phenophase at a randomly observed time)
		condition = (observed>t_start & observed<t_end)
		Ts[condition] = observed[condition]

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
#' @description Simulate the phenophase onset, duration, and cessation times for individuals of a population. Phenological extremes (first onset, last cessation) are also provided. No covariates are included.
#'
#' Data can be simulated under either the beta onset, beta duration model (BB) or under the Gaussian process model (GP). 
#' 
#' Note that for the BB model, the number of observed times will be less than the number of individuals in the population because individuals are viewed at random times, and not all individuals are in the phenophase at the randomly observed time. Individuals with longer phenophases are more likely to be recorded. Only simulated specimens in the phenophase are recorded. Individuals that are not recorded occur as an NA value in the returned vector of observed collection times (Ts in the output list). This is not the case for the GP model, because under the GP model, every individual is equally likely to be recorded, so the observed time of each individual can be randomly sampled between the onset time and the cessation time without biasing the resulting distribution of observed collection times. 
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
#' \donttest{
#' #Set the parameters
#' N = 1000000
#' mu_O = 30
#' sigma_O = 10
#' mu_D_raw = 40
#' sigma_D = 14
#' #Simulate the data under the beta onset, beta duration (BB) model
#' data = simulatePopulation(N=N, mu_O=mu_O, sigma_O=sigma_O, mu_D_raw=mu_D_raw
#'                           , sigma_D=sigma_D, type="BB")
#' #Plot histograms of the phenological values
#' xlim = c(min(data$O, data$C), max(data$O, data$C))
#' breaks = seq(xlim[1],xlim[2], length.out=100)
#' hist(data$O, col=rgb(1,0,0,0.3), xlab="Day of year", probability=TRUE, breaks=breaks
#'      , xlim=xlim, main=NULL) #Onset
#' #Observed collection times
#' hist(na.omit(data$Ts), col=rgb(1,0,1,0.3), probability=TRUE, breaks=breaks, add=TRUE)
#' hist(data$C, col=rgb(0,0,1,0.3), probability=TRUE, breaks=breaks, add=TRUE) #Cessation
#' abline(v=data$Ok1, col="yellow") #First onset for the population as yellow vertical line
#' abline(v=data$CkN, col="cyan") #Last cessation for the population as cyan vertical line
#' #overlay theoretical density curve for onset
#' curve(dO(x,mu_O=mu_O, sigma_O=sigma_O, type="BB"), add=TRUE, col="red") 
#' #overlay theoretical density curve for observed times
#' curve(dT(x,mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D_raw, sigma_D=sigma_D, type="BB")
#'       , add=TRUE, col="purple") 
#' #overlay theoretical density curve for cessation
#' curve(dC(x,mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D_raw, sigma_D=sigma_D, type="BB")
#'       , add=TRUE, col="blue") 
#' #Use the above parameter values, but simulate data under the GP model (default)
#' #!!Note the "wrapping around" effect for which values below 0 get wrapped to the end 
#' #     of the "previous" time period 
#' data = simulatePopulation(N=N, mu_O=mu_O, sigma_O=sigma_O, mu_D_raw=mu_D_raw)
#' #Plot histograms of the phenological values
#' dev.new()
#' xlim = c(min(data$O, data$C), max(data$O, data$C))
#' breaks = seq(xlim[1],xlim[2], length.out=100)
#' hist(data$O, col=rgb(1,0,0,0.3), xlab="Day of year", probability=TRUE, breaks=breaks
#'      , xlim=xlim, main=NULL) #Onset
#' #Observed collection times
#' hist(data$Ts, col=rgb(1,0,1,0.3), probability=TRUE, breaks=breaks, add=TRUE)
#' #Cessation
#' hist(data$C, col=rgb(0,0,1,0.3), probability=TRUE, breaks=breaks, add=TRUE) 
#' #First onset for the population as yellow vertical line
#' abline(v=data$Ok1, col="yellow")
#' #Last cessation for the population as cyan vertical line
#' abline(v=data$CkN, col="cyan") 
#' #overlay theoretical density curve for onset
#' curve(dO(x,mu_O=mu_O, sigma_O=sigma_O), add=TRUE, col="red") 
#' #overlay theoretical density curve for observed times
#' curve(dT(x,mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D_raw, sigma_D=sigma_D)
#'       , add=TRUE, col="purple") 
#' #overlay theoretical density curve for cessation
#' curve(dC(x,mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D_raw, sigma_D=sigma_D)
#'       , add=TRUE, col="blue") 
#' }
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
