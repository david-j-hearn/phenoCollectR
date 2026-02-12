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
#' @export
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
	X = simulateCorrelatedCovariateData(n=n, cov_names=covs_all, means=covariateMeans, Sigma = covarianceMatrix, R = correlationMatrix, sds = covariateStandardDeviations, seed = seed)  
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
#' @export
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
#' sds = c(1,2,4)
#' names(sds) = covariate_names
#' n=1000
#' X = simulateCorrelatedCovariateData(n=n, cov_names=covariate_names, means=means, R=correlation_matrix, sds=sds)
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

#' Simulate linearly correlated covariate data
#'
#' @description Simulate linearly correlated covariate data for any positive number of covariates. 
#'
#' @param n Sample size
#' @param cov_names Vector of the names of the covariates
#' @param means Vector of the means of the covariates; must be of same length as beta. (default: all 0 mean)
#' @param Sigma The covariance matrix. Must be real positive semidefinite and of dimension length(beta) X length(beta) (default: identity matrix if no correlation matrix provided)
#' @param R The correlation matrix. Diagonal must be all ones, and off diagonal entries must vary between -1 and 1. Needs value of sds parameter to be applied. (default: identity matrix)
#' @param sds A vector of standard deviations of each covariate. 
#' @param seed An optional seed for the random number generator.
#'
#' @return A data frame with columns labeled with the covariate names and n rows representing simulation replicates.
#' @export
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
#' sds = c(1,2,4)
#' names(sds) = covariate_names
#' n=1000
#' #Simulate the data
#' X = simulateCorrelatedCovariateData(n=n, cov_names=covariate_names, means=means, R=correlation_matrix, sds=sds)
#' #Make a scatter plot of the simulated data
#' plot(X$x1, X$x2, main=NULL, xlab="X1", ylab="X2")
#' #Correlation matrix should be close to the true correlation matrix
#' cor(X)
#' }
simulateCorrelatedCovariateData = function (n, cov_names, means=NULL, Sigma = NULL, R = NULL, sds = NULL, seed = NULL)  {

	if(is.null(cov_names)) {
		return(rep(0,n))
	}

#Rename n
	N = n

#Set seed
	if (!is.null(seed)) set.seed(seed)

# Load required package
	if (!requireNamespace("MASS", quietly = TRUE)) {
		stop("Please install the 'MASS' package with install.packages('MASS')")
	}

  # ---- all covariates (union) ----
  all_covs = cov_names
  p <- length(all_covs)

  # ---- checks ----
  if (is.null(Sigma) && (is.null(R) || is.null(sds))) {
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
    # R + sds path
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

    # sds can be named or unnamed
    if (is.null(names(sds))) {
      if (length(sds) != nrow(R)) {
        stop("If sds is unnamed, it must have length equal to nrow(R).")
      }
      sds_named <- sds
      names(sds_named) <- rownames(R)
    } else {
      sds_named <- sds
    }

    missing_sds <- setdiff(all_covs, names(sds_named))
    if (length(missing_sds) > 0) {
      stop("sds is missing covariates: ", paste(missing_sds, collapse = ", "))
    }

    sds_use <- sds_named[all_covs]
    p <- length(sds_use)
	D <- diag(as.numeric(sds_use), nrow = p, ncol = p)
    Sigma_use <- diag(D) %*% R_use %*% diag(D)
  }

  # ---- means ----
  if (is.null(means)) {
    means_use <- rep(0, p)
    names(means_use) <- all_covs
  } else {
    if (is.null(names(means))) {
      stop("means must be a named vector with names matching covariates.")
    }
    missing_means <- setdiff(all_covs, names(means))
    if (length(missing_means) > 0) {
      stop("means is missing covariates: ", paste(missing_means, collapse = ", "))
    }
    means_use <- means[all_covs]
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
#' @param cov_names Vector of the names of the covariates
#' @param mu Vector of the means of the covariates; must be of same length as beta. (default: all 0 mean)
#' @param response_name The name of the response variable (default: "Y")
#' @param Sigma The covariance matrix. Must be of dimension length(beta) X length(beta) (default: identity matrix). Either this or the correlation matrix R with covariate standardiviations (sds) should be provided, and if none, the identity matrix is default. Default: identity matrix
#' @param R The correlation matrix. Either with the standard deviations of each covariate (sds) or the Sigma covariance matrix should be given, and if none, the identity matrix is default. Default: identity matrix
#' @param anchor Marginal mean value of the response variable
#' @param noise_sd Standard deviation of the noise in the response variable
#'
#' @return A data frame with columns labeled with the variable names and rows representing simulation replicates.
#' @export
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
#' simulated_data = simulateCorrelatedCovariateAndResponseData(n=n, beta=slopes, cov_names=covariate_names
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
simulateCorrelatedCovariateAndResponseData = function(n, beta, cov_names, cov_means = NULL, response_name = "Y", Sigma = NULL, R = NULL, sds = NULL, anchor = 0, noise_sd = 1) {

# Simulate covariates
	X = simulateCorrelatedCovariateData(n=n, cov_names=cov_names, means=cov_means, Sigma = Sigma, R = R, sds = sds)  

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
			print(c(alpha_s,beta_s,alpha_d,beta_d))
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
#' @description Simulate the phenophase onset, duration, and cessation times for individuals of a population. Phenological extremes (first onset, last cessation) are also provided. 
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

#' @importFrom copula getSigma
#' @importFrom MASS mvrnorm
#' @importFrom stats ecdf
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
