#' Functions that calculate the expected value of the associated phenological random variable:
#'
#' The random variables include the following:
#' O: phenological onset times
#' D: phenophase durations
#' C: phenological cessation times
#' Ok1: first onset (order k = 1)
#' CkN: last cessation (order k = N, the population size)
#' T: observed specimen collection times, e.g., day of year for a yearly time period
#' Most functions share a subset of the same input parameters, but if you use a function that does not require all the parameters, provide only the needed parameters.
#' The function name is created by prepending "E." in front of one of the above variables names. 
#'
#' @rdname expectation_functions
#' @param mu_O Mean onset time
#' @param sigma_O Standard deviation for the onset time distribution
#' @param mu_D Mean phenophase duration length
#' @param sigma_D Standard deviation for the phenophase duration distribution
#' @param minResponse Minimum value of the response (e.g., day of year); must be set to 0 under current implementation (default = 0)
#' @param maxResponse Maximum value of the response (e.g., day of year); typically 365 for Gregorian calendar (default = 365)
#' @param N The population size for estimation of extreme events
#' @param type The model type, either BB (beta onset, beta duration) or GP (Gaussian process with a shared standard deviation for onset and cessation and a constant duration) (default = "GP")
#' @param threshApprox An error threshold set to use approximation schemes when numerical integration fails. Safe to leave at default, and deprecated. (default: NA)
#'
#' @return The expected value of the associated random variable and input parameters
#' @examples
#' \donttest{
#' #Set the mean onset time (day of year)
#' mu_O = 100
#' #Set the standard deviation of the distribution of onset times
#' sigma_O = 7
#' #set the mean duration of the phenophase (days)
#' mu_D = 30
#' #Set the standard deviation of the distribution of durations
#' sigma_D = 7
#' #Set the population size
#' N = 100000
#'
#' #Calculate the expected cessation time for the beta onset, beta duration (BB) model
#' eC = E.C(mu_O=mu_O, mu_D=mu_D, type="BB")
#'
#' #Calculate the expected last cessation time under the BB model
#' eCkN = E.CkN(N=N, mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, type="BB")
#' }
#'
#' @name phenologyExpectationFunctions
NULL

#' @rdname expectation_functions
#' @export
E.C = function(mu_O, mu_D=NA, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="BB") {
		if(is.na(mu_D)) {
			stop("Mean duration must be provided as input.")
		}
		E.C.BB(mu_O=mu_O, mu_D=mu_D, minResponse=minResponse, maxResponse=maxResponse)
	}
	else if(type=="GP") {
		mu_C = mu_O + mu_D
		E.C.GP(mu_C=mu_C)
	}
	else {
		stop(paste("Unsupported type: ", type))
	}
}

#' @rdname expectation_functions
#' @export
E.CkN = function(N, mu_O, sigma_O, mu_D, sigma_D=NA, minResponse=0, maxResponse=365, type=c("GP","BB"), threshApprox=NA) {
	type = match.arg(type)
	if(type=="BB") {
		if(is.na(sigma_D)) {
			stop("The standard deviation for the duration must be provided.")
		}
		vals = E.CkN.BB(N=N, mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse, threshApprox = threshApprox)
	}
	else if(type=="GP") {
		mu_C = mu_O + mu_D
		sigma_C = sigma_O
		vals = E.CkN.GP(N=N, mu_C=mu_C, sigma=sigma_C, minResponse=minResponse, maxResponse=maxResponse, threshApprox = threshApprox)
	}
	else {
		stop(paste("Unsupported type: ", type))
	}
	return(vals)
}

#' @rdname expectation_functions
#' @export
E.D = function(mu_O, sigma_O, mu_D, sigma_D=NA, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="BB") {
		if(is.na(sigma_D)) {
			stop("The standard deviation for the duration must be provided.")
		}
		E.D.BB(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	}
	else if(type=="GP") {
		mu_C = mu_O + mu_D
		E.D.GP(mu_O=mu_O, mu_C=mu_C)
	}
	else {
		stop(paste("Unsupported type: ", type))
	}
}

#' @rdname expectation_functions
#' @export
E.O = function(mu_O, sigma_O, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="BB") {
		E.O.BB(mu_O=mu_O, sigma_O=sigma_O, minResponse=minResponse, maxResponse=maxResponse)
	}
	else if(type=="GP") {
		E.O.GP(mu_O=mu_O)
	}
	else {
		stop(paste("Unsupported type: ", type))
	}
}

#' @rdname expectation_functions
#' @export
E.Ok1 = function(N, mu_O, sigma_O, minResponse=0, maxResponse=365, type=c("GP","BB"), threshApprox=NA) {
	type = match.arg(type)
	if(type=="BB") {
		vals = E.Ok1.BB(N=N, mu_O=mu_O, sigma_O=sigma_O, minResponse=minResponse, maxResponse=maxResponse,threshApprox=threshApprox)

	}
	else if(type=="GP") {
		vals = E.Ok1.GP(N=N, mu_O=mu_O, sigma=sigma_O, minResponse=minResponse, maxResponse=maxResponse, threshApprox=threshApprox )
	}
	else {
		stop(paste("Unsupported type: ", type))
	}
	return(vals)
}

#' @rdname expectation_functions
#' @export
E.T = function(mu_O, mu_D, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="BB") {
		E.T.BB(mu_O=mu_O, mu_D=mu_D, minResponse=minResponse, maxResponse=maxResponse)
	}
	else if(type=="GP") {
		E.T.GP(mu_O=mu_O, mu_D=mu_D)
	}
	else {
		stop(paste("Unsupported type: ", type))
	}
}


#' Functions that calculate the frequentist probability interval (PI) for an associated phenological random variable:
#'
#' O: phenological onset times
#' D: phenophase durations
#' C: phenological cessation times
#' Ok1: first onset (order k = 1)
#' CkN: last cessation (order k = N, the population size)
#' T: observed specimen collection times, e.g., day of year for a yearly time period
#' PNt: Proportion of the population of size N in the phenophase at time t; only the unnormalized (for efficiency) 'd' and 'r' functions are currently implemented
#'
#' Most functions share a subset of the same input parameters, but if you use a function that does not require all the parameters, provide only the needed parameters.
#'
#' The function name is PI.* where * is one of the above random variables.
#'
#' @rdname pi_functions
#' @param mu_O Mean onset time
#' @param sigma_O Standard deviation for the onset time distribution
#' @param mu_D Mean phenophase duration length
#' @param sigma_D Standard deviation for the phenophase duration distribution
#' @param minResponse Minimum value of the response (e.g., day of year); must be set to 0 under current implementation (default = 0)
#' @param maxResponse Maximum value of the response (e.g., day of year); typically 365 for Gregorian calendar (default = 365)
#' @param N The population size for estimation of extreme events
#' @param type The model type, either BB (beta onset, beta duration) or GP (Gaussian process with a shared standard deviation for onset and cessation and a constant duration) (default = "GP")
#' @param alpha The alpha level. For example, alpha = 0.05 calculates the 95% probability interval
#'
#' @return A vector with the lower and upper values of the probability interval.
#' @examples
#' \donttest{
#' #Set the mean onset time
#' mean_onset = 100
#' #Set the onset time standard deviation, sigma
#' sigma_onset = 10
#' #Set the duration of the phenophase
#' duration = 50
#' #Calculate the 90% probability interval for the observed times under the GP model (default)
#' observed_t_PI = PI.T(mu_O = mean_onset, sigma_O=sigma_onset, mu_D=duration, alpha=0.1)
#' }
#'
#' @name phenologyProbabilityIntervalFunctions
NULL

#' @rdname pi_functions
#' @export
PI.C = function(mu_O, sigma_O, mu_D, sigma_D=NA, minResponse=0, maxResponse=365, alpha=0.05, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="BB") {
		if(is.na(sigma_D)) {
			stop("The standard deviation for the duration must be provided.")
		}
		PI.C.BB(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse, alpha=alpha)
	}
	else if(type=="GP") {
		mu_C = mu_O + mu_D
		PI.C.GP(mu_C=mu_C, sigma=sigma_O, alpha=alpha)
	}
	else {
		stop(paste("Unsupported type: ", type))
	}
}

#' @rdname pi_functions
#' @export
PI.CkN = function(N, mu_O, sigma_O, mu_D, sigma_D=NA, minResponse=0, maxResponse=365, alpha=0.05, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="BB") {
		if(is.na(sigma_D)) {
			stop("The standard deviation for the duration must be provided.")
		}
		PI.CkN.BB(N=N, mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse, alpha=alpha)
	}
	else if(type=="GP") {
		mu_C = mu_O + mu_D
		PI.CkN.GP(N=N, mu_C=mu_C, sigma=sigma_O, alpha=alpha)
	}
	else {
		stop(paste("Unsupported type: ", type))
	}
}

#' @rdname pi_functions
#' @export
PI.D = function(mu_O, sigma_O, mu_D, sigma_D=NA, minResponse=0, maxResponse=365, alpha=0.05, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="BB") {
		if(is.na(sigma_D)) {
			stop("The standard deviation for the duration must be provided.")
		}
		PI.D.BB(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse, alpha=alpha)
	}
	else if(type=="GP") {
		mu_C = mu_O + mu_D
		PI.D.GP(mu_O=mu_O, mu_C=mu_C)
	}
	else {
		stop(paste("Unsupported type: ", type))
	}
}

#' @rdname pi_functions
#' @export
PI.O = function(mu_O, sigma_O, minResponse=0, maxResponse=365, alpha=0.05, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="BB") {
		PI.O.BB(mu_O=mu_O, sigma_O=sigma_O, minResponse=minResponse, maxResponse=maxResponse, alpha=alpha)
	}
	else if(type=="GP") {
		PI.O.GP(mu_O=mu_O, sigma=sigma_O, alpha=alpha)
	}
	else {
		stop(paste("Unsupported type: ", type))
	}
}

#' @rdname pi_functions
#' @export
PI.Ok1 = function(N, mu_O, sigma_O, minResponse=0, maxResponse=365, alpha=0.05, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="BB") {
		PI.Ok1.BB(N=N, mu_O=mu_O, sigma_O=sigma_O, minResponse=minResponse, maxResponse=maxResponse, alpha=alpha)
	}
	else if(type=="GP") {
		PI.Ok1.GP(N=N, mu_O=mu_O, sigma=sigma_O, alpha=alpha)
	}
	else {
		stop(paste("Unsupported type: ", type))
	}
}

#' @rdname pi_functions
#' @export
PI.T = function(mu_O, sigma_O, mu_D, sigma_D=NA, minResponse=0, maxResponse=365, alpha=0.05, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="BB") {
		if(is.na(sigma_D)) {
			stop("The standard deviation for the duration must be provided.")
		}
		PI.T.BB(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse, alpha=alpha)
	}
	else if(type=="GP") {
		mu_C = mu_O + mu_D
		PI.T.GP(mu_O=mu_O, mu_C=mu_C, sigma=sigma_O, alpha=alpha)
	}
	else {
		stop(paste("Unsupported type: ", type))
	}
}

#' Functions that calculate the standard deviation (square root of variance) for an associated phenological random variable:
#'
#' O: phenological onset times
#' D: phenophase durations
#' C: phenological cessation times
#' Ok1: first onset (order k = 1)
#' CkN: last cessation (order k = N, the population size)
#' T: observed specimen collection times, e.g., day of year for a yearly time period
#' PNt: Proportion of the population of size N in the phenophase at time t; only the unnormalized (for efficiency) 'd' and 'r' functions are currently implemented
#'
#' Most functions share a subset of the same input parameters, but if you use a function that does not require all the parameters, provide only the needed parameters.
#'
#' The function name is SD.* where * is one of the above random variables.
#'
#' @rdname sd_functions
#' @param mu_O Mean onset time
#' @param sigma_O Standard deviation for the onset time distribution
#' @param mu_D Mean phenophase duration length
#' @param sigma_D Standard deviation for the phenophase duration distribution
#' @param minResponse Minimum value of the response (e.g., day of year); must be set to 0 under current implementation (default = 0)
#' @param maxResponse Maximum value of the response (e.g., day of year); typically 365 for Gregorian calendar (default = 365)
#' @param N The population size for estimation of extreme events
#' @param type The model type, either BB (beta onset, beta duration) or GP (Gaussian process with a shared standard deviation for onset and cessation and a constant duration) (default = "GP")
#' @param intFailLow Lowest acceptable integration value (optional; leave at default for most causes)
#' @param intFailHigh Highest acceptable integration value (optional; leave at default for most causes)
#'
#' @return The standard deviation of the associated random variable for the input parameter values.
#' @examples
#' \donttest{
#' #Set the mean onset time
#' mean_onset = 100
#' #Set the onset time standard deviation
#' sigma_onset = 10
#' #Set the mean duration of the phenophase
#' mean_duration = 50
#' #Set the duration standard deviation
#' sigma_duration = 7
#' #Set the population size
#' N=1000
#' #Calculate the standard deviation of the distribution of last cessation times under the beta 
#' #    onset, beta duration (BB) model
#' sdCkN = SD.CkN(N=N, mu_O=mean_onset, sigma_O=sigma_onset, mu_D=mean_duration
#'                , sigma_D=sigma_duration, type="BB")
#' }
#'
#' @name phenologyStandardDeviationFunctions
NULL

#' @rdname sd_functions
#' @export
SD.C =  function(mu_O=NA, sigma_O, mu_D=NA, sigma_D=NA, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="BB") {
		# if(is.na(sigma_D) || is.na(mu_O) || is.na(mu_C)) {
	  if(is.na(sigma_D) || is.na(mu_O)) {
			stop("The standard deviation for the duration, mean onset, and mean duration must be provided.")
		}
		SD.C.BB(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	}
	else if(type=="GP") {
		SD.C.GP(sigma=sigma_O)
	}
	else {
		stop(paste("Unsupported type: ", type))
	}
}

#' @rdname sd_functions
#' @export
SD.CkN = function(N, mu_O, sigma_O, mu_D, sigma_D=NA, minResponse=0, maxResponse=365, type=c("GP","BB"), intFailLow=NA, intFailHigh=NA) {
	type = match.arg(type)
	if(type=="BB") {
		if(is.na(sigma_D)) {
			stop("The standard deviation for the duration must be provided.")
		}
		SD.CkN.BB(N=N, mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	}
	else if(type=="GP") {
		mu_C = mu_O + mu_D
		SD.CkN.GP(N=N, mu_C=mu_C, sigma=sigma_O, minResponse=minResponse, maxResponse=maxResponse, intFailLow=intFailLow, intFailHigh=intFailHigh)
	}
	else {
		stop(paste("Unsupported type: ", type))
	}
}

#' @rdname sd_functions
#' @export
SD.D = function(mu_O=NA, sigma_O=NA, mu_D=NA, sigma_D=NA, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="BB") {
		if(is.na(sigma_D) || is.na(mu_O) || is.na(mu_D) || is.na(sigma_O)) {
			stop("The standard deviation for the duration, SD for onset, mean onset and mean duration must be provided.")
		}
		SD.D.BB(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	}
	else if(type=="GP") {
		SD.D.GP()
	}
	else {
		stop(paste("Unsupported type: ", type))
	}
}

#' @rdname sd_functions
#' @export
SD.O = function(sigma_O, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="BB") {
		SD.O.BB(sigma_O=sigma_O)
	}
	else if(type=="GP") {
		SD.O.GP(sigma=sigma_O)
	}
	else {
		stop(paste("Unsupported type: ", type))
	}
}

#' @rdname sd_functions
#' @export
SD.Ok1 = function(N, mu_O, sigma_O, minResponse=0, maxResponse=365, intFailLow=NA, intFailHigh=NA, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="BB") {
		SD.Ok1.BB(N=N, mu_O=mu_O, sigma_O=sigma_O, minResponse=minResponse, maxResponse=maxResponse)
	}
	else if(type=="GP") {
		SD.Ok1.GP(N=N, mu_O=mu_O, sigma=sigma_O, minResponse=minResponse, maxResponse=maxResponse, intFailLow=intFailLow, intFailHigh=intFailHigh)
	}
	else {
		stop(paste("Unsupported type: ", type))
	}
}

#' @rdname sd_functions
#' @export
SD.T = function(mu_O, sigma_O, mu_D, sigma_D=NA, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="BB") {
		if(is.na(sigma_D)) {
			stop("The standard deviation for the duration must be provided.")
		}
		SD.T.BB(mu_O=mu_O, sigma_O = sigma_O, mu_D=mu_D, sigma_D = sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	}
	else if(type=="GP") {
		mu_C = mu_O + mu_D
		SD.T.GP(mu_O=mu_O, mu_C=mu_C, sigma=sigma_O, minResponse=minResponse, maxResponse=maxResponse)
	}
	else {
		stop(paste("Unsupported type: ", type))
	}
}

#' Maximum a posteriori (MAP) point estimate
#' 
#' Finds the MAP under the beta onset, beta duration model (BB) using numerical optimization rather than MCMC techniques. The method does not include covariates.
#' 
#' This function may be useful to obtain estimates when the histogram of observed collection dates is highly skewed. The Gaussian process (GP) model implemented in Stan is limited to the inference of symmetric  distributions.
#'
#' @param responseData A vector of the observed collection times (e.g., day of year)
#' @param minResponse Minimum value of the response (e.g., day of year); must be set to 0 under current implementation (default = 0)
#' @param maxResponse Maximum value of the response (e.g., day of year); typically 365 for Gregorian calendar (default = 365)
#' @param minS Minimum allowable beta shape parameter value; must be positive and less than maxS (default = 1)
#' @param maxS Maximum allowable beta shape parameter value; must be positive and greater than minS (default = 3000)
#' @param init_params A vector of four elements that represent the parameter values of the starting state of the numerical optimization procedure. The four elements must be in this order: mean onset, standard deviation of onset, mean duration, standard deviation of duration. (default = c(180,20,60,7); default values appear to work well for most phenophases within a year time period)
#' @param hyperparameters A vector of 8 elements that represent the mean and standard deviation (sd) hyperparameter values for the prior distributions of the mean_onset, sd_onset, mean_duration, sd_duration parameters.  The eight elements must be in this order: mean mean_onset, sd mean_onset, mean sd_onset, sd sd_onset, mean mean_duration, sd mean_duration, mean sd_duration, sd sd_duration. (default: c(100, 7, 60, 6, 24, 12, 24, 12); these default values are somewhat arbitrarily set and should be set to values that are appropriate for your model)
#' @param type Currently, there is only one option: "BB". (default: "BB")
#'
#' @return List including the results of the R optim function with results min-max scaled, estimates of the four parameters at the original scale (par_orig_scale), a Boolean error status (error), and an error message when there is an error (error_m)
#' @export
#'
#' @examples
#' \donttest{
#' #Set the mean onset time
#' mean_onset = 100
#' #Set the onset time standard deviation
#' sigma_onset = 10
#' #Set the mean duration of the phenophase
#' mean_duration = 50
#' #Set the duration standard deviation
#' sigma_duration = 7
#' #Set the sample size
#' n=500
#' #Simulate observed collection times under the beta onset, beta duration (BB) model
#' ts = rT(n=n, mu_O=mean_onset, sigma_O=sigma_onset, mu_D=mean_duration
#'         , sigma_D=sigma_duration, type="BB")
#'
#' #estimate the MAP based on the simulated data and default initialization and hyperparameters.
#' map = getMAP(responseData=ts)
#' }
getMAP = function(responseData, minResponse=0, maxResponse=365,minS=1, maxS=3000,  init_params = c(180,20,60,7), hyperparameters = c(100, 7, 60, 6, 24, 12, 24, 12), type="BB") {
	if(type == "BB") {
		getMAP.T.BB(fileOrData=responseData, minResponse=minResponse, maxResponse=maxResponse, minS=minS, maxS=maxS,  init_params = init_params, hyperparameters = hyperparameters, dataProvided=TRUE)
	}
	else {
		stop(paste("Only model type beta onset, beta duration (BB) is implemented."))
	}
}

#' Maximum likelihood estimate of phenological parameters
#'
#' @description Numerically optimizes the likelihood function to find the combination of parameter values that result in the highest probability of the data, the maximum likelihood estimate (MLE). 
#' 
#' This method returns only a point estimate with no estimates of parameter uncertainty, and there is no use of prior knowledge. 
#' 
#' MLE can provide a useful check to determine if model identifiability is an issue. When a model is not identifiable, multiple combinations of parameter values result in the same or similar probability of the data. 
#' 
#' Currently only implemented for the beta onset, beta duration (BB) model. 
#'
#' @param responseData A vector of the observed collection times (e.g., day of year)
#' @param minResponse The minimum possible observed time. Current implementation requires this to be set to 0. (default: 0)
#' @param maxResponse The maximum possible observed time. For the Gregorian calendar, this is 365 (or 366 during leap years). (default: 365)
#' @param minS The smallest permissible value for a shape parameter for the beta distribution. This should be positive and less than maxS. (default: 1)
#' @param maxS The largest permissible value for a shape parameter fo the beta distribution. This should be positive and greater than minS. (default: 3000)
#' @param init_params A vector of four elements that represent the parameter values of the starting state of the numerical optimization procedure. The four elements must be in this order: mean onset, standard deviation of onset, mean duration, standard deviation of duration. Default values appear to work well for most phenophases within a year time period. (default = c(180,20,60,7)
#' @param type Currently, there is only one option: "BB". (default: "BB")
#'
#' @return A list that includes the estimate of the beta shape parameters from the R optim function (par), a vector of parameter values in the original scale (par_orig_scale), and a Boolean error flag
#' 
#' @export
#'
#' @examples
#' \donttest{
#' #Set the mean onset time
#' mean_onset = 100
#' #Set the onset time standard deviation
#' sigma_onset = 10
#' #Set the mean duration of the phenophase
#' mean_duration = 50
#' #Set the duration standard deviation
#' sigma_duration = 7
#' #Set the sample size
#' n=500
#' #Simulate observed collection times under the beta onset, beta duration (BB) model
#' ts = rT(n=n, mu_O=mean_onset, sigma_O=sigma_onset, mu_D=mean_duration, sigma_D=sigma_duration
#'         , type="BB")
#'
#' #estimate the MLE based on the simulated data and default initialization values
#' mle = getMLE(responseData=ts)
#' }
getMLE = function(responseData, minResponse=0, maxResponse=365, minS=1, maxS=3000, init_params = c(180, 20, 60, 7), type="BB") {
	if(type == "BB") {
		getMLE.T.BB(fileOrData=responseData, minResponse=minResponse, maxResponse=maxResponse, minS=minS, maxS=maxS, init_params = init_params, dataProvided=TRUE)
	}
	else {
		stop(paste("Only model type beta onset, beta duration (BB) is implemented."))
	}
}

#' Estimate the time of peak phenophase for a population
#'
#' Uses numerical optimization (BB model) or symmetry (GP model) to estimate the time when the population will have the highest percentage of individuals in the phenopase. 
#'
#' @param mu_O Mean onset time
#' @param sigma_O Standard deviation of the onset time
#' @param mu_D Mean duration time
#' @param sigma_D Standard deviation of the duration time
#' @param minResponse The minimum possible observed time. Current implementation requires this to be set to 0. (default: 0)
#' @param maxResponse The maximum possible observed time. For the Gregorian calendar, this is 365 (or 366 during leap years). (default: 365)
#' @param type The model type, either BB (beta onset, beta duration) or GP (Gaussian process with a shared standard deviation for onset and cessation and a constant duration) (default = "GP")
#'
#' @return Estimate of the peak phenophase time. 
#' @export
#'
#' @examples
#' \donttest{
#' #Set the mean onset time
#' mean_onset = 180
#' #Set the onset time standard deviation
#' sigma_onset = 5
#' #Set the mean duration of the phenophase
#' mean_duration = 60
#' #Set the duration standard deviation
#' sigma_duration = 40
#'
#' #estimate the peak phenophase under the BB model
#' peak = getPeakPhenophase(mu_O=mean_onset, sigma_O = sigma_onset, mu_D = mean_duration
#'                          , sigma_D = sigma_duration, type="BB")
#' #Create theoretical distribution of the observed times
#' curve(dT(x, mu_O=mean_onset, sigma_O = sigma_onset, mu_D = mean_duration
#'          , sigma_D = sigma_duration, type="BB"),from=0, to=365, n = 1000)
#' #Plot the estimate of the peak phenophase
#' points(peak, dT(peak, mu_O=mean_onset, sigma_O = sigma_onset, mu_D = mean_duration
#'       , sigma_D = sigma_duration, type="BB"), pch=16, col="yellow")
#' }
getPeakPhenophase = function(mu_O, sigma_O=NA, mu_D, sigma_D=NA, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="BB") {
		if(is.na(sigma_D) || is.na(sigma_O)) {
			stop("The standard deviations for the duration and onset must be provided.")
		}
		getPeak.T.BB(mu_O=mu_O, sigma_O = sigma_O, mu_D=mu_D, sigma_D = sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	}
	else if(type=="GP") {
		mu_C = mu_O + mu_D
		getPeak.T.GP(mu_O=mu_O, mu_C=mu_C)
	}
	else {
		stop(paste("Unsupported type: ", type))
	}
}

#' Fit the Weibull distribution to the distributions of phenological extremes
#'
#' @description Fits a Weibull distribution to the distribution of the first onset and the distribution of the last cessation. The Fisher–Tippett–Gnedenko theorem (AKA extreme value theorem) suggests that for distributions with bounded support, the distribution of the minimum (or the maximum) of a finite random sample will converge to the Weibull distribution under the appropriate scaling and translation. 
#'
#' This function is a wrapper for the ForestFit package function fitWeibull.
#'
#' @details
#' The resulting parameter estimates can be used with R's built-in \{r,d,q,p\} weibull functions. The x values need to be shifted by the location parameter (the third element of the returned onset or cessation vector). For the first onset, the values are shifted as follows: x = x - onsetLocation. For the last cessation, the values are shifted as follows: x = -(x - cessationLocation). As implemented in R, the Weibull shape parameter is the first element of the returned onset or cessation vector, and the Weibull scale parameter is the second element.
#'
#' Setting verbose to TRUE provides a print out of the specific call to dweibull using the estimated parameters. 
#'
#' @param N The population size for estimation of extreme events
#' @param mu_O Mean onset time
#' @param sigma_O Standard deviation for the onset time distribution
#' @param mu_D Mean phenophase duration length
#' @param sigma_D Standard deviation for the phenophase duration distribution
#' @param minResponse Minimum value of the response (e.g., day of year); must be set to 0 under current implementation (default = 0)
#' @param maxResponse Maximum value of the response (e.g., day of year); typically 365 for Gregorian calendar (default = 365)
#' @param type The model type, either BB (beta onset, beta duration) or GP (Gaussian process with a shared standard deviation for onset and cessation and a constant duration) (default = "GP")
#' @param precision Sample size to use to simulate data to which the Weibull is fit. 
#' @param includePlot Boolean whether to include a plot of the estimated Weibull distribution over simulated first onset and last cessation histograms. (default: FALSE)
#' @param verbose When set TRUE, provides detailed information about how to use the output from this function to obtain the Weibull estimates of the extreme value distributions. (default: TRUE)
#' @return A list with two vectors. The first vector provides the shape, scale, and location parameters of the Weibull for the first onset, whereas the second vector provides the corresponding estimates for the last cessation.
#' @export
#' @importFrom ForestFit fitWeibull
#' @importFrom graphics hist curve
#' @importFrom grDevices rgb
#' @importFrom stats dweibull
#'
#' @examples
#' \donttest{
#' #Set the mean onset time
#' mean_onset = 180
#' #Set the onset time standard deviation
#' sigma_onset = 5
#' #Set the mean duration of the phenophase
#' mean_duration = 60
#' #Set the duration standard deviation
#' sigma_duration = 40
#' #Set the sample size
#' precision = 10000
#' #Set the population size
#' N = 500
#' #Fit the Weibull distribution to the phenological extremes
#' fitWB = fitWeibullExtremes(N=N, mu_O=mean_onset, sigma_O=sigma_onset
#'                            , mu_D=mean_duration, sigma_D=sigma_duration, type="BB"
#'                            , precision=precision, includePlot=TRUE)
#' #Create histograms of random deviates based on the estimated parameter values
#' dev.new()
#' par(mfrow = c(2,1))
#' #Simulate draws from the first onset Weibull estimate, and overlay the theoretical curve
#' #      onto the histogram
#' firstOnset = rweibull(n=precision, shape = fitWB$Ok1Params[1]
#'                       , scale = fitWB$Ok1Params[2]) + fitWB$Ok1Params[3]
#' hist(firstOnset, probability=TRUE, xlab="First onset day of year", col="yellow")
#' curve(dweibull(x - fitWB$Ok1Params[3], shape = fitWB$Ok1Params[1]
#'                , scale = fitWB$Ok1Params[2]), add = TRUE)
#' #Simulate draws from the last cessation Weibull estimate, and overlay the theoretical curve
#' #       onto the histogram.
#' #Note that the mirror image is needed for the last cessation distribution.
#' lastCessation = -rweibull(n=precision, shape = fitWB$CkNParams[1]
#'                          , scale = fitWB$CkNParams[2]) + fitWB$CkNParams[3]
#' hist(lastCessation, probability=TRUE, xlab="Last cessation day of year", col="cyan")
#' curve(dweibull(-(x - fitWB$CkNParams[3]), shape = fitWB$CkNParams[1]
#'                , scale = fitWB$CkNParams[2]), add = TRUE)
#' }
fitWeibullExtremes = function(N, mu_O, sigma_O, mu_D, sigma_D=NA, minResponse=0, maxResponse=365, type=c("GP","BB"), precision=10000, includePlot=FALSE, verbose=TRUE) {
	type = match.arg(type)

	if(type=="GP") {
		mu_C = mu_O + mu_D
		Ok1 = rOk1.GP(n=precision, N=N, mu_O=mu_O, sigma=sigma_O)
		CkN = rCkN.GP(n=precision, N=N, mu_C=mu_C, sigma=sigma_O)
	}
	else if(type=="BB") {
		if(is.na(sigma_D)) {
			stop("Please provide the standard deviation of the phenophase duration distribution.")
		}
	  ## DANIEL: Changed sigma argument to sigma_O on rOk1.BB call.
		Ok1 = rOk1.BB(n=precision, N=N, mu_O=mu_O, sigma_O=sigma_O, minResponse=minResponse, maxResponse=maxResponse)
		CkN = rCkN.BB(n=precision, N=N, mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	}
	else {
		stop(paste("Unsupported model type: ", type))
	}

	starts <- c(2,2,0)
	mOk1 = mean(Ok1)
	mCkN = mean(CkN)
	fitOk1 <- fitWeibull((Ok1-mOk1), location = TRUE, method = "mps", starts = starts)
	fitCkN <- fitWeibull(-(CkN-mCkN), location = TRUE, method = "mps", starts = starts)

	if(includePlot) {
		xlim = range(Ok1, CkN)

		nBin = round(3*100/2)
		inc = (xlim[2]-xlim[1])/((2/3)*nBin)
		x = seq(xlim[1]-inc,xlim[2]+inc, inc)

		hist(Ok1, breaks=x, col=rgb(1,1,0,0.5), xlim=c(xlim[1]-inc, xlim[2]+inc), probability=TRUE, xlab="Collection date", main=NULL)
		hist(CkN, breaks=x, col=rgb(0,1,1,0.5),  add=TRUE, probability=TRUE)

		fitOk1$estimate[3] = fitOk1$estimate[3] + mOk1
		fitCkN$estimate[3] = mCkN - fitCkN$estimate[3]
		curve(dweibull(x - fitOk1$estimate[3], shape = fitOk1$estimate[1], scale = fitOk1$estimate[2]), add = TRUE, col = "black", lwd = 2, lty = 2, n = precision)
		curve(dweibull(-(x - fitCkN$estimate[3]), shape = fitCkN$estimate[1], scale = fitCkN$estimate[2]), add = TRUE, col = "black", lwd = 2, lty = 2, n = precision)
	}

	cat("To obtain the fit density values for your data, use the following command ('res' is the output of this function, and x is a vector of your x values):\n")
	cat("First onset:\n\tdweibull(x - res$Ok1Params[3], shape = res$Ok1Params[1], scale = res$Ok1Params[2])\n")
	cat(paste0("\tspecifically: dweibull(x - ", fitOk1$estimate[3], ", shape = ", fitOk1$estimate[1], ", scale = ", fitOk1$estimate[2], ")\n"))
	cat("Last cessation:\n\tdweibull(-(x - res$CkNParams[3], shape = res$CkNParams[1], scale = res$CkNParams[2])\n")
	cat(paste0("\tspecifically: dweibull(x - ", fitCkN$estimate[3], ", shape = ", fitCkN$estimate[1], ", scale = ", fitCkN$estimate[2], ")\n"))

	return(list( CkNParams = fitCkN$estimate, Ok1Params = fitOk1$estimate ))
}

#' Run a Bayesian analysis of phenological events using biocollection data
#'
#' @description
#' This function is the core of the phenoCollectR package. It wraps many other functions to make the Bayesian analysis of phenological events using biocollection data as streamlined as possible, yet flexible if need be. In particular, it aims to connect variation in life history timing to environmental variation. 
#'
#' Models with or without covariates can be used. Covariates can be climate variables, physical environment variables, physiological states, genetic states, and so forth - anything that may contribute to phenological variation. However, as currently implemented, only continuous covariates are supported.
#' 
#' Mean phenophase onset time and mean phenophase duration are modeled as linear functions of the covariates, if supplied. Otherwise, just the mean phenophase onset time and the mean phenophase duration are modeled without covariates. Separate covariates can be used to model the mean onset time and the mean phenophase duration. However, some downstream methods for summaries and plotting expect the covariates to be the same for both models.
#' 
#' The helper function 'preparePhenologyData' can be used to prepare the response data and the covariate data for input to this function.
#' 
#' The theory developed by Hearn et al. (XXXX) unifies the models for onset and duration to derive the theoretical distribution of observed collection times. Bayesian inference is carried out using the probabilistic programming language Stan, which draws samples from the posterior distribution of the model parameters. These samples can be used for downstream analyses. For example, to estimate the marginal posterior probability that a slope coefficient associated with a covariate is negative, one simply counts the number of posterior draws where the coefficient is negative and divides by the total number of draws from the posterior distribution.
#'
#' As part of a Bayesian analysis, and to reduce model identifiability issues with presence-only data, a prior distribution is needed for the analysis. Flat priors will not work for this analysis. Because the setting of priors can be challenging, this function can optionally implement methods that will automate the specification of priors for you. Advanced users can manually tailor prior hyperparameters and achieve weakly informative priors, if desired. Under the current implementation, the functional form of the priors is pre-set. 
#' 
#' Users can also optionally choose to analyze phenological extremes. To do so requires an estimate of the population size from with biocollection records were obtained. The concept of "population" in this case can be quite loose and can be the metapopulation consisting of all populations from which samples were selected. 
#' 
#' The parameters to be estimated for the 'full' model include the marginal mean onset time (\eqn{\mu_O}), the marginal mean duration (\eqn{\mu_D}), one parameter (\eqn{\beta_O} coefficient) for each covariate used to model the mean onset time, one parameter (\eqn{\beta_D} coefficient) for each covariate used to model the mean duration time, and one \eqn{\sigma_O} parameter which models both the standard deviation of the distribution of onset times and the standard deviation of the distribution of durations. If there are \eqn{K_O} onset covariates, and \eqn{K_D} duration covariates, then there is a total of 3 + \eqn{K_O} + \eqn{K_D} free parameters (including the two \eqn{\alpha} intercept terms. 
#'
#' The parameters to be estimated for the 'intercept-only' model include the marginal mean onset time (\eqn{\mu_O}), the marginal mean duration (\eqn{\mu_D}), and the sigma parameter (\eqn{\sigma_O}). 
#'
#' Other parameters associated with phenophase cessation events and phenological extremes are deterministic functions of the free parameters based on the implemented Gaussian process model.
#'
#' A summary of the Stan run is available through the helper function 'summarizePhenologyResults', or advanced users can access the Stan sample directly from the output of this function for customized downstream analyses using other packages designed for the analysis of samples from posterior distributions, such as 'posterior' and 'bayesplot'. 
#' 
#' \bold{\emph{BACKGROUND FRAMEWORK AND THEORY}}
#'
#' \bold{DATA TYPES}
#' 
#' There are two main types of datasets that can be analyzed. When only observation times of individuals in the stage of interest are available, these are presence-only data (PO). When observation times are available before, during, and after the stage of interest, these are before-during-after data (BDA). Theory for these two types of datasets differs builds on itself.
#' 
#' \bold{NOTATION}
#' 
#' Both visible and latent events are modeled with random variables (RVs). 
#' 
#' RVs are symbolized with the first letter of the state of interest: Onset (O), Duration (D), Cessation (C), observed collection Times (T). 
#' 
#' RVs representing extremes are subscripted. For example, first onset is \eqn{O_1}, and last onset is \eqn{O_N} where \eqn{N} is the population size. 
#' 
#' Parameters are subscripted by the associated random variable. For example, mean duration is \eqn{\mu_D} and the standard deviation of the last cessation is \eqn{\sigma_{C_N}}. In the linear models, intercept parameters are represented by \eqn{\alpha} and slope coefficients are \eqn{\beta}. So, for the model of mean onset, the intercept is \eqn{\alpha_O} and a slope coefficient for the \eqn{i^{th}} covariate is \eqn{\beta_{O,i}}.
#' 
#' The model implemented in Stan is based on the properties of a degenerate Gaussian process that samples lines with intercepts whose values are normally distributed. When a line passes a first threshold at \eqn{\mu_O}, the process begins, and when that line passes a second threshold at \eqn{\mu_D}, the process ends. This process results in normally-distributed onsets: \eqn{O \sim N(\mu_O, \sigma^2_O)}, a constant duration \eqn{D}, and a normally-distribed cessation: \eqn{C \sim N(\mu_O + D, \sigma^2_O)}.
#' 
#' The process of interest occurs in a time period (TP), such as a year for cyclical annual processes, or a day for cyclical diel processes, that has a minimum time of \eqn{m} and a maximum time of \eqn{M}. For a year, the day of year time is between \eqn{m=0} and \eqn{M=365}. 
#' 
#' The dataset consists of a population of individuals undergoing the process of interest contained within the TP. Onset and duration times can vary for each individual in the population. 
#' 
#' Onset times have intrinsic noise due to \eqn{\sigma_O}, whereas durations do not. Both onset times and duration lengths can vary in relation to internal states such as energy store level or organismal size, or to external conditions such as climate variation (temperature summaries for spring, etc.) or environmental variation (soil type, moisture levels, light levels, etc.). Any number of covariates can be used with standard cautions about overly complex models (overfitting, model identifiability, parameter uncertainty, collinearity, etc.).
#' 
#' The sample design is as follows. Individuals are randomly sampled from a population at random points in time during the time period.  
#' 
#' The only RVs with directly observed states are the observed times (\eqn{T}) and the observed states (\eqn{S}). All the other RVs are latent (hidden) and estimated from the observed times and associated stages.
#' 
#' For the theory, we simplify things by examining a time period (TP) starting at time 0 and of length 1 time unit. Data sampled from other length time periods can be translated to start at time 0 by subtracting \eqn{m} from all times and scaled to a length of 1 by dividing all sampled times by \eqn{M-m}. In what follows, we assume such transformations have already been made. These transformations are handled automatically by this runStanPhenology function. If users wish to transform their data in other ways prior to analysis, we leave this up to the user. 
#' 
#' \bold{THEORY: BEFORE-DURING-AFTER {BDA}}
#' 
#' For \bold{before-during-after (BDA)} data, all observed times are kept and the stage is recorded.
#' 
#' For the core theory, we need to determine the probabilities of each stage given the sampled time and model parameters \eqn{\theta \equiv (D,\mu_O,\sigma_O)}. 
#' 
#' The probability that a time is before the onset is \deqn{P_{before}(O>t|\theta)=1-\Phi(\frac{t-\mu_O}{\sigma_O})} where \eqn{\Phi(z) \equiv P(Z<z)} is the cumulative distribution function of a standard normal random variable \eqn{Z}. 
#' 
#' With a duration \eqn{D}, the probability that a sampled time \eqn{t} occurs after the cessation \eqn{C=O+D} is \deqn{P_{after}(t \ge O+D | \theta) = P_{after}(O \le t-D | \theta) = \Phi(\frac{t-D-\mu_O}{\sigma_O})}
#' 
#' The "during event" is the complement of the "after event" unioned with the "before event", since before, after, and during are exhaustive and mutually exclusive events. So, the probability that the sampled time \eqn{t} occurs between onset \eqn{O} and cessation \eqn{O+D} follows from the complement rule of probability: \deqn{P_{during}(O \le t \lt O+D|\theta) = 1 - (P_{after}+P_{before})=\Phi(\frac{t-\mu_O}{\sigma_O})-\Phi(\frac{t-D-\mu_O}{\sigma_O})}
#' 
#' For stages before, during, and after the process of interest, symbolized as 1, 2, and 3, respectively, the likelihood function for \eqn{n} randomly sampled times \eqn{\vec{t}} with associated stages \eqn{\vec{s}} is then \deqn{P(\vec{s}|\vec{t},\theta)=\prod_{i=1}^{i=n}P_{before}^{I(s_i=1)}P_{during}^{I(s_i=2)}P_{before}^{I(s_i=3)}} where the indicator function \eqn{I(s_i=j)} returns 1 whenever the stage for the \eqn{i^{th}} individual equals \eqn{j}.
#' 
#' \bold{THEORY: PRESENCE-ONLY (PO)}
#' 
#' For \bold{presence-only (PO)} data, only the observed times \emph{during} the process of interest are recorded and all other times are discarded. 
#' 
#' The theory for presence only is a bit more complicated than the theory for before-during-after data. The ultimate goal for PO analysis is to infer the parameters that describe the distribution of observed sampled times. The other distributions of interest, including extremes, onsets, peak, and cessation, are functions of these parameters with a caveat for extremes. 
#' 
#' The distribution whose parameters we want to infer is the probability density of the observed times given the individual was sampled during the process of interest (state \eqn{s=1})and given the parameters of the onset distribution and duration: \deqn{p(t|s=1,D,\mu_O,\sigma_O)} For simplicity, let \eqn{\theta \equiv (D,\mu_O,\sigma_O)}. Using Bayes rule, \deqn{p(t|s=1,\theta) = \frac{p(s=1|t,\theta)p(t|\theta)}{p(s=1|\theta)}}
#' 
#' The probability that the stage is during the process given the sampled time and parameters was derived above as \eqn{P_{during}}. 
#'
#' The sampling strategy defines \eqn{p(t|\theta)}. In the current implementation, the sampling of times is at random and independent of \eqn{\theta}, i.e., \eqn{p(t|\theta)=\frac{1}{M-m}}. Since all observed times are scaled to a TP of length 1 already, \eqn{p(t|\theta)=1}. 
#' 
#' Finally, \eqn{p(s=1|\theta)} is the length of the process of interest when the TP is length 1 and the process occurs entirely within the TP. By the total law of probability this is \deqn{\int_0^1 p(s=1,t|\theta)dt=\int_0^1 p(s=1|t,\theta)p(t|\theta)dt=\int_0^1 p(s=1|\theta)dt=\int_0^1p_{present}dt=\int_0^1 (\Phi(\frac{t-\mu_O}{\sigma_O})-\Phi(\frac{t-D-\mu_O}{\sigma_O}))dt} Hearn et al. (XXXX) in their Appendix 1 provide a closed form for this integral.
#' 
#' \bold{THEORY: MODELING COVARIATES}
#' 
#' In the current implementation, the mean onset time and the duration time are modeled as linear functions of covariates. Different covariates can be used for the onset model and for the duration model, but here, for simplicity, they are assumed to be the same for both models. If we let \eqn{\vec{X}} represent the vector of covariates and \eqn{\vec{\beta}} represent the vector of coefficients, then the following are the models for onset and duration: \deqn{O \sim N(\alpha_O + \vec{\beta_O} \cdot \vec{X}, \sigma_O^2)} and \deqn{D \equiv \mu_D = \alpha_D+\vec{\beta_D} \cdot \vec{X}}
#' 
#' These models of mean onset and duration make the model hierarchical, so that with the inclusion of covariates, everywhere in the likelihood function where \eqn{\mu_O} or \eqn{D} appear, they are replaced by the above linear functions. If there are \eqn{K_O} covariates for the onset model and \eqn{K_D} covariates for the duration model, along with the two \eqn{\alpha} intercepts and the \eqn{\sigma_O} standard deviation parameter, the fully specified model has \eqn{2 + 1 + K_O + K_D} estimable parameters.
#' 
#' \bold{THEORY: BAYESIAN ANALYSIS AND PRIORS}
#' 
#' In a Bayesian analysis the probability density of the parameters given the data, which is called the posterior distribution, is proportional to the product of the likelihood, which is the probability of the data given the parameters, and the prior distribution, which is the marginal probability of the parameters: \eqn{Posterior \propto Likehood * Prior}. 
#' 
#' The above sections provided derivations of the likelihood functions for both PO datasets and BDA datasets, so the only remaining thing needed to carry out a Bayesian analysis is to specify the prior distributions. The priors are defined for the \eqn{D,\beta_O,\alpha_O, \beta_D, \alpha_D, \sigma_O} parameters. In the current implementation, the priors can be flat, or they can be normally distributed. 
#' 
#' Set the \emph{priorLevel} function parameter (described below) to 0 for a flat prior. This is only recommended for use with BDA datasets, as identifiability for PO sets necessitates an informative prior. 
#' 
#' Set the \emph{priorLevel} to any value greater than 1 for normally distributed priors. A normally distributed prior requires two hyperparameters for each parameter: the mean of the parameter and the standard deviation of the parameter. Users can manually input these hyperparameter values or can use automated prior specification. 
#' 
#' Normal priors are used primarily for convenience and interpretability. Under standard regularity conditions, the Bernstein–von Mises theorem implies that as sample size increases, the posterior distribution becomes approximately multivariate normal, centered near the maximum likelihood estimate, with covariance given by the inverse Fisher information scaled by \eqn{1/n}. In this large-sample regime, posterior uncertainty is therefore well summarized by normal approximations and depends only weakly on the particular prior, provided the prior assigns positive density in a neighborhood of the true parameter value. This asymptotic behavior makes normal priors a natural default choice in many settings. Moreover, in sequential Bayesian analyses where researchers use a previous large-sample posterior as the prior for a subsequent analysis, the inherited prior will itself be approximately normal.
#' 
#'
#' @param type The type of model to be implemented. Four options are available. Option 'full' will model mean onset and mean duration as a linear function of covariates for presence-only data, option 'interval-full' will model mean onset and mean duration as a linear function of covariates for before-during-after data, option 'intercept-only' will model the distribution of onset times and the distribution of durations with no consideration of covariates for presence-only data, and option 'interval' will model the distribution of onset times and durations with no consideration of covariates for before-during-after data. The Examples section below provides example usage of each scenario (default: intercept-only)
#' @param responseData A vector of the response data. The response data are simply the observed times of collection, often the day of year for biocollection data. This vector is prepared by the 'preparePhenologyData' function, or advanced users can prepare it using whatever tools they choose. The data should be in the original scale unless advanced users have a reason to use preprocessing transformations. (default: NULL)
#' @param stage Only for use when type is 'interval' or 'interval-full'. A vector of integers the same length as the number of observations. For each observation, the stage will be a 1 if the observation is before the phenophase, a 2 if it is during the phenophase, and a 3 if it is after the phenophase. Any other values will produce errors. (default: NULL)
#' @param hyperparams_noCovariates For use with 'intercept-only' or 'interval' models. A vector of six elements that represent the mean and standard deviation (sd) hyperparameter values for the prior distributions of the mean onset, the mean duration, and the sigma parameters, in that order. Items in the vector can be named or unnamed. (default: NULL)
#' @param onsetCovariateData For use with '*full' models. A data frame with each column representing a covariate and each row representing a specimen observation. The number of rows must match the number of items in the responseData vector, and specimens are in the same order as they are in the response data vector. The function 'preparePhenologyData' can prepare this data frame, or advanced users can prepare the data frame using their own tools. Data should be in the original scale, or advanced users can transform the data as deemed appropriate. (default: NULL)
#' @param durationCovariateData For use with '*full' models. A data frame with each column representing a covariate and each row representing a specimen observation. The number of rows must match the number of items in the responseData vector, and specimens are in the same order as they are in the response data vector. The function 'preparePhenologyData' can prepare this data frame, or advanced users can prepare the data frame using their own tools. Data should be in the original scale, or advanced users can transform the data as deemed appropriate. (default: NULL)
#' @param onsetHyperBeta For use with '*full' models. A data frame with two columns to provide hyperparameter values for the covariate slope coefficients for the onset model. The first column is the the mean hyperparameter, and the second is the standard deviation parameter. Each row corresponds with a separate covariate. The first row in the data frame corresponds to the covariate in the first column of the covariate data frame (onsetCovariateData), the second row with the second covariate, and so forth. If there are \eqn{K_O} columns in the covariate data frame, there should be \eqn{K_O} rows in this data frame. Can be left at default when priorLevel is set to 0, in which case Stan default priors are used (not recommended). (default: NULL) 
#' @param onsetHyperAnchor For use with '*full' models. A two-element vector with the mean of the prior and the standard deviation of the prior for the onset anchor (the marginal mean onset time value). Can be left at default when priorLevel is set to 0, in which case Stan default priors are used (not recommended). (default: NULL) 
#' @param durationHyperBeta For use with '*full' models. A data frame with two columns to provide hyperparameter values for the covariate slope coefficients for the duration model. The first column is the the mean hyperparameter, and the second is the standard deviation parameter. Each row corresponds with a separate covariate. The first row in the data frame corresponds to the covariate in the first column of the covariate data frame (durationCovariateData), the second row with the second covariate, and so forth. If there are \eqn{K_D} columns in the covariate data frame, there should be \eqn{K_D} rows in this data frame. Can be left at default when priorLevel is set to 0, in which case Stan default priors are used (not recommended). (default: NULL) 
#' @param durationHyperAnchor For use with '*full' models. A two-element vector with the mean of the prior and the standard deviation of the prior for the duration anchor (the marginal mean phenophase duration value). Can be left at default when priorLevel is set to 0, in which case Stan default priors are used (not recommended). (default: NULL) 
#' @param sigmaHyper For use with '*full' models. A two-element vector with the mean and standard deviation of the prior distribution for the sigma parameter. Sigma is the standard deviation of the onset distribution and is also the standard deviation of the cessation distribution. Can be left at default when priorLevel is set to 0, in which case Stan default priors are used (not recommended). (default: NULL) 
#' @param minResponse Minimum value of the response (e.g., day of year); must be set to 0 (i.e., at default) under current implementation (default = 0)
#' @param maxResponse Maximum value of the response (e.g., day of year); typically 365 for Gregorian calendar (default = 365)
#' @param keepScale Leave the response data at original scale during Stan run. Only implemented for type "interval". Default is usually fine. Use keepScale if you need the Stan sampler to initialize the sampling at the hyperparameter values and keep the original scale. Best to keep at default. (default: FALSE, which min-max scales the data)
#' @param maxDiv The maximum number of divergences to be tolerated during the Stan posterior sampling. This should be 0 unless a biased sample from the posterior is acceptable. (default: 0)
#' @param setStringent Boolean flag to indicate more stringent sampling during Stan runs. Specifically, adapt_delta = 0.99, max_treedepth = 15. Setting to TRUE reduces the chances of divergences, but run time is slower. Usually the difference in run time is not of practical concern, so the default is fine. (default: TRUE)
#' @param runMAP Boolean flag to indicate if Stan should be used to estimate the maximum a posteriori (MAP) values. (default: TRUE)
#' @param processExtremes Boolean flag to indicate if parameters for the phenological extreme distributions should be estimated during the analysis. If set to TRUE, asymptotic estimates of the expected first onset and expected last cessation are sampled. Because asymptotic estimates are provided, estimates are made very quickly. A limitation is that an estimate of the population size is needed for this analysis. (default: TRUE)
#' @param N The population size. Needed when processExtremes is set to TRUE. (default: 500)
#' @param partitionDataForPriors Boolean flag to indicate if methods to automate the specification of prior hyperparameters should be used. In particular, 30% of the data are partitioned into a set to estimate the hyperparameters, and 70% of the data are used for the Bayesian analysis using Stan. If set to TRUE, all other user-provided hyperparameter values (e.g., hyperparams_noCovariates, onsetHyperBeta, onsetHyperAnchor, etc.) are ignored. This is not recommended if accurate estimates of duration and onset parameters are needed, but estimates of covariate coefficients tend to be accurate. (default: FALSE)
#' @param maximizeSampleSize Boolean flag for use when partitionDataForPriors is set to TRUE. If a user wants 100% of the data to be used for the Bayesian analysis, set this to TRUE. This is statistically invalid because the same data are used to estimate the prior hyperparameters and carry out the full Bayesian analysis. When the sample size is very small, and no prior information is available, the error on the estimates may be unacceptably large unless the fuller dataset is used for inference. Recommended to keep at default. (default: FALSE)
#' @param threshApprox An error threshold set to use approximation schemes when numerical integration fails. Safe to leave at default. (default: NULL)
#' @param byPassChecks Boolean flag indicating whether to bypass checks for Stan. If Stan is not located, an attempt will be made to install Stan. Recommended to keep at default. (default: FALSE)
#' @param priorLevel Specifies types of priors by an integer. Setting to 0 will use flat priors. Above 0 will use the user-input prior hyperparameters with normally distributed priors. For presence-only data, leave at default. (default: 2)
#' @param ... Parameters to be input into the Stan sample function. Not currently implemented. Do not use.
#'
#' @return Details of what is returned will depend on which 'type' of model is used. See 'type' parameter.
#'
#' In the case of \bold{'*full' type}, a list with the following items is returned:
#'
#' data: the data that were passed to Stan. These include the scaled and translated response and covariate data.
#'
#' responseData: the original response data
#'
#' onsetCovariateData: the original covariate data for the onset model
#'
#' durationCovariateData: the original covariate data for the duration model
#'
#' onsetHyperBeta: the mean and standard deviation hyperparameter values for the onset model covariate coefficients
#'
#' onsetHyperAnchor: the mean and standard deviation hyperparameter values for the onset anchor (marginal mean onset)
#'
#' durationHyperBeta: the mean and standard deviation hyperparameter values for the duration model covariate coefficients
#'
#' durationHyperAnchor: the mean and standard deviation hyperparameter values for the duration anchor (marginal mean duration)
#'
#' sigmaHyper: the mean and standard deviation hyperparameter values for sigma
#'
#' result: the sample from the Stan sampling
#'
#' model: the Stan model
#'
#' maxDiv: the maximum allowed number of divergences for the run
#'
#' minResponse: the minimum possible response value
#'
#' maxResponse: the maximum possible response value
#'
#' setStringent: whether setStringent was set to TRUE.
#'
#' error: whether an error was caught during the Stan run. 
#'
#' error_m: an error message, possibly indicating no error occurred.
#'
#' In the case of \bold{'intercept-only' or 'interval'} type, a list with the following items is returned:
#'
#' data: the data that were passed to Stan. These include the scaled and translated response and covariate data.
#'
#' responseData: the original response data
#'
#' hyperparameters: the hyperparameters used by the Stan model
#'
#' sample: the sample from Stan
#'
#' model: the Stan model
#'
#' runMap: Boolean indicates of MAP estimate was made
#'
#' MAP: Stan's MAP estimate, if selected. NA if runMap is FALSE
#' 
#' N: input population size
#'
#' minResponse: the minimum possible response value
#'
#' maxResponse: the maximum possible response value
#'
#' setStringent: whether setStringent was set to TRUE.
#'
#' maxDiv: the maximum allowed number of divergences for the run
#'
#' error: whether an error was caught during the Stan run. 
#'
#' error_m: an error message, possibly indicating no error occurred.
#' 
#' threshApprox: error threshold value for integration
#' @export
#'
#' @examples
#' \donttest{
#' ##########################################################################################################################################################
#' ##Run an example with EMPIRICAL COLLECTION TIMES of individuals in the phenophase (PRESENCE ONLY - PO) - type 'full'.
#'
#' ##DEFAULT PRIORS are used and are only reliable for slope coefficient estimations. 
#' #!!*Use a prior based on biological data for the duration if accurate duration estimation is a priority.*!!
#' #
#' #Get the FILE NAME with data for the blood root plant. Data files for 12 other species are also available
#' file  =  getDatasetPath("Sanguinaria_canadensis")
#' #
#' ## ACCESS MORE SPECIES
#' help(getDatasetPath)
#' #
#' #COVARIATE NAMES - remove up to all but 1
#' vars = c("Latitude", "Year", "Elevation", "AnnualMonthlyAverageTemp"
#'          , "SpringMonthlyAverageTemp", "FirstQuarterMonthlyAverageTemp")
#' #
#' #EXTRACT PHENOLOGY DATA
#' data  =  preparePhenologyData(dataFile=file, responseVariableName="DOY"
#'                               , onsetCovariateNames=vars, durationCovariateNames=vars
#'                               , taxonName="Sanguinaria_canadensis", removeOutliers=TRUE)
#'
#' \bold{#EMPIRICAL DATA STAN RUN WITH COVARIATES}
#' stanResult  =  runStanPhenology(type="full", responseData = data$responseData
#'                                 , onsetCovariateData = data$onsetCovariateData
#'                                 , durationCovariateData = data$durationCovariateData
#'                                 , partitionDataForPriors = TRUE)
#' #
#' #SUMMARIZE 
#' stanSummary  =  summarizePhenologyResults(stanRunResult = stanResult
#'                                           , taxonName = "Sanguinaria_canadensis"
#'                                           ,standardLinearModel = TRUE)
#' #WRITE EMPIRICAL RESULTS
#' write.csv(stanSummary,"SCanadensis.Empirical.wCovariates.csv")
#' ##
#' ##
#' ##########################################################################################################################################################
#' ##Run more complicated examples, this time with SIMULATED PHENOLOGY and SIMULATED COVARIATES 
#' #
#' #COVARIATE MODEL
#' covariate_namesOnset = c("v1", "v2")
#' covariate_namesDuration = c("v1", "v3")
#' cov_all = union(covariate_namesOnset, covariate_namesDuration)
#' covarianceMatrix = matrix(c( 1.0, 0.5, 0.3,
#' 				      0.5, 2.0, 0.4,
#' 				      0.3, 0.4, 1.5), nrow = 3, byrow = TRUE)
#' rownames(covarianceMatrix) = cov_all
#' colnames(covarianceMatrix) = cov_all
#' covariateMeans = c(50,60,80)
#' names(covariateMeans) = cov_all
#' #
#' #ONSET PARAMETERS
#' slopesOnset = c(1,-2)
#' names(slopesOnset) = covariate_namesOnset
#' mean_responseOnset = 150
#' noiseOnset = 3
#' #
#' #DURATION PARAMETERS
#' slopesDuration = c(3,1)
#' names(slopesDuration) = covariate_namesDuration
#' mean_responseDuration = 30
#' #
#' #SAMPLE SIZE
#' n=3000
#' #
#' #SIMULATE PHENOLOGY TIMES, OBSERVATION TIMES, AND COVARIATE DATA
#' simulated_data = simulatePopulationLatentIntervalStates(n=n,
#'                        betaOnset=slopesOnset, betaDuration=slopesDuration,
#'                covariateNamesOnset=covariate_namesOnset, covariateNamesDuration=covariate_namesDuration,
#'                covariateMeans = covariateMeans,
#'                covarianceMatrix = covarianceMatrix,
#'                meanOnset = mean_responseOnset, meanDuration = mean_responseDuration,
#'                sigmaOnset = noiseOnset)
#' #
#' #_________________________________________________________________________________________________________________________________________________________
#' #
#' #HYPERPARAMETERS
#' onsetHyperBeta = data.frame(slopesOnset,c(1,1))
#' onsetHyperAnchor = c(150, 5)
#' durationHyperBeta = data.frame(slopesDuration,c(1,1))
#' durationHyperAnchor = c(30, 3)
#' sigmaHyper = c(3,1)
#' #
#' #_________________________________________________________________________________________________________________________________________________________
#' #
#' #BDA DATA (Before, during, after)
#' responseDataBDA = simulated_data$observedTime
#' onsetCovariateDataBDA = simulated_data[,covariate_namesOnset, drop=FALSE]
#' durationCovariateDataBDA = simulated_data[,covariate_namesDuration, drop=FALSE]
#' stageBDA = simulated_data$stage
#' #
#' \bold{#BDA STAN RUN - type 'interval-full'}
#' stanResultBDA = runStanPhenology(
#' 			type="interval-full",
#' 			responseData=responseDataBDA,
#' 			stage=stageBDA,
#' 			onsetCovariateData=onsetCovariateDataBDA, durationCovariateData=durationCovariateDataBDA,
#' 			onsetHyperBeta=onsetHyperBeta,onsetHyperAnchor=onsetHyperAnchor,
#' 			durationHyperBeta=durationHyperBeta,durationHyperAnchor=durationHyperAnchor,
#' 			sigmaHyper=sigmaHyper,
#' 			minResponse=0,maxResponse=365,
#' 			runMAP=TRUE,processExtremes=TRUE,N=500,
#' 			partitionDataForPriors=FALSE,
#' 			maxDiv=0,setStringent=TRUE,
#' 			priorLevel=2
#' 			)
#' #
#' ##SUMMARIZE
#' stanSummaryBDA  =  summarizePhenologyResults(stanRunResult = stanResultBDA
#'                                           , taxonName = "Simulated"
#'                                           ,standardLinearModel = FALSE)
#' ##WRITE BDA RESULTS 
#' write.csv(stanSummaryBDA,"Simulated.BeforeDuringAfter.informativePriors.csv")
#' #
#' #_________________________________________________________________________________________________________________________________________________________
#' #
#' #PO DATA 
#' responseData = simulated_data$observedTime[simulated_data$stage==2]
#' onsetCovariateData = data.frame(simulated_data$v1[simulated_data$stage==2],simulated_data$v2[simulated_data$stage==2])
#' durationCovariateData = data.frame(simulated_data$v1[simulated_data$stage==2],simulated_data$v3[simulated_data$stage==2])
#' colnames( onsetCovariateData ) = covariate_namesOnset
#' colnames( durationCovariateData ) = covariate_namesDuration
#' #
#' \bold{#PO STAN RUN - type 'full'}
#' stanResultPO = runStanPhenology(
#' 			type="full",
#' 			responseData=responseData,
#' 			onsetCovariateData=onsetCovariateData, durationCovariateData=durationCovariateData,
#' 			onsetHyperBeta=onsetHyperBeta,onsetHyperAnchor=onsetHyperAnchor,
#' 			durationHyperBeta=durationHyperBeta,durationHyperAnchor=durationHyperAnchor,
#' 			sigmaHyper=sigmaHyper,
#' 			minResponse=0,maxResponse=365,
#' 			runMAP=TRUE,processExtremes=TRUE,N=500,
#' 			partitionDataForPriors=FALSE,
#' 			maxDiv=0,setStringent=TRUE,
#' 			priorLevel=2
#' 			)
#' ##summarize the Stan run
#' stanSummaryPO  =  summarizePhenologyResults(stanRunResult = stanResultPO
#'                                           , taxonName = "Simulated"
#'                                           ,standardLinearModel = FALSE)
#' ##write results to file
#' write.csv(stanSummaryPO,"Simulated.PresenceOnly.informativePriors.csv")
#' ##
#' ##
#' ##
#' ##########################################################################################################################################################
#' ##Conduct an analysis with NO COVARIATES (NC) using SIMULATED DATA
#' #
#' #SIMULATION PARAMETERS 
#' mean_responseOnset = 150
#' sdOnset = 3
#' mean_responseDuration = 30
#' #
#' #SAMPLE SIZE
#' n=5000
#' #
#' #SIMULATE
#' simulated_data = simulatePopulationLatentIntervalStates(n=n,
#'                meanOnset = mean_responseOnset, meanDuration = mean_responseDuration,
#'                sigmaOnset = sdOnset)
#' #
#' #_________________________________________________________________________________________________________________________________________________________
#' #
#' #HYPERPARAMETERS (biased): mean meanOnset, sd meanOnset, mean meanDuration, sd meanDuration, mean sigma, sd sigma)
#' hyperparameters = c(120, 30, 35, 5, 5, 2)
#' #
#' #_________________________________________________________________________________________________________________________________________________________
#' #
#' #BDA DATA
#' responseDataBDA.NC = simulated_data$observedTime
#' stageBDA.NC = simulated_data$stage
#' \bold{#STAN RUN BDA - type 'interval'}
#' stanResultBDA.NC = runStanPhenology(
#' 			type="interval",
#' 			responseData=responseDataBDA.NC,
#' 			stage=stageBDA.NC,
#'			hyperparams_noCovariates=hyperparameters,
#' 			minResponse=0,maxResponse=365,
#' 			runMAP=TRUE,processExtremes=TRUE,N=500,
#' 			partitionDataForPriors=FALSE,
#' 			maxDiv=0,setStringent=TRUE,
#' 			priorLevel=2
#' 			)
#' #
#' #SUMMARIZE BDA RUN
#' print(stanResultBDA.NC$sample, max_rows = 15)
#' #
#' #_________________________________________________________________________________________________________________________________________________________
#' #
#' #PO DATA
#' responseData.NC = simulated_data$observedTime[simulated_data$stage==2]
#' #
#' \bold{#STAN RUN PO - type 'intercept-only'}
#' stanResult.NC = runStanPhenology(
#' 			type="intercept-only",
#' 			responseData=responseData.NC,
#'			hyperparams_noCovariates=hyperparameters,
#' 			minResponse=0,maxResponse=365,
#' 			runMAP=TRUE,processExtremes=TRUE,N=500,
#' 			partitionDataForPriors=FALSE,
#' 			maxDiv=0,setStringent=TRUE,
#' 			priorLevel=2
#' 			)
#' #
#' #SUMMARIZE PO RUN
#' print(stanResult.NC$sample, max_rows = 15)
#' ##########################################################################################################################################################
#' }
runStanPhenology = function(type=c("intercept-only","interval","interval-full","full"), responseData=NULL, stage=NULL, hyperparams_noCovariates=NULL, onsetCovariateData=NULL, durationCovariateData=NULL, onsetHyperBeta=NULL, onsetHyperAnchor=NULL, durationHyperBeta=NULL, durationHyperAnchor=NULL, sigmaHyper=NULL, minResponse=0, maxResponse=365, maxDiv=0, setStringent=TRUE, runMAP=TRUE, processExtremes=TRUE, N=500, keepScale=FALSE, partitionDataForPriors=FALSE, maximizeSampleSize=FALSE, byPassChecks=FALSE,priorLevel=2, threshApprox=NULL, ...) {

  ## ###########################################################################
	## CHECK STAN BLOCK
  if( byPassChecks == FALSE ){
    passed <- makeSTANpassChecks() # Quick silent check.
    if( passed == FALSE ){
      ## Run the interactive check that can fix the dependencies:
      ## This function can stop the code if dependencies are not detected.
      makeSTANchecks()
    }
  }
  ## ###########################################################################

  type = match.arg(type)
	cat("Running a Stan analysis.\n")

	if(!is.vector(responseData)) {
		stop("Expecting a vector of real numeric values of the collection times (e.g., day of year (DOY) of specimens in the phenophase.")
	}
	if(length(responseData)<10) {
		warning("Ten or fewer data items is a very small sample size and will likely result in inaccurate inferences and a high divergence rate during Bayesian inference. The sample size should be at least 60, especially when covariates are used.")
	}

	if(!(type=="intercept-only" || type=="full" || type=="interval-full" || type=="interval") ) {
		cat(paste("Unsupported type: ", type, "\nType should be 'intercept-only' or 'interval' or 'full' or 'interval-full'.\n"))
		stop("Unsupported type error.")
	}

	if(partitionDataForPriors) {
		cat("The data will be partitioned into two sets with 30% and 70% of the data. \n\n30% will be used to carry out a preliminary analysis using quantiles to estimate the prior distribution hyperparameter values. \n\n70% of the data will be used to carry out a Stan Bayesian analysis to obtain the posterior distributions of parameters.\n\nIf other hyperparameter information was provided as input, it will be ignored. \n\nThe calculated values based on quantiles are approximate; you may need to use other sources of data to get better estimates of prior hyperparameter values, especially if the Stan run results in divergences or other poor diagnostics.\n\n")

		prop = 0.3
		scale = (maxResponse-minResponse) / (365 - 0)

		if(type=="intercept-only") {
			warning("Automated hyperparameters are not recommended for data without covariates. Estimates of duration are likely to be inaccurate.")
			if(partitionDataForPriors) {
			partition = partitionResponseData(responseData = responseData, prop = prop)
			responseData = partition$dataForInference
			hyperparams_noCovariates = getHyperparametersViaQuantiles(responseDataForPrior = partition$dataForPrior, scale = scale)
			}
		}
		else if(type=="interval-full" || type=="interval") {
			if(partitionDataForPriors) {
				stop("Partitioning data for priors is not supported for interval estimation.")
			}
		}
		else if(type=="full") {
			if(partitionDataForPriors) {
			#partition data
			partition = partitionResponseCovariateData(responseData=responseData, onsetCovariateData=onsetCovariateData, durationCovariateData=durationCovariateData, prop=prop)

			#get the data for inference
			#if(length(partition$responseDataForInference)<150 && maximizeSampleSize) {
			if(maximizeSampleSize) {
				   #don't reduce the amount of data - this is statistically invalid because the same data are used to estimate the prior hyperparameters and carry out the full Bayesian analysis, but when the sample size is too small, the error on the estimates will be unacceptably large unless a fuller dataset is used.
			}
		else {
			responseData = partition$responseDataForInference
			onsetCovariateData = partition$onsetCovariateDataForInference
			durationCovariateData = partition$durationCovariateDataForInference
		}

			#get the data for prior and set the prior hyperparameters
			prior = getHyperparametersViaQuantileRegression(responseDataForPrior=partition$responseDataForPrior, onsetCovariateDataForPrior=partition$onsetCovariateDataForPrior, durationCovariateDataForPrior=partition$durationCovariateDataForPrior, lowerQuantile=0.1, upperQuantile=0.9)

			#set the prior hyperparameters
			onsetHyperBeta = prior$onsetHyperBeta
			onsetHyperAnchor = prior$onsetHyperAnchor
			durationHyperBeta = prior$durationHyperBeta
			durationHyperAnchor = prior$durationHyperAnchor
			#cessationHyperAnchor = prior$cessationHyperAnchor
			sigmaHyper = prior$sigmaHyper
		}
		}
	}

	if(type == "intercept-only" || type == "interval") {
		cat("No covariates will be included in this analysis.\n")
		if(sum(is.na(hyperparams_noCovariates)) || length(hyperparams_noCovariates) != 6) stop("Expecting six hyperparameter values (mean and sd for mean onset, mean and sd for mean duration, mean and sd for sigma. Or, if you want hyperparameter values to be estimated for you, set 'partitionDataForPriors' to TRUE.")
		if(type == "intercept-only") {
			return(runStan.NoCovariates.T.GP(fileOrData=responseData, minResponse=minResponse, maxResponse=maxResponse, hyperparameters = hyperparams_noCovariates, dataProvided=TRUE, runMAP=runMAP, setStringent=setStringent, maxDiv=maxDiv, processExtremes=processExtremes, N=N, threshApprox=threshApprox, ...))
		}
		else if(type=="interval") {
			if(keepScale) {
				return(runStan.NoCovariates.Interval.T.GP(responseData=responseData,stage=stage, minResponse=minResponse, maxResponse=maxResponse, hyperparameters = hyperparams_noCovariates, dataProvided=TRUE, runMAP=runMAP, keepScale=TRUE, setStringent=setStringent, maxDiv=maxDiv, processExtremes=processExtremes, N=N, threshApprox=threshApprox, ...))
			}
			else {
				return(runStan.NoCovariates.Interval.T.GP(responseData=responseData,stage=stage, minResponse=minResponse, maxResponse=maxResponse, hyperparameters = hyperparams_noCovariates, dataProvided=TRUE, runMAP=runMAP, keepScale=FALSE, setStringent=setStringent, maxDiv=maxDiv, processExtremes=processExtremes, N=N, threshApprox=threshApprox, ...))
			}
		}
	}
	else if(type == "full" || type == "interval-full") {
		#cat("Checking conditions to run a Stan analysis with covariates A.\n")
		if(sum(is.na(onsetCovariateData)) || !is.data.frame(onsetCovariateData) ) {
			cat("Please remove all NA values, and provide a data frame of the onset model covariate (predictor variable) data. \nThese might be temperature or precipitation data at the specimen collection sites, for example.\nEach column of the data frame should be named with the predictor variable name (e.g. 'meanAnnualTemperature').\nThe order of the rows should correspond to the order of the elements in the response variable data vector (i.e., the first element in the response vector corresponds with the first row in the onset covariate data frame, the second element with the second row, and so forth).\nPlease see documentation for additional information and examples.")
			stop("Please provide appropriate inputs")
		}
		#cat("Checking conditions to run a Stan analysis with covariates B.\n")
		if(sum(is.na(durationCovariateData)) || !is.data.frame(durationCovariateData) ) {
			cat("Please remove all NA values, and provide a data frame of the duration model covariate (predictor variable) data. \nThese might be temperature or precipitation data at the specimen collection sites, for example.\nEach column of the data frame should be named with the predictor variable name (e.g. 'meanAnnualTemperature').\nThe order of the rows should correspond to the order of the elements in the response variable data vector (i.e., the first element in the response vector corresponds with the first row in the duration covariate data frame, the second element with the second row, and so forth).\nPlease see documentation for additional information and examples.")
			stop("Please provide appropriate inputs")
		}
		#cat("Checking conditions to run a Stan analysis with covariates C.\n")
		if(!is.data.frame(onsetHyperBeta) || !is.data.frame(durationHyperBeta) || !is.vector(onsetHyperAnchor) || !is.vector(durationHyperAnchor) || !is.vector(sigmaHyper)) {
			cat("Expecting the following:\n\tonsetHyperBeta:\n\t\tA data frame with two columns. The first column is the mean hyperparameter for the onset slope coefficient for each covariate. The second column is the standard deviation of the onset slope coefficient for each covariate. The first row in the hyperparameters data frame is the first covariate corresponding to the first column in the covariate file, the second row with the second covariate, and so forth. \n\tdurationHyperBeta:\n\t\tA data frame with two columns. The first column is the mean hyperparameter for the duration slope coefficient for each covariate. The second column is the standard deviation of the duration slope coefficient for each covariate. The first row in the hyperparameters data frame is the first covariate corresponding to the first column in the covariate file, the second row with the second covariate, and so forth. \nonsetHyperAnchor:\n\t\tA two-element vector with the mean of the prior and the standard deviation of the prior for the onset anchor (the mean onset value when no covariate data are included).\n\tdurationHyperAnchor:\n\t\tA two-element vector with the mean of the prior and the standard deviation of the prior for the duration anchor (the mean duration value when no covariate data are included).\n\tsigmaHyper:\n\t\tA two-element vector with the mean of the prior and the standard deviation of the prior for the sigma model parameter (variation in onset times and variation in cessation times).\nSee documentation for additional information and examples")
			stop("Please provide appropriate inputs")
		}
		cat(paste0(ncol(onsetCovariateData), " onset covariates and ", ncol(durationCovariateData), " duration covariates will be included in this analysis. Continuing.\n"))
		if(type=="full") {
			#if(!is.vector(cessationHyperAnchor)) {
				#cat("Cessation anchor mean and standard deviation hyperparameter values are needed for the full model and should be provided as a two-element vector.\n")
				#stop("Please provide appropriate inputs")
			#}
			cat("Calling specialized runStan functions for presence-only data.\n")
		  ## DANIEL: Changed argument response to responseData in runStan.WithCovariates.T.GP
		return(runStan.WithCovariates.T.GP(responseData=responseData, minResponse=minResponse, maxResponse=maxResponse, onsetCovariateData=onsetCovariateData, durationCovariateData=durationCovariateData, onsetHyperBeta=onsetHyperBeta, onsetHyperAnchor=onsetHyperAnchor, durationHyperBeta=durationHyperBeta, durationHyperAnchor=durationHyperAnchor, sigmaHyper=sigmaHyper, setStringent=setStringent, dataProvided=TRUE, priorLevel=priorLevel))
		}
		else if(type == "interval-full") {
			cat("Calling specialized runStan functions for presence-absence data.\n")
			if(length(stage) != length(responseData))
			{
				stop("The vector of stages must be the same length as the vector of response data. For each individual sampled, there should be a response time and there should be a stage associated with the time. If the time is before the phenophase started, the stage is 1, if during the phenophase, the stage is 2, and if after the phenophase ended, the stage is 3.")
			}

			return(runStan.WithCovariates.Interval.T.GP(responseData=responseData,stage=stage, minResponse=minResponse, maxResponse=maxResponse, onsetCovariateData=onsetCovariateData, durationCovariateData=durationCovariateData, onsetHyperBeta=onsetHyperBeta, onsetHyperAnchor=onsetHyperAnchor, durationHyperBeta=durationHyperBeta, durationHyperAnchor=durationHyperAnchor, sigmaHyper=sigmaHyper, setStringent=setStringent, dataProvided=TRUE, priorLevel=priorLevel))
		}
	}
}

