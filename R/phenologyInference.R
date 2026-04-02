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
#' @noRd
NULL

#' @rdname expectation_functions
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
NULL

#' @rdname pi_functions
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
NULL

#' @rdname sd_functions
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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

#' Run a Bayesian analysis of phenological events using biocollection data (Core Package Function!)
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
#' Users can also optionally choose to analyze phenological extremes. To do so requires an estimate of the population size from which biocollection records were obtained. The concept of "population" in this case can be quite loose and can be the metapopulation consisting of all populations from which samples were selected. 
#' 
#' The parameters to be estimated for the 'full' model include the marginal mean onset time (\eqn{\mu_O}), the marginal mean duration (\eqn{\mu_D}), one parameter (\eqn{\beta_O} coefficient) for each covariate used to model the mean onset time, one parameter (\eqn{\beta_D} coefficient) for each covariate used to model the mean duration time, and one \eqn{\sigma_O} parameter which models both the standard deviation of the distribution of onset times and the standard deviation of the distribution of durations. If there are \eqn{K_O} onset covariates, and \eqn{K_D} duration covariates, then there is a total of 3 + \eqn{K_O} + \eqn{K_D} free parameters (including the two \eqn{\alpha} intercept terms and \eqn{\sigma_O}). 
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
#' There are two main types of datasets that can be analyzed. When only observation times of individuals in the stage of interest are available, these are presence-only data (PO). When observation times are available for multiple stages of interest, these are multistage data. Theory for these two types of datasets differs.
#'
#' The following sections can safely be skipped if not interested in theory. 
#' 
#' \bold{NOTATION}
#' 
#' Both visible and latent events are modeled with random variables (RVs). 
#' 
#' RVs are symbolized with the first letter of the state of interest: Onset (O), Duration (D), Cessation (C), observed collection Times (T). 
#' 
#' RVs representing extremes are subscripted. For example, first onset in a population is \eqn{O_1}, and last cessation is \eqn{C_N} where \eqn{N} is the population size. 
#' 
#' Parameters are subscripted by the associated random variable. For example, mean duration is \eqn{\mu_D} and the standard deviation of the last cessation is \eqn{\sigma_{C_N}}. In the linear models, intercept parameters are represented by \eqn{\alpha} and slope coefficients are \eqn{\beta}. So, for the model of mean onset, the intercept is \eqn{\alpha_O} and a slope coefficient for the \eqn{i^{th}} covariate is \eqn{\beta_{O,i}}.
#' 
#' The model implemented in Stan is based on the properties of a degenerate Gaussian process that samples lines with intercepts whose values are normally distributed. When a line passes a first threshold at \eqn{\mu_O}, the process begins, and when that line passes a second threshold at \eqn{\mu_D}, the process ends. This process results in normally-distributed onsets: \eqn{O \sim N(\mu_O, \sigma^2_O)}, a constant duration \eqn{D}, and a normally-distribed cessation: \eqn{C \sim N(\mu_O + D, \sigma^2_O)}.
#' 
#' The process of interest occurs in a time period (TP), such as a year for cyclical annual processes, or a day for cyclical diel processes, that has a minimum time of \eqn{m=0} and a maximum time of \eqn{M}. For a year, the day of year time is between \eqn{m=0} and \eqn{M=365}. 
#' 
#' The dataset consists of a population of individuals undergoing the process of interest contained within the TP. Onset and duration times can vary for each individual in the population. 
#' 
#' Onset times have intrinsic noise due to \eqn{\sigma_O}, whereas durations do not. Both onset times and duration lengths can vary in relation to internal states such as energy store level or organismal size, or to external conditions such as climate variation (temperature summaries for spring, etc.) or environmental variation (soil type, moisture levels, light levels, etc.). Any number of covariates can be used with standard cautions about overly complex models (overfitting, model identifiability, parameter uncertainty, collinearity, etc.).
#' 
#' The sample design is as follows. Individuals are randomly sampled from a population at random points in time during the time period.  
#' 
#' The only RVs with directly observed states are the observed times (\eqn{T}) and the observed phenological stages (\eqn{S}). All the other RVs are latent (hidden) and estimated from the observed times and associated stages.
#' 
#' For the theory, we simplify things by examining a time period (TP) starting at time 0 and of length 1 time unit. Data sampled from other length time periods can be translated to start at time 0 by subtracting \eqn{m} from all times and scaled to a length of 1 by dividing all sampled times by \eqn{M-m}. In what follows, we assume such transformations have already been made. These transformations are handled automatically by this runStanPhenology function. If users wish to transform their data in other ways prior to analysis, we leave this up to the user. 
#' 
#'  
#'  \bold{THEORY: MULTISTAGE}
#'  
#'  For \bold{multistage} data, all observed times are kept and the stage (e.g., flowering, fruiting, active growth, dormancy) is recorded with an integer label.
#'  For the core theory, we need to determine the probabilities of each stage given the sampled time and model parameters. The core parameters are the mean onset of the first full stage in the TP, the mean durations, and the variation in onset times. 
#'
#'  Here, the notation is slightly different that that presented above. \eqn{O_1} represents the onset time of the first stage, rather than the first onset phenological extreme.
#'
#'  If you are interested in doing an analysis of data before a phenophase of interest, during a phenophase of interest, and after a phenophase of interest (BDA analysis), code the data with the following stages: 1 = before, 2 = during, and 3 = after.
#'  
#'  The probability that an individual is in a stage that is before the first full stage begins (ie, stage 1) is \deqn{P_{1|j}(O_1>t_j|\theta)=1-\Phi(z_{1|j})} where \eqn{t_j} is the observed time of collection of the \eqn{j^{th}} individual, \eqn{\Phi(z) \equiv P(Z \le z)} is the cumulative distribution function of a standard normal random variable, and z values are the standardized times \eqn{z = \frac{t-\mu}{\sigma}}. 
#'  
#'  The probability that an individual \eqn{j} is in stage \eqn{s} is the probability that the sampled time is after the onset of stage \eqn{s} and before the onset of the next stage \eqn{s+1}. Using De Morgan's law for probabilities, \deqn{P_{s|j}(O_s \leq t_j, O_{s+1} > t_j|\theta) = \Phi(z_{s|j})-\Phi(z_{s+1|j})}
#'  
#'  Last, the probability that an individual is in the last stage is the probability that the sampled time is after the onset of the last stage S: \deqn{P_{S|j}(O_S \leq t_j|\theta) = \Phi(z_{S|j})}
#'  
#'  The likelihood follows a categorical random variable: \deqn{L(\theta) = \prod_{j=1}^{j=N}\prod_{s=1}^{s=S}P_{s|j}^{I(s_j = s)}} where the indicator function \eqn{I(s_j=s)} returns 1 whenever the stage for the \eqn{j^{th}} individual equals \eqn{s} and returns 0 otherwise.
#'  
#'  \bold{THEORY: OVERLAPPING MULTISTAGE}
#'  
#'  Please note that implementation of the overlapping multistage model is in progress.  You may run the model and expect fairly accurate estimates, but the uncertainty estimates are unlikely to be accurate. Please use with discretion.
#'
#'  The multistage theory above is a special case of a more general multistage theory in which each individual produces one or more phenological units that can be at different developmental stages on the same individual. For example, a plant may produce multiple floral units in a season. In this general case, the categorical model for the single individual in one of multiple stages generalizes to a multinomial model where counts are provided for each stage.
#'
#'  For \bold{individuals with multiple units of organization each with multiple stages}, such as an individual plant with multiple floral units with stages 'bud', 'flower', 'immature fruit', and 'mature fruit', the counts of units in each stage are provided as input along with the time the individual was observed. 
#'  
#'  The probability that a unit is in a particular stage \eqn{s} given the time of observation, \eqn{P_{s|j}}, are the same probabilities as for the multistage theory above. 
#'  
#'  In principle, there are many options to model the distribution of stages of units within an individual at a particular time of observation. One could model the developmental dynamics of the units and use this model to generate the distribution of stages. This is ideal, but requires detailed knowledge of the developmental dynamics of individuals. As an idealization, here, the probability mass function of the stages of units within an individual is modeled as multinomial, which models each unit as developing independently of the other units within the individual once the unit itself has begun development (The units are synchronized with noise within an individual). The stage probabilities \eqn{P_{s|j}} provide the probabilities of occuring in each stage for the multinomial model. 
#'
#'  We noted that this overlapping multistage model is a generalization of the multistage model presented above. However, there are stages when units are not yet present (pre) or when units have abscised (post) from the individual. In such cases, it is impossible to count the number of units in these latent pre and post stages. Additional information in the way of total expected units produced by the individual is therefore needed when latent pre and post stages are present. 
#'
#'  As a generalization of the multistage model, a special case of the overlapping multistage model reduces to the multistage model. When an individual is the single unit of phenological activity, then a 0 is scored for all stages that the individual is not in, and a 1 is scored for the stage the individual is in. In this context, the multinomial probability reduces to the categorical probabilties defined for the multistage analysis above. 
#' 
#'  \bold{THEORY: PRESENCE-ONLY (PO)}
#' 
#'  For \bold{presence-only (PO)} data, only the observed times \emph{during} the process of interest are recorded and all other times are discarded. 
#' 
#'  The theory for presence only is a bit more complicated than the theory for multistage data because marginalization is required to normalize the probability of being in a stage. The ultimate goal for PO analysis is to infer the parameters that describe the distribution of observed sampled times. The other distributions of interest, including extremes, onsets, peak, and cessation, are functions of these parameters (with a caveat for extremes). 
#' 
#'  The distribution whose parameters we want to infer is the probability density of the observed times given the individual was sampled during the phenological stage of interest (stage \eqn{s_j=s}) and given the parameters of the onset distribution and duration: \deqn{p_{j|s} \equiv p(t_j|s_j=s,D,\mu_O,\sigma_O)} For simplicity, let \eqn{\theta \equiv (D,\mu_O,\sigma_O)}. Using Bayes rule, \deqn{p(t_j|s_j=s,\theta) = \frac{P(s_j=s|t_j,\theta)p(t_j|\theta)}{P(s_j=s|\theta)}}
#' 
#'  The probability that a sampled individual \eqn{j} is in the focal stage of interest, \eqn{s}, given the sampled time and parameters for individual \eqn{j}, was defined above as \eqn{P_{s|j}}. 
#'
#' The sampling strategy defines \eqn{p(t_j|\theta)}. In the current implementation, the sample of times is assumed to be a random sample and independent of \eqn{\theta}, i.e., \eqn{p(t_j|\theta)=\frac{1}{M-m}}. Since all observed times are scaled to a TP of length 1 already, \eqn{p(t_j|\theta)=1}. 
#' 
#'  Finally, \eqn{P(s_j=s|\theta)} is the expected length of the stage of interest when the TP is length 1 and the stage occurs entirely within the TP. By the law of total probability, this is \deqn{\int_0^1 P(s_j=s,t_j|\theta)dt_j=\int_0^1 P(s_j=s|t_j,\theta)p(t_j|\theta)dt_j=\int_0^1 P(s_j=s|t_j,\theta)dt_j=\int_0^1P_{s|j}dt_j=\int_0^1 (\Phi(\frac{t-\mu_O}{\sigma_O})-\Phi(\frac{t-D-\mu_O}{\sigma_O}))dt_j} Hearn et al. (XXXX) in their Supporting Information Notes S1 provide a solution for this integral in terms of \eqn{\Phi} and \eqn{\phi} (the density function for the standard normal).
#'  
#'  Defining the normalization constant, \eqn{c_j}, as \eqn{\frac{P(s_j=s|\theta)}{p(t_j|\theta)} = P(s_j=s|\theta)}, then the likelihood function for presence-only data is: \deqn{\prod_{j=1}^{j=N} \frac{p_{j|s}}{c_j}} In this case, the normalization constant is different for each individual once predictors for mean onset and duration are defined in terms of covariates, as each individual has different covariate values. So, it doesn't cancel. The next section develops this aspect of the model. 
#' 
#' \bold{THEORY: MODELING COVARIATES}
#' 
#' In the current implementation, the mean onset time and the duration times are modeled as linear functions of covariates. Different covariates can be used for the onset model and for the duration model(s), but for the current implementation of multistage inference, the covariates must be the same for onset and duration model(s). If we let \eqn{\vec{X}} represent the vector of covariates and \eqn{\vec{\beta}} represent the vector of coefficients, then the following are the models for mean onset and mean duration: \deqn{O \sim N(\alpha_O + \vec{\beta_O} \cdot \vec{X}, \sigma_O^2)} and \deqn{D \equiv \mu_D = \alpha_D+\vec{\beta_D} \cdot \vec{X}} For multistage inference, there are multiple duration models, one for each stage, each with a different intercept \eqn{\alpha} and coefficients \eqn{\vec{\beta}}. 
#' 
#' These models of mean onset and duration make the model hierarchical, so that with the inclusion of covariates, everywhere in the likelihood function where \eqn{\mu_O} or \eqn{D} appear, they are replaced by the above linear functions. If there are \eqn{K_O} covariates for the onset model and \eqn{K_D} (\eqn{K_D} = \eqn{K_O} for multistage models) covariates for the duration model(s), along with the two \eqn{\alpha} intercepts and the \eqn{\sigma_O} standard deviation parameter, the fully specified model has \eqn{2 + K_O + K_D + 1} estimable parameters, or S + S*K_O + 1 parameters for the multistage model.
#' 
#' \bold{THEORY: BAYESIAN ANALYSIS AND PRIORS}
#' 
#' In a Bayesian analysis, the probability density of the parameters given the data (which is called the posterior distribution) is proportional to the product of the likelihood (which is the probability of the data given the parameters) and the prior distribution (which is the marginal probability of the parameters): \eqn{Posterior \propto Likehood * Prior}. 
#' 
#' The above sections provided derivations of the likelihood functions for both PO datasets and multistage datasets, so the only remaining thing needed to carry out a Bayesian analysis is to specify the prior distributions. The priors are defined for the \eqn{D,\beta_O,\alpha_O, \beta_D, \alpha_D, \sigma_O} parameters. In the current implementation, the priors can be flat, or they can be normally distributed. 
#' 
#' Set the \emph{priorLevel} function parameter (described below) to 0 for a flat prior. Identifiability for PO analyses necessitates an informative prior.  Priors for the multistage model are automatically set to weakly informative priors if none are provided by the user.
#' 
#' Set the \emph{priorLevel} to any value greater than or equal to 1 for normally distributed priors. A normally distributed prior requires two hyperparameters for each parameter: the mean of the parameter and the standard deviation of the parameter. Users can manually input these hyperparameter values or can use automated prior specification. 
#' 
#' Normal priors are used primarily for convenience and interpretability. Under standard regularity conditions, the Bernstein–von Mises theorem implies that as sample size increases, the posterior distribution becomes approximately multivariate normal, centered near the maximum likelihood estimate, with covariance given by the inverse Fisher information scaled by \eqn{1/n}. In this large-sample regime, posterior uncertainty is therefore well summarized by normal approximations and depends only weakly on the particular prior, provided the prior assigns positive density in a neighborhood of the true parameter value. This asymptotic behavior makes normal priors a natural default choice in many settings. Moreover, in sequential Bayesian analyses where researchers use a previous large-sample posterior as the prior for a subsequent analysis, the inherited prior will itself be approximately normal.
#' 
#'
#' @param type The type of model to be implemented. Five options are available. Option 'full' will model mean onset and mean duration as a linear function of covariates for presence-only data, 'intercept-only' will model the distribution of onset times and the distribution of durations with no consideration of covariates for presence-only data. Option 'multistage-full' includes randomly sampled times of individuals, their phenological stage, and covariate data, with one linear model per stage, including for the interval between the start of a time period and the first stage in that time period. Option 'multistage-overlap-full' is the most general option. This option models individuals with 1 or more units that can occur in different stages. For example, an individual plant can have multiple floral units in different stages: bud, open flower, immature fruit, mature fruit. There is no multistage option for no covariates, but if you want to run such a model, set up a dummy covariate with a constant value and ignore the parameters for the linear model in the output. The Examples section below provides example usage of each scenario (default: intercept-only)
#' @param responseData A vector of the response data. The response data are simply the observed times of collection, often the day of year for biocollection data. This vector is prepared by the 'preparePhenologyData' function, or advanced users can prepare it using whatever tools they choose. The data should be in the original scale unless advanced users have a reason to use preprocessing transformations. (default: NULL)
#' @param stage Only for use with type 'multistage-full'. A vector of integers the same length as the number of observations. For each observation, the stage is the numeric value associated with each stage. The stage is 1 if the observation is in the interval between the start of the time period and before the first phenophase in the time period. (default: NULL)
#' @param stageCounts Only for use with type 'multistage-overlap-full'. A matrix of intergers with columns representing stages and rows representing individuals. Each row represents an individual, and the counts of units (e.g. floral units) in each stage corresponding to the columns is provided. If the stage is before a unit has developed or after units have senesced, provide a count for that stage indicating the average number of units per individual when the units are active (e.g., expected number of flowers for the plant). If the individual represents the unit and stages are sequential (non-overlapping within the unit), then there is one unit, and the counts should be 0 for all stages, except for the stage the individual is in, which should have a count of 1. For example, if each season, the individual produces a single flower and that flower is open, then the stage counts for pre-reproductive, in bud, in flower, in immature fruit, in fruit, post-reproductive would be: 0, 0, 1, 0, 0, 0. This will perform the identical analysis as type 'multistage-full'.
#' @param hyperparams_noCovariates For use with 'intercept-only' models. A vector of six elements that represent the mean and standard deviation (sd) hyperparameter values for the prior distributions of the mean onset, the mean duration, and the sigma parameters, in that order. Items in the vector can be named or unnamed. (default: NULL)
#' @param onsetCovariateData For use with '*full' models. A data frame with each column representing a covariate and each row representing a specimen observation. The number of rows must match the number of items in the responseData vector, and specimens are in the same order as they are in the response data vector. The function 'preparePhenologyData' can prepare this data frame, or advanced users can prepare the data frame using their own tools. Data should be in the original scale, or advanced users can transform the data as deemed appropriate. (default: NULL)
#' @param durationCovariateData For use with '*full' models. A data frame with each column representing a covariate and each row representing a specimen observation. The number of rows must match the number of items in the responseData vector, and specimens are in the same order as they are in the response data vector. In the 'multistage-full' model, the same covariates are used for onset and duration models. The function 'preparePhenologyData' can prepare this data frame, or advanced users can prepare the data frame using their own tools. Data should be in the original scale, or advanced users can transform the data as deemed appropriate. (default: NULL)
#' @param onsetHyperBeta For use with '*full' models. A data frame with two columns to provide hyperparameter values for the covariate slope coefficients for the onset model. The first column is the the mean hyperparameter, and the second is the standard deviation parameter. Each row corresponds with a separate covariate. The first row in the data frame corresponds to the covariate in the first column of the covariate data frame (onsetCovariateData), the second row with the second covariate, and so forth. If there are \eqn{K_O} columns in the covariate data frame, there should be \eqn{K_O} rows in this data frame. Can be left at default when priorLevel is set to 0, in which case Stan default priors are used (not recommended). (default: NULL) 
#' @param onsetHyperAnchor For use with '*full' models. A two-element vector with the mean of the prior and the standard deviation of the prior for the onset anchor (the marginal mean onset time value). Can be left at default when priorLevel is set to 0, in which case Stan default priors are used (not recommended). (default: NULL) 
#' @param durationHyperBeta For use with '*full' models. A data frame with two columns to provide hyperparameter values for the covariate slope coefficients for the duration model. The first column is the the mean hyperparameter, and the second is the standard deviation parameter, or for the 'multistage-full' model, a table with three columns, with the third indicating which stage the parameters are for. Each row corresponds with a separate covariate. The first row in the data frame corresponds to the covariate in the first column of the covariate data frame (durationCovariateData), the second row with the second covariate, and so forth. If there are \eqn{K_D} columns in the covariate data frame, there should be \eqn{K_D} rows in this data frame. Can be left at default when priorLevel is set to 0, in which case Stan default priors are used (not recommended). (default: NULL) 
#' @param durationHyperAnchor For use with '*full' models. A two-element vector with the mean of the prior and the standard deviation of the prior for the duration anchor (the marginal mean phenophase duration value), or for the 'multistage-full' model, a data frame with two columns, the first being the mean and the second being the standard deviation, as before, but each row is for a different stage. Can be left at default when priorLevel is set to 0, in which case Stan default priors are used (not recommended). (default: NULL) 
#' @param sigmaHyper For use with '*full' models. A two-element vector with the mean and standard deviation of the prior distribution for the sigma parameter, or for the 'multistage-full' model, a data frame with two columns, the first being the mean and the second being the standard deviation, as before, but each row is for a different stage. Sigma is the standard deviation of the onset distribution and is also the standard deviation of the cessation distribution. Can be left at default when priorLevel is set to 0, in which case Stan default priors are used (not recommended). (default: NULL) 
#' @param minResponse Minimum value of the response (e.g., day of year); must be set to 0 (i.e., at default) under current implementation (default = 0)
#' @param maxResponse Maximum value of the response (e.g., day of year); typically 365 for Gregorian calendar (default = 365)
#' @param maxDiv The maximum number of divergences to be tolerated during the Stan posterior sampling. This should be 0 unless a biased sample from the posterior is acceptable. (default: 0)
#' @param setStringent Boolean flag to indicate more stringent sampling during Stan runs. Specifically, adapt_delta = 0.99, max_treedepth = 15. Setting to TRUE reduces the chances of divergences, but run time is slower. Usually the difference in run time is not of practical concern, so the default is fine. (default: TRUE)
#' @param runMAP Boolean flag to indicate if Stan should be used to estimate the maximum a posteriori (MAP) values. (default: TRUE)
#' @param processExtremes Boolean flag to indicate if parameters for the phenological extreme distributions should be estimated during the analysis. If set to TRUE, asymptotic estimates of the expected first onset and expected last cessation are sampled. Because asymptotic estimates are provided, estimates are made very quickly. A limitation is that an estimate of the population size is needed for this analysis. (default: TRUE)
#' @param N The population size. Needed when processExtremes is set to TRUE. (default: 500)
#' @param partitionDataForPriors Boolean flag to indicate if methods to automate the specification of prior hyperparameters should be used. In particular, 30% of the data are partitioned into a set to estimate the hyperparameters, and 70% of the data are used for the Bayesian analysis using Stan. If set to TRUE, all other user-provided hyperparameter values (e.g., hyperparams_noCovariates, onsetHyperBeta, onsetHyperAnchor, etc.) are ignored. This is *not* recommended if accurate estimates of duration and onset parameters are needed, but estimates of covariate coefficients tend to be accurate. (default: FALSE)
#' @param maximizeSampleSize Boolean flag for use when partitionDataForPriors is set to TRUE. If a user wants 100% of the data to be used for the Bayesian analysis, set this to TRUE. This is statistically invalid because the same data are used to estimate the prior hyperparameters and carry out the full Bayesian analysis. When the sample size is very small, and no prior information is available, the error on the estimates may be unacceptably large unless the fuller dataset is used for inference. Deprecated and no longer implemented. (default: FALSE)
#' @param byPassChecks Boolean flag indicating whether to bypass checks for Stan. If Stan is not located, an attempt will be made to install Stan. Recommended to keep at default. (default: FALSE)
#' @param priorLevel Specifies types of priors by an integer. Setting to 0 will use flat priors. Above 0 will use the user-input prior hyperparameters with normally distributed priors. For presence-only data, leave at default. (default: 2)
#' @param calculatePPD For use with 'multistage-full' analyses only. Boolean indicating whether the posterior predictive distribution should be estimated. (default: FALSE)
#' @param nXs For use with 'multistage-full' analyses only. Number of increments along the x-axis for posterior predictive plotting. (default: 101)
#' @param nReps For use with 'multistage-full' analysis only. Number of replicates per x-axis increment to sample for posterior predictive estimates. (default: 100)
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
#'
#' In the case of \bold{'intercept-only'} type, a list with the following items is returned:
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
#' 
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
#'                                           )
#' #WRITE EMPIRICAL RESULTS
#' write.csv(stanSummary,"SCanadensis.Empirical.wCovariates.csv")
#' ##
#' ##
#' ##########################################################################################################################################################
#' ##Run more complicated examples, this time with SIMULATED PHENOLOGY and SIMULATED COVARIATES
#' #
#' #COVARIATE MODEL
#' covariate_namesOnset = c("v1", "v2")	
#' covariate_namesDuration = c("v1", "v3")	#Note that for presence-only data for one stage, different covariates can be used
#' cov_all = union(covariate_namesOnset, covariate_namesDuration)
#' #covariance matrix must have all variables (v1, v2, and v3 in this case)
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
#' #	The function simulatePopulationLatentIntervalStates can handle different covariates for onset and for duration of the stage of interest
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
#' #HYPERPARAMETERS - Note these are unbiased. Consider changing them from the true values set above
#' onsetHyperBeta = data.frame(slopesOnset,c(1,1))
#' onsetHyperAnchor = c(150, 5)
#' durationHyperBeta = data.frame(slopesDuration,c(1,1))
#' durationHyperAnchor = c(30, 3)
#' sigmaHyper = c(3,1)
#' #
#' #_________________________________________________________________________________________________________________________________________________________
#' #
#' #PO DATA - Extract data from one stage only
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
#'                                           )
#' ##write results to file
#' write.csv(stanSummaryPO,"Simulated.PresenceOnly.informativePriors.csv")
#' ##
#' ##
#' ##########################################################################################################################################################
#' ##Run a full multistage analysis with SIMULATION PARAMETERS themselves simulated and SIMULATED COVARIATES with MULTIPLE STAGES
#' #
#' library(dplyr)
#' library(posterior)
#' library(ggplot2)
#'  #Set the parameters 
#'  nStages = 4         #Number of stages to be simulated (will be 1 more than this once "unrolled" from cyclical time period
#'  nCovariates = 3     #Number of covariates to be simulated
#'  n = 500             #Sample size
#'  minResponse = 0 		#The minimum observed time should always be at least 0
#'  maxResponse = 365		#Days in the year for yearly cycles
#' 
#'  #Automated - many options available for more control
#'  simulatedData = simulateMultistageData(n=n, nStages=nStages,nCovariates=nCovariates)
#'  
#'  #Plotting 
#'  	#Set which covariate is the x-axis in the plot
#'  targetCovariateIndex = 1
#'  
#'      #Set colors
#'  stageColors = RColorBrewer::brewer.pal(nStages+1, "Set1")
#' 
#'  #Plot the sampled data color-coded by stage 
#'  trueModels = plotMultistageSimulation(simulatedData=simulatedData,
#'                   targetCovariateIndex=targetCovariateIndex,
#'                   stageColors=stageColors,
#'                   drawLatentOnset=FALSE,
#'                   drawObserved=TRUE,
#'                   drawTrueModel=FALSE,
#'                   drawInferredModel=FALSE,
#'                   shadeStage=TRUE,
#'                   minResponse=minResponse, maxResponse=maxResponse)
#'  
#'  #Run inference with Stan - THIS MAY TAKE A WHILE!
#'  stanResults =runStanPhenology(
#'      type="multistage-full",                             #model with many stages and covariates
#'      responseData=simulatedData$outputData$sampledTime,  #observed collection times
#'      stage=simulatedData$outputData$sampledStage,        #observed stage at collection time
#'	    nStages=nStages,				            #number of stages - needed, as not all stages are necessarily sampled
#'      onsetCovariateData=simulatedData$X,                 #covariate data (same as for duration)
#'      maxDiv=4000,					    #should be set to 0, but in case of the stray divergence, set high for this example
#'      calculatePPD=TRUE
#'  )
#'
#' #Plot the posterior predictive distribution - THIS TAKES SOME TIME
#' #BLACK dotted lines = true marginal mean onsets (marginalized across other covariates)
#' #COLORED solid lines = posterior predictive mean onset (no marginalizing across other covariates)
#' #RED dotted lines = mean posterior marginal mean onsets (marginalized across other covariates)
#' p = makeMultistagePosteriorPredictivePlot(stanResult=stanResults,
#'    responseData=simulatedData$outputData$sampledTime,
#'    targetCovariateName="cov1",
#'    covariateData=simulatedData$X,
#'    stageData=simulatedData$outputData$sampledStage,
#'    nStages=nStages,    #with stage coding, an extra stage is added
#'    y_pred=TRUE)          #the predicted responses are calculated when calculatePPD is set to TRUE during Stan run (above)
#'  
#'  #Plot the true mean stage boundaries (based on used simulation parameters) with black dotted lines
#'  for(i in 1:(nStages-1)) {
#'    p = p + geom_abline(intercept = trueModels$trueIntercepts[i], slope = trueModels$trueSlopes[i], color = "black", linewidth = 0.5, linetype="dashed")
#'   }
#'
#'  #Extract basic summary data
#'  probs = c(2.5, 97.5)
#'  measures=c("mean", "median", "sd", "mad")
#'  
#'  summary = print(
#'    stanResults$result$summary(
#'      variables = c(
#'        "sigma","anchor_d", "beta_d", "alpha_d", "anchor_o", "beta_o", "alpha_o"
#'       ),
#'      quantiles = ~ quantile2(., probs=probs/100),
#'      measures
#'    ),
#'    n = Inf
#'  ) %>% as.data.frame()
#'  
#'  
#'   #Set up vectors to store mean durations (means) and mean onsets (means_o) with low and high bounds of credible intervals.
#'  
#'   #Extract slopes and intercepts (anchors) and overlay the inferred lines and 95% BCI's onto earlier plot as colored lines and shaded regions
#'   for(j in 1:(nStages)) {
#'           beta = rep(0,nCovariates)
#'           beta_low = rep(0,nCovariates)
#'           beta_high = rep(0,nCovariates)
#'           varO = paste0("anchor_o[", j, "]")
#'           alpha = summary[summary$variable==varO, "mean"]
#'           alpha_low = summary[summary$variable==varO, paste0("q",probs[1])]
#'           alpha_high = summary[summary$variable==varO, paste0("q",probs[2])]
#'  
#'          for(i in 1:nCovariates) {
#'                   varO = paste0("beta_o[", j, ",", i, "]")
#'                   beta[i] = summary[summary$variable==varO, "mean"]
#'                   beta_low[i] = summary[summary$variable==varO, paste0("q",probs[1])]
#'                   beta_high[i] = summary[summary$variable==varO, paste0("q",probs[2])]
#'           }
#'          model = phenoCollectR:::true_marginal_line(alpha=alpha, beta=beta, mu=simulatedData$covariateMeans, Sigma=simulatedData$Sigma, j=targetCovariateIndex)
#'          model_low = phenoCollectR:::true_marginal_line(alpha=alpha_low, beta=beta_low, mu=simulatedData$covariateMeans, Sigma=simulatedData$Sigma, j=targetCovariateIndex)
#'          model_high = phenoCollectR:::true_marginal_line(alpha=alpha_high, beta=beta_high, mu=simulatedData$covariateMeans, Sigma=simulatedData$Sigma, j=targetCovariateIndex)
#'    p = p + geom_abline(intercept = model$intercept, slope = model$slope, color = "red", linewidth = 0.5, linetype="dashed")
#'   }
#'
#'  print(p)
#'  
#'  
#' ##
#' ##
#' ##
#' ##########################################################################################################################################################
#' ##Conduct an analysis of individuals with multiple phenological units and each unit having one of multiple possible stages.
#' #(For example: an individual plant with multiple flower units, each in a different stage (budding, flowering, fruiting) with pre-reproductive and post-reproductive stages possible with no flowers. 
#' ###SIMULATION PARAMETERS
#'  n = 500            #individuals sampled
#'  meanFlowers = 20   #mean number of flowers per individual
#'
#'  #There can be at most 1 stage before units are visible (pre) or 1 stage after units have abscised, but not both. If there are no stages with visible units, then there can be multiple pre stages only. 
#'  nPre = 1              #number of stages before visible units develop
#'  nVis = 2              #number of stages with visible units
#'  nPost = 0             #number of stages after visible units
#'
#'  nCovariates = 2     #number of simulated environmental factors
#'  nStages = nPre+nVis+nPost     #total number of stages: due to a bug, this must be 3 or greater
#'
#'  maxResponse = 365   #maximum day of year of observation
#'  minResponse = 0     #minimum day of year of observation - must stay at 0
#'
#'  scale = 50    #for transforming the scale within Stan (should be around 10 to 50 - default is fine)
#'
#'  meanOnsetSpread = sample(10:20, size=nStages, replace=TRUE)       #mean weight for each stage duration (Dirichlet)
#'
#'    minStageVariance = (1/10)*(maxResponse*(10/sum(meanOnsetSpread)))^2  #set smallest simulation stage onset standard deviation to expected smallest duration - plenty of overlap!
#'  maxStageVariance = 2*minStageVariance
#'
#'  # Simulate the data
#'  sim = phenoCollectR:::simulateMultistageOverlapData(n=n, meanUnits=meanFlowers, nPre=nPre, nVis=nVis, nPost=nPost, nCovariates=nCovariates, minStageVariance = minStageVariance, maxStageVariance=maxStageVariance, meanOnsetSpread=meanOnsetSpread)
#'
#'  # Plot the simulated data
#'  par(mfrow = c(1, 2))  # 1 row, 2 columns
#'  plotMultistageSimulation(simulatedData = sim, includeOverlap=TRUE, drawInferredModel=FALSE, main="X: times of observation; solid lines: true model of mean onset times", drawLatentOnset=FALSE)
#'
#'  #Get the unit count matrix (e.g., number of floral units in floral stages)
#'  stageCountMatrix <- t(
#'  sapply(sim$simulatedIndividuals, `[[`, "stageCounts")
#'  )
#'
#'  #Get the true total counts for each individual (including missing ones in pre or post)
#'  #Currently treated as data (no prior)
#'  #Not needed if all stages are stages that produce visible units (vis stages)
#'  trueTotCounts = sapply(sim$simulatedIndividuals, `[[`, "nUnits")
#'
#'  #Get the observed times
#'  t_obs <- sapply(sim$simulatedIndividuals, `[[`, "sampledTime")
#'
#'  #Print simulated data
#'  print(cbind(stageCountMatrix,t_obs))
#'
#'  stanResults =runStanPhenology(
#'    type="multistage-overlap-full",       #model with many stages and covariates
#'    responseData=t_obs,  #observed collection times
#'    stageCounts=stageCountMatrix,        #observed stage counts at collection time
#'    preN = nPre,
#'    visN = nVis,
#'    postN = nPost,
#'    onsetCovariateData=sim$simulatedData$X,           #covariate data (same as for duration)
#'    totCounts = trueTotCounts,  #currently using true simulated times
#'    fixedCounts = TRUE,		#not used
#'    minResponse = minResponse,
#'    maxResponse = maxResponse,
#'    scale = scale,	    	#scale for transforming data in Stan - default is usually fine
#'    bypass=T,			#force running of multistage overlap, acknowledging still in development
#'    maxDiv=4000             #should be set to 0, but in case of the stray Stan divergence, set high for this example
#'  )
#'
#'          #Extract basic summary data
#'          stanSummary = phenoCollectR:::getMultistageSummary(stanResult=stanResults)
#'          nDivergent = sum(stanResults$result$diagnostic_summary()$num_divergent)
#'          nMaxTreeDepth = sum(stanResults$result$diagnostic_summary()$num_max_treedepth)
#'
#'          stageColors = rainbow(nStages)
#'          trueModels = plotMultistageSimulation(simulatedData=sim, stageColors=stageColors, includeOverlap=TRUE, shadeStage=FALSE, drawLatentOnset=FALSE, drawTrueModel=FALSE, drawInferredModel=FALSE, main="Dotted red: true model; Dotted black: mean posterior estimate")
#'
#'          targetCovariateIndex = 1
#'          # Add estimated mean stage boundaries
#'	    #Black dotted lines are the mean posterior estimates
#'	    #Red dotted lines are the true model lines
#'          for(i in 2:nStages) {
#'            model = phenoCollectR:::true_marginal_line(alpha=stanSummary$mean_o[i-1], beta=stanSummary$beta_o[i-1,], mu=sim$simulatedData$covariateMeans, Sigma=sim$simulatedData$Sigma, j=targetCovariateIndex)
#'            abline(a=model$intercept, b=model$slope, col="black", lwd=3, lty="dotted")
#'            abline(a=trueModels$trueIntercepts[i-1], b=trueModels$trueSlopes[i-1], col="red", lwd=3, lty="dotted")
#'          }
#'	    print("True models (a = intercepts / stage durations / anchors, b = slopes)")
#'	    print(trueModels)
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
#' simulated_data = simulatePopulation(N=n, mu_O=mean_responseOnset, sigma_O=sdOnset, mu_D_raw=mean_responseDuration) 
#' #
#' #
#' #_________________________________________________________________________________________________________________________________________________________
#' #
#' #HYPERPARAMETERS (biased): 
#' hyperparameters = c(
#'  		#mean onset, sd mean onset
#'  		120, 30, 
#'  		#mean duration, sd mean duration
#'  		35, 5, 
#'  		#mean sigma, sd mean sigma
#'  		5, 2)
#' #
#' #_________________________________________________________________________________________________________________________________________________________
#' #
#' #PO DATA
#' responseData.NC = simulated_data$Ts
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
#runStanPhenology = function(type=c("intercept-only","full","multistage-full"), responseData=NULL, stage=NULL, hyperparams_noCovariates=NULL, onsetCovariateData=NULL, durationCovariateData=NULL, onsetHyperBeta=NULL, onsetHyperBetaMean=NULL, onsetHyperBetaSD=NULL, onsetHyperAnchor=NULL, durationHyperBeta=NULL, durationHyperBetaMean=NULL, durationHyperBetaSD=NULL, durationHyperAnchor=NULL, sigmaHyper=NULL, minResponse=0, maxResponse=365, maxDiv=0, setStringent=TRUE, runMAP=TRUE, processExtremes=TRUE, N=500, keepScale=FALSE, partitionDataForPriors=FALSE, maximizeSampleSize=FALSE, byPassChecks=FALSE,priorLevel=2, threshApprox=NULL, debug=0, ...) {
#runStanPhenology = function(type=c("intercept-only","full","multistage-full", "multistage-overlap-full"), responseData=NULL, stage=NULL, stageCounts=NULL, preN=NULL, visN=NULL, postN=NULL, nStages=NULL, hyperparams_noCovariates=NULL, onsetCovariateData=NULL, nOnsetCovariates=NULL, durationCovariateData=NULL, nDurationCovariates=NULL, onsetHyperBeta=NULL, onsetHyperBetaMean=NULL, onsetHyperBetaSD=NULL, onsetHyperAnchor=NULL, durationHyperBeta=NULL, durationHyperBetaMean=NULL, durationHyperBetaSD=NULL, durationHyperAnchor=NULL, sigmaHyper=NULL, nTotMean=NULL, nTotSD=NULL, minResponse=0, maxResponse=365, maxDiv=0, setStringent=TRUE, runMAP=FALSE, processExtremes=FALSE, N=500, partitionDataForPriors=FALSE, byPassChecks=FALSE, priorLevel=2, debug=0, nXs=101, nReps=100, calculatePPD=FALSE, ...) {
runStanPhenology = function(type=c("intercept-only","full","multistage-full", "multistage-overlap-full"), responseData=NULL, stage=NULL, stageCounts=NULL, preN=NULL, visN=NULL, postN=NULL, nStages=NULL, hyperparams_noCovariates=NULL, onsetCovariateData=NULL, nOnsetCovariates=NULL, durationCovariateData=NULL, nDurationCovariates=NULL, onsetHyperBeta=NULL, onsetHyperBetaMean=NULL, onsetHyperBetaSD=NULL, onsetHyperAnchor=NULL, durationHyperBeta=NULL, durationHyperBetaMean=NULL, durationHyperBetaSD=NULL, durationHyperAnchor=NULL, sigmaHyper=NULL, totCounts=NULL, fixedCounts=FALSE, minResponse=0, maxResponse=365, scale=50, maxDiv=0, setStringent=TRUE, runMAP=FALSE, processExtremes=FALSE, N=NULL, partitionDataForPriors=FALSE, byPassChecks=FALSE, priorLevel=2, debug=0, nXs=101, nReps=100, calculatePPD=FALSE, bypass=FALSE, ...) {

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

	if(type=="multistage-overlap-full" && !bypass) {
		stop("The multistage overlap model is not currently fullly implemented.")
	}

  #Allow users to input one set of covariate data for multistage analyses
  if(type=="multistage-full" || type=="multistage-overlap-full") {
    if(is.null(onsetCovariateData) && !is.null(durationCovariateData)) {
      onsetCovariateData = durationCovariateData
    }
    if(!is.null(onsetCovariateData) && is.null(durationCovariateData)) {
      durationCovariateData = onsetCovariateData
    }
  nDurationCovariates = ncol(durationCovariateData)
  nOnsetCovariates = ncol(onsetCovariateData)
    if(is.null(nStages)) {
      nStages = ncol(stageCounts)
    }
  }

	#Check data to make sure it is all as expected
	checkInput(
		   type=type, 
		   responseData=responseData,
		   onsetCovariateData=onsetCovariateData,
		   nOnsetCovariates=nOnsetCovariates,
		   durationCovariateData=durationCovariateData,
		   nDurationCovariates=nDurationCovariates,
		   stage=stage,
       		   stageCounts=stageCounts,
		   nStages=nStages,
		   preN=preN,
		   visN=visN,
		   postN=postN,
		   minResponse=minResponse,
		   maxResponse=maxResponse,
		   scale=scale,
		   maxDiv=maxDiv,
		   fixedCounts=fixedCounts,
		   N=N,
		   processExtremes=processExtremes
	)

	checkedPriors = checkPriors (
				     type=type,
				     nStages=nStages,
				     preN=preN,
				     visN=visN,
				     postN=postN,
				     responseData=responseData,
				     hyperparams_noCovariates=hyperparams_noCovariates,
				     onsetCovariateData=onsetCovariateData,
				     durationCovariateData=durationCovariateData,
				     onsetHyperBeta=onsetHyperBeta,
				     onsetHyperBetaMean=onsetHyperBetaMean,
				     onsetHyperBetaSD=onsetHyperBetaSD,
				     onsetHyperAnchor=onsetHyperAnchor,
				     durationHyperBeta=durationHyperBeta,
				     durationHyperBetaMean=durationHyperBetaMean,
				     durationHyperBetaSD=durationHyperBetaSD,
				     durationHyperAnchor=durationHyperAnchor,
				     sigmaHyper=sigmaHyper,
				     stageCounts=stageCounts,
				     totCounts=totCounts,
				     #nTotMean=nTotMean,
				     #nTotSD=nTotSD,
				     partitionDataForPriors=partitionDataForPriors,
				     minResponse=minResponse,
				     maxResponse=maxResponse)

	#update any variables that may have been altered by checkPriors function
	responseData=checkedPriors$responseData
	hyperparams_noCovariates=checkedPriors$hyperparams_noCovariates
	onsetCovariateData=checkedPriors$onsetCovariateData
	durationCovariateData=checkedPriors$durationCovariateData
	onsetHyperAnchor=checkedPriors$onsetHyperAnchor
	onsetHyperBeta=checkedPriors$onsetHyperBeta
	onsetHyperBetaMean=checkedPriors$onsetHyperBetaMean
	onsetHyperBetaSD=checkedPriors$onsetHyperBetaSD
	durationHyperAnchor=checkedPriors$durationHyperAnchor
	durationHyperBeta=checkedPriors$durationHyperBeta
	durationHyperBetaMean=checkedPriors$durationHyperBetaMean
	durationHyperBetaSD=checkedPriors$durationHyperBetaSD
	sigmaHyper=checkedPriors$sigmaHyper
	#nTotMean=checkedPriors$nTotMean
	#nTotSD=checkedPriors$nTotSD

	cat("Running a Stan analysis.\n")

	if(type == "intercept-only") {
		cat("No covariates will be included in this analysis.\n")
		cat("Calling specialized runStan function for presence-only data with no covariates.\n")
		return(runStan.NoCovariates.T.GP(fileOrData=responseData,
			 minResponse=minResponse,
			 maxResponse=maxResponse,
			 hyperparameters = hyperparams_noCovariates,
			 dataProvided=TRUE,
			 runMAP=runMAP,
			 setStringent=setStringent,
			 maxDiv=maxDiv,
			 processExtremes=processExtremes,
			 N=N,
			 ...))
	}

	#Analyses with covariate data
	cat(paste0(ncol(onsetCovariateData), " onset covariates and ", ncol(durationCovariateData), " duration covariates will be included in this analysis. Continuing.\n"))

	if(type=="full") {
		cat("Calling specialized runStan functions for presence-only data with covariates.\n")
		## DANIEL: Changed argument response to responseData in runStan.WithCovariates.T.GP
		return(runStan.WithCovariates.T.GP(responseData=responseData,
			   minResponse=minResponse,
			   maxResponse=maxResponse,
			   onsetCovariateData=onsetCovariateData,
			   durationCovariateData=durationCovariateData,
			   onsetHyperAnchor=onsetHyperAnchor,
			   onsetHyperBeta=onsetHyperBeta,
			   durationHyperAnchor=durationHyperAnchor,
			   durationHyperBeta=durationHyperBeta,
			   sigmaHyper=sigmaHyper,
			   setStringent=setStringent,
			   dataProvided=TRUE,
			   priorLevel=priorLevel))
	}
	else if(type == "multistage-full") {
		cat("Calling specialized runStan functions for multistage data.\n")
		return(
		       runStan.WithCovariates.Multistage.durations.GP(
			      responseData=responseData,
			      stage=stage,
			      nStages=nStages,
			      minResponse=minResponse,
			      maxResponse=maxResponse,
			      covariateData=onsetCovariateData,
			      nCovariates=nOnsetCovariates,
			      onsetHyperAnchor=onsetHyperAnchor,
			      onsetHyperBeta=onsetHyperBeta,
			      durationHyperAnchor=durationHyperAnchor,
			      durationHyperBetaMean=durationHyperBetaMean,
			      durationHyperBetaSD=durationHyperBetaSD,
			      sigmaHyper=sigmaHyper,
			      setStringent=setStringent,
			      maxDiv=maxDiv,
			      nXs=nXs,
			      nReps=nReps,
			      calculatePPD=calculatePPD
			      #priorLevel=priorLevel,
			      #processExtremes=processExtremes,
			      #N=N,
		       )
		)
	}
	else if(type == "multistage-overlap-full") {
		cat("Calling specialized runStan functions for within individual overlapping multistage data.\n")
		return(
		       runStan.WithCovariates.Multistage.overlap.GP(
			      responseData=responseData,
			      stageCounts=stageCounts,
			      preN=preN,
			      visN=visN,
			      postN=postN,
			      nCovariates=nOnsetCovariates,
			      minResponse=minResponse,
			      maxResponse=maxResponse,
			      scale=scale,
			      covariateData=onsetCovariateData,
			      onsetHyperBeta=onsetHyperBeta,
			      onsetHyperAnchor=onsetHyperAnchor,
			      durationHyperBetaMean=durationHyperBetaMean,
			      durationHyperBetaSD=durationHyperBetaSD,
			      durationHyperAnchor=durationHyperAnchor,
			      sigmaHyper=sigmaHyper,
			      totCounts=totCounts,
			      fixedCounts=fixedCounts,
			      #nTotMean=nTotMean,
			      #nTotSD=nTotSD,
			      setStringent=setStringent,
			      maxDiv=maxDiv,
			      debug=FALSE
		       )
		)
	}
}
