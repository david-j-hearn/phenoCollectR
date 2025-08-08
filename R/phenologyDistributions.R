#' Functions that model phenological distributions.
#'
#' @description Methods to model and simulate phenological distributions. The R convention of 'd' (density or mass function), 'r' (sampling), 'q' (quantile), and 'p' (cumulative distribution function) prefixes the function name. 
#'
#' Eight phenological properties are modeled by random variables. These include the following:
#'
#' O: phenological onset times 
#'
#' D: phenophase durations 
#'
#' C: phenological cessation times
#'
#' Ok1: first onset (order k = 1)
#'
#' CkN: last cessation (order k = N, the population size)
#'
#' T: observed specimen collection times, e.g., day of year for a yearly time period
#'
#' PNt: Proportion of the population of size N in the phenophase at time t; only the unnormalized 'd' and random sampling 'r' functions are currently implemented
#'
#' R: The total time range that at least one individual in a population is in the phenophase; only the 'd' and 'r' functions are currently implemented 
#'
#' The function name is constructed by combining one of 'd', 'r', 'q', and 'p' as a prefix with one of the above variable labels. For example, rT refers to the function that generates random samples from the distribution of observed collection times, T. 
#'
#' Most functions share a subset of input parameters. Provide only the needed parameters for the function. 
#'
#' @rdname dist_family
#' @param x A vector of the response data (e.g. day of year values)
#' @param q The quantile
#' @param p The lower tail area (probability) 
#' @param n The sample size
#' @param t The time at which the proportion of the population in the phenophase is to be evaluated
#' @param mu_O Mean onset time
#' @param sigma_O Standard deviation for the onset time distribution
#' @param mu_D Mean phenophase duration length
#' @param sigma_D Standard deviation for the phenophase duration distribution
#' @param minResponse Minimum value of the response (e.g., day of year); must be set to 0 under current implementation (default = 0)
#' @param maxResponse Maximum value of the response (e.g., day of year); typically 365 for Gregorian calendar (default = 365)
#' @param N The population size
#' @param type The model type, either BB (beta onset, beta duration) or GP (Gaussian process with a shared standard deviation for onset and cessation and a constant duration) (default = "GP")
#'
#' @return A vector, following standard R conventions for the 'd', 'p', 'r', and 'q' functions
#' @examples 
#' \donttest{
#' #Set the sample size
#' n=100000
#' #Set the mean onset time
#' mean_onset = 100
#' #Set the onset time standard deviation, sigma
#' sigma_onset = 10
#' #Set the duration of the phenophase
#' duration = 50
#' #Sample the observed collection times
#' observed_t = rT(n=n, mu_O = mean_onset, sigma_O=sigma_onset, mu_D=duration)
#' #Make a histogram of the observed collection times
#' hist(observed_t, probability=TRUE, xlab="Simulated observed collection times")
#' #Overlay the theoretical curve on the histogram
#' curve(dT(x,mu_O=mean_onset,sigma_O=sigma_onset, mu_D=duration),col="red",add=TRUE)
#'
#' #Or, for a more "extreme" distribution that illustrates the flexibility of the model's shape:
#' #Set the mean onset time
#' mean_onset = 180
#' #Set the onset time standard deviation, sigma
#' sigma_onset = 100
#' #Set the duration of the phenophase
#' mean_duration = 50
#' sigma_duration = 10
#' #Sample the observed collection times
#' observed_t = rT(n=n, mu_O = mean_onset, sigma_O=sigma_onset, mu_D=mean_duration
#'                 , sigma_D=sigma_duration, type="BB")
#' #Make a histogram of the observed collection times
#' hist(observed_t, probability=TRUE
#'      , xlab="Simulated observed collection times (extreme BB model)")
#' #Overlay the theoretical curve on the histogram
#' curve(dT(x,mu_O=mean_onset,sigma_O=sigma_onset, mu_D=mean_duration
#'       , sigma_D=sigma_duration, type="BB"),col="red",add=TRUE,from=0, to=365)
#' }
#' @name phenologyDistributionFunctions
NULL

#' @rdname dist_family
#' @export
dO = function(x, mu_O, sigma_O, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") { dO.GP(x, mu_O, sigma_O) }
	else { dO.BB(x, mu_O, sigma_O, minResponse, maxResponse) }
}

#' @rdname dist_family
#' @export
pO = function(q, mu_O, sigma_O, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") { pO.GP(q, mu_O, sigma_O) }
	else { pO.BB(q, mu_O, sigma_O, minResponse, maxResponse) }
}

#' @rdname dist_family
#' @export
qO = function(p, mu_O, sigma_O, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") { qO.GP(p, mu_O, sigma_O) }
	else { qO.BB(p, mu_O, sigma_O, minResponse, maxResponse) }
}

#' @rdname dist_family
#' @export
rO = function(n, mu_O, sigma_O, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") { rO.GP(n, mu_O, sigma_O) }
	else { rO.BB(n, mu_O, sigma_O, minResponse, maxResponse) }
}

#' @rdname dist_family
#' @export
dD = function(x, mu_O=NA, sigma_O=NA, mu_D, sigma_D=NA, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		dD.GP(x=x, mu_D=mu_D)
	}
	else {
		dD.BB(x=x, mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	}
}

#' @rdname dist_family
#' @export
pD = function(q, mu_O=NA, sigma_O=NA, mu_D, sigma_D=NA, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		pD.GP(q=q, mu_D=mu_D)
	}
	else {
		pD.BB(q=q, mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	}
}

#' @rdname dist_family
#' @export
rD = function(n, mu_O=NA, sigma_O=NA, mu_D, sigma_D=NA, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		rD.GP(n=n, mu_D=mu_D)
	}
	else {
		rD.BB(n=n, mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	}
}

#' @rdname dist_family
#' @export
qD = function(p, mu_O=NA, sigma_O=NA, mu_D, sigma_D=NA, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		qD.GP(p=p, mu_D=mu_D)
	}
	else {
		qD.BB(p=p, mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	}
}

#' @rdname dist_family
#' @export
dOk1 = function(x, N, mu_O, sigma_O, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") { dOk1.GP(x, N, mu_O, sigma_O) }
	else { dOk1.BB(x, N, mu_O, sigma_O, minResponse, maxResponse) }
}

#' @rdname dist_family
#' @export
qOk1 = function(p, N, mu_O, sigma_O, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") { qOk1.GP(p, N, mu_O, sigma_O) }
	else { qOk1.BB(p, N, mu_O, sigma_O, minResponse, maxResponse) }
}

#' @rdname dist_family
#' @export
pOk1 = function(q, N, mu_O, sigma_O, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") { pOk1.GP(q, N, mu_O, sigma_O) }
	else { pOk1.BB(q, N, mu_O, sigma_O, minResponse, maxResponse) }
}

#' @rdname dist_family
#' @export
rOk1 = function(n, N, mu_O, sigma_O, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") { rOk1.GP(n, N, mu_O, sigma_O) }
	else { rOk1.BB(n, N, mu_O, sigma_O, minResponse, maxResponse) }
}

#' @rdname dist_family
#' @export
dC = function(x, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		mu_C = mu_O + mu_D
		sigma = sigma_O
		dC.GP(x, mu_C, sigma)
	}
	else { dC.BB(x, mu_O, sigma_O, mu_D, sigma_D, minResponse, maxResponse) }
}

#' @rdname dist_family
#' @export
pC = function(q, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		mu_C = mu_O + mu_D
		sigma = sigma_O
		pC.GP(q, mu_C, sigma)
	}
	else { pC.BB(q, mu_O, sigma_O, mu_D, sigma_D, minResponse, maxResponse) }
}

#' @rdname dist_family
#' @export
rC = function(n, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		mu_C = mu_O + mu_D
		sigma = sigma_O
		rC.GP(n, mu_C, sigma)
	}
	else { rC.BB(n, mu_O, sigma_O, mu_D, sigma_D, minResponse, maxResponse) }
}

#' @rdname dist_family
#' @export
qC = function(p, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		mu_C = mu_O + mu_D
		sigma = sigma_O
		qC.GP(p, mu_C, sigma)
	}
	else { qC.BB(p, mu_O, sigma_O, mu_D, sigma_D, minResponse, maxResponse) }
}

Pt = function(x, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		mu_C = mu_O + mu_D
		sigma = sigma_O
		Pt.GP(x, mu_O, mu_C, sigma, minResponse, maxResponse)
	}
	else { Pt.BB(x, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) }
}

#' @rdname dist_family
#' @export
dT = function(x, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		mu_C = mu_O + mu_D
		sigma = sigma_O
		dT.GP(x, mu_O, mu_C, sigma, minResponse=0, maxResponse=365, nc=NULL)
	}
	else { dT.BB(x, mu_O, sigma_O, mu_D, sigma_D, minResponse, maxResponse, nc=NULL) }
}

#' @rdname dist_family
#' @export
pT = function(q, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		mu_C = mu_O + mu_D
		sigma = sigma_O
		pT.GP(q, mu_O, mu_C, sigma, minResponse, maxResponse, nc=NULL)
	}
	else { pT.BB(q, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, nc=NULL) }
}

#' @rdname dist_family
#' @export
qT = function(p, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		mu_C = mu_O + mu_D
		sigma = sigma_O
		qT.GP(p, mu_O, mu_C, sigma, minResponse, maxResponse, nc=NULL)
	}
	else { qT.BB(p, mu_O, sigma_O, mu_D, sigma_D, minResponse, maxResponse, nc=NULL) }
}

#' @rdname dist_family
#' @export
rT = function(n, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		mu_C = mu_O + mu_D
		sigma = sigma_O
		rT.GP(n=n, mu_O=mu_O, mu_C=mu_C, sigma=sigma, minResponse=minResponse, maxResponse=maxResponse)
	}
	else { rT.BB(n=n, mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse) }
}

#' @rdname dist_family
#' @export
rCkN = function(n, N, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		mu_C = mu_O + mu_D
		sigma = sigma_O
		rCkN.GP(n, N, mu_C, sigma)
	}
	else { rCkN.BB(n, N, mu_O, sigma_O, mu_D, sigma_D, minResponse = minResponse, maxResponse = maxResponse) }
}

#' @rdname dist_family
#' @export
dCkN = function(x, N, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		mu_C = mu_O + mu_D
		sigma = sigma_O
		dCkN.GP(x, N, mu_C, sigma)
	}
	else { dCkN.BB(x, N, mu_O, sigma_O, mu_D, sigma_D, minResponse=minResponse, maxResponse=maxResponse) }
}

#' @rdname dist_family
#' @export
pCkN = function(q, N, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		mu_C = mu_O + mu_D
		sigma = sigma_O
		pCkN.GP(q, N, mu_C, sigma)
	}
	else { pCkN.BB(q, N, mu_O, sigma_O, mu_D, sigma_D, minResponse=minResponse, maxResponse=maxResponse) }
}

#' @rdname dist_family
#' @export
qCkN = function(p, N, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		mu_C = mu_O + mu_D
		sigma = sigma_O
		qCkN.GP(p, N, mu_C, sigma)
	}
	else { qCkN.BB(p, N, mu_O, sigma_O, mu_D, sigma_D, minResponse=minResponse, maxResponse=maxResponse) }
}

#' @rdname dist_family
#' @export
dR = function(x, N, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		mu_C = mu_O + mu_D
		sigma = sigma_O
		dR.GP(x, N, mu_O, mu_C, sigma, minResponse=minResponse, maxResponse=maxResponse)
	}
	else { dR.BB(x, N, mu_O, sigma_O, mu_D, sigma_D, minResponse, maxResponse) }
}

#' @rdname dist_family
#' @export
rR = function(n, N, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		mu_C = mu_O + mu_D
		sigma = sigma_O
		rR.GP(n, N, mu_O, mu_C, sigma)
	}
	else { rR.BB(n, N, mu_O, sigma_O, mu_D, sigma_D, minResponse, maxResponse) }
}

#' @rdname dist_family
#' @export
dPNt = function(x, t, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		mu_C = mu_O + mu_D
		sigma = sigma_O
		dPNt.GP(x, t, mu_O, mu_C, sigma, minResponse, maxResponse, res = 1000)
	}
	else { dPNt.BB(x, t, mu_O, sigma_O, mu_D, sigma_D, minResponse=minResponse, maxResponse=maxResponse, res = 1000) }
}

#' @rdname dist_family
#' @export
rPNt = function(n, t, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") {
		mu_C = mu_O + mu_D
		sigma = sigma_O
		rPNt.GP(n, t, mu_O, mu_C, sigma, minResponse, maxResponse, res=1000, asPercent=TRUE)
	}
	else { rPNt.BB(n, t, mu_O, sigma_O, mu_D, sigma_D, minResponse, maxResponse, res=1000) }
}
