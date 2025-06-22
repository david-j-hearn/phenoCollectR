#' Distribution functions
#'
#' These are distribution functions with 'd' (density), 'r' (sampling), 'q' (quantiles), and 'p' (probability) options. Most functions share the same input parameters.
#' @rdname dist_family
#' @param x data
#' @param q the quantile
#' @param p the probability
#' @param n the number of samples
#' @param t the parameter t
#' @param mu_O mu 0
#' @param sigma_O sigma 0
#' @param mu_D mu D
#' @param sigma_D sigma D 
#' @param minResponse min of the response
#' @param maxResponse max of the response
#' @param N the population size.
#' @param type the type
#'
#' @returns The results.
#' @export
#'
#' @examples
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

#' @rdname dist_family
#' @export
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
		rT.GP(n, mu_O, mu_C, sigma, minResponse, maxResponse)
	}
	else { rT.BB(n, mu_O, sigma_O, mu_D, sigma_D, minResponse=minResponse, maxResponse=maxResponse) }
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
