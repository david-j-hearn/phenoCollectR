#' @export
dO = function(x, mu_O, sigma_O, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") { dO.GP(x, mu_O, sigma_O) }
	else { dO.BB(x, mu_O, sigma_O, minResponse, maxResponse) }
}

#' @export
pO = function(q, mu_O, sigma_O, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") { pO.GP(q, mu_O, sigma_O) }
	else { pO.BB(q, mu_O, sigma_O, minResponse, maxResponse) }
}

#' @export
qO = function(p, mu_O, sigma_O, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") { qO.GP(p, mu_O, sigma_O) }
	else { qO.BB(p, mu_O, sigma_O, minResponse, maxResponse) }
}

#' @export
rO = function(n, mu_O, sigma_O, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") { rO.GP(n, mu_O, sigma_O) }
	else { rO.BB(n, mu_O, sigma_O, minResponse, maxResponse) }
}

#' @export
dOk1 = function(x, N, mu_O, sigma_O, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") { dOk1.GP(x, N, mu_O, sigma_O) }
	else { dOk1.BB(x, N, mu_O, sigma_O, minResponse, maxResponse) }
}

#' @export
qOk1 = function(p, N, mu_O, sigma_O, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") { qOk1.GP(p, N, mu_O, sigma_O) }
	else { qOk1.BB(p, N, mu_O, sigma_O, minResponse, maxResponse) }
}

#' @export
pOk1 = function(q, N, mu_O, sigma_O, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") { pOk1.GP(q, N, mu_O, sigma_O) }
	else { pOk1.BB(q, N, mu_O, sigma_O, minResponse, maxResponse) }
}


#' @export
rOk1 = function(n, N, mu_O, sigma_O, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="GP") { rOk1.GP(n, N, mu_O, sigma_O) }
	else { rOk1.BB(n, N, mu_O, sigma_O, minResponse, maxResponse) }
}

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
