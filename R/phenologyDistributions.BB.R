dD.BB = function(x, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	dataO = convert_beta_check(x, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check(x, mu_D, sigma_D, minResponse, maxResponse)
	alpha_s = dataO$alpha_s
	beta_s = dataO$beta_s
	alpha_d = dataD$alpha_s
	beta_d = dataD$beta_s

	return(dD.BB.s(dataO$scaled_x, alpha_s, beta_s, alpha_d, beta_d) / (maxResponse-minResponse))
}

rD.BB = function(n, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, n=n, minResponse=minResponse, maxResponse=maxResponse)
	dataO = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check((maxResponse-minResponse)/2, mu_D, sigma_D, minResponse, maxResponse)
	t_start = rbeta(n, dataO$alpha_s, dataO$beta_s)
	raw_duration = rbeta(n, dataD$alpha_s, dataD$beta_s)
	duration = raw_duration * (1 - t_start)
	return(duration*(maxResponse-minResponse))
}

pD.BB = function(q, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	dataO = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check((maxResponse-minResponse)/2, mu_D, sigma_D, minResponse, maxResponse)
	alpha_s = dataO$alpha_s
	beta_s = dataO$beta_s
	alpha_d = dataD$alpha_s
	beta_d = dataD$beta_s
	return(pD.BB.s(dataO$scaled_x, alpha_s, beta_s, alpha_d, beta_d) / (maxResponse-minResponse))
}


qD.BB = Vectorize(function(p, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse)

					   if (is.na(p) || p < 0 || p > 1) {
							 stop("Error in qD.BB: p values must be in the interval (0,1)")
  warning("p must be between 0 and 1")
  return(NA_real_)
 }

	dataO = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check((maxResponse-minResponse)/2, mu_D, sigma_D, minResponse, maxResponse)
	alpha_s = dataO$alpha_s
	beta_s = dataO$beta_s
	alpha_d = dataD$alpha_s
	beta_d = dataD$beta_s


 q = tryCatch(
			  {
				  minResponse + uniroot(
  function(d) pD.BB.s(d, alpha_s, beta_s, alpha_d, beta_d) - p,
  lower = 0,
  upper = 1,
  tol = 1e-6
 )$root * (maxResponse-minResponse)
			  },
			  error = function(e) {
				     warning("uniroot failed: ", e$message)
				     stop("uniroot failed: ", e$message)
   NA_real_
			  }
			  )

	return(q)
})

dO.BB = function(x, mu_O, sigma_O, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, minResponse=minResponse, maxResponse=maxResponse)
	data = convert_beta_check(x, mu_O, sigma_O, minResponse, maxResponse)
	return(dO.BB.s(data$scaled_x, data$alpha_s, data$beta_s)/(maxResponse-minResponse) )
}

pO.BB = function(q, mu_O, sigma_O, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, minResponse=minResponse, maxResponse=maxResponse)
	data = convert_beta_check(q, mu_O, sigma_O, minResponse, maxResponse)
	return(pO.BB.s(data$scaled_x, data$alpha_s, data$beta_s))
}

#' @importFrom stats qbeta
qO.BB = function(p, mu_O, sigma_O, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, minResponse=minResponse, maxResponse=maxResponse)
	data = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	return(ifelse(p < 0, NA, ifelse(p>1, NA, minResponse + (maxResponse - minResponse) * qbeta(p, data$alpha_s, data$beta_s))))
}

rO.BB = function(n, mu_O, sigma_O, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, n=n, minResponse=minResponse, maxResponse=maxResponse)
	data = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	return(minResponse + rbeta(n, data$alpha_s, data$beta_s)*(maxResponse-minResponse))
}

dOk1.BB = function(x, N, mu_O, sigma_O, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, N=N, minResponse=minResponse, maxResponse=maxResponse)
	data = convert_beta_check(x, mu_O, sigma_O, minResponse, maxResponse)
	return(dOk1.BB.s(data$scaled_x, N, data$alpha_s, data$beta_s)/(maxResponse-minResponse) )
}

#' @importFrom stats qbeta
qOk1.BB = function(p, N, mu_O, sigma_O, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, N=N, minResponse=minResponse, maxResponse=maxResponse)
	data = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	q_fst = qbeta(1 - (1-p)^(1/N), data$alpha_s, data$beta_s)
	return(minResponse + q_fst * (maxResponse - minResponse))
}

#' @importFrom stats pbeta pnorm
pOk1.BB = function(q, N, mu_O, sigma_O, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, N=N, minResponse=minResponse, maxResponse=maxResponse)
	data = convert_beta_check(q, mu_O, sigma_O, minResponse, maxResponse)
	Fx = pbeta(q=data$scaled_x, data$alpha_s, data$beta_s)
	return(1 - (1 - Fx)^N)
}

rOk1.BB = function(n, N, mu_O, sigma_O, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, n=n, N=N, minResponse=minResponse, maxResponse=maxResponse)
	data = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	return(minResponse + replicate(n, min(rbeta(N, data$alpha_s, data$beta_s)))*(maxResponse-minResponse))
}

dC.BB = function(x, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	dataO = convert_beta_check(x, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check(x, mu_D, sigma_D, minResponse, maxResponse)
	return(dC.BB.s(dataO$scaled_x, dataO$alpha_s, dataO$beta_s, dataD$alpha_s, dataD$beta_s)/(maxResponse-minResponse))
}

pC.BB = function(q, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	dataO = convert_beta_check(q, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check(q, mu_D, sigma_D, minResponse, maxResponse)
	pC.BB.s(dataO$scaled_x, dataO$alpha_s, dataO$beta_s, dataD$alpha_s, dataD$beta_s)
}

#' @importFrom stats rbeta
rC.BB = function(n, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, n=n, minResponse=minResponse, maxResponse=maxResponse)
	dataO = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check((maxResponse-minResponse)/2, mu_D, sigma_D, minResponse, maxResponse)
	t_start = rbeta(n, dataO$alpha_s, dataO$beta_s)
	raw_duration = rbeta(n, dataD$alpha_s, dataD$beta_s)
	duration = raw_duration * (1 - t_start)
	t_end = t_start + duration
	return(minResponse + t_end*(maxResponse-minResponse))
}

qC.BB = function(p, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	dataO = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check((maxResponse-minResponse)/2, mu_D, sigma_D, minResponse, maxResponse)
	#creates interpolated qfunc based on n equi-distant data values
	qfunc = make_quantile_function_from_pdf(ddist= function(x){ dC.BB.s(x,dataO$alpha_s,dataO$beta_s,dataD$alpha_s,dataD$beta_s) }, support = c(0, 1), n = 100)
	return(minResponse + qfunc(p) * (maxResponse - minResponse))
}

Pt.BB = function(x, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	dataO = convert_beta_check(x, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check(x, mu_D, sigma_D, minResponse, maxResponse)
	return(Pt.BB.s(dataO$scaled_x, dataO$alpha_s, dataO$beta_s, dataD$alpha_s, dataD$beta_s))
}

dT.BB = function(x, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, nc=NULL) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	if(is.null(nc)) {
		nc = Pt.nc.BB(mu_O, mu_D, minResponse, maxResponse)
	}
	return((Pt.BB(x, mu_O, sigma_O, mu_D, sigma_D, minResponse, maxResponse) / nc) / (maxResponse-minResponse) )
}

pT.BB = function(q, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, nc=NULL) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	dataO = convert_beta_check(q, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check(q, mu_D, sigma_D, minResponse, maxResponse)
	if(is.null(nc)) {
		nc = Pt.nc.BB(mu_O, mu_D, minResponse, maxResponse)
	}
	return(pT.BB.s(dataO$scaled_x, dataO$alpha_s, dataO$beta_s, dataD$alpha_s, dataD$beta_s, nc))
}

#' @importFrom stats uniroot
qT.BB = function(p, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, nc=NULL) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	dataO = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check((maxResponse-minResponse)/2, mu_D, sigma_D, minResponse, maxResponse)
	if(is.null(nc)) {
		nc = Pt.nc.BB(mu_O, mu_D, minResponse, maxResponse)
	}
	minResponse + sapply(p, function(prob) {
		    if (prob <= 0) return(0)
		    if (prob >= 1) return(1)
		    # Root-finding: Find q such that pobserved(q) = prob
			  tryCatch( {
		    uniroot(
			    function(q) pT.BB.s(q, dataO$alpha_s, dataO$beta_s, dataD$alpha_s, dataD$beta_s,nc) - prob,
			    lower = 0,
			    upper = 1,
			    tol = 1e-6
			    )$root
			  },
			  error = function(e) {
    stop(sprintf("uniroot failed for p=%.4f: %s", prob, e$message))
    NA_real_
   }
  )	}) * (maxResponse - minResponse)
}

rT.BB = function(n, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, overSample=5) {
	nO = n*overSample
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, n=nO, minResponse=minResponse, maxResponse=maxResponse)
	dataO = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check((maxResponse-minResponse)/2, mu_D, sigma_D, minResponse, maxResponse)
	alpha_s = dataO$alpha_s
	beta_s = dataO$beta_s
	alpha_d = dataD$alpha_s
	beta_d = dataD$beta_s

	observed = numeric(0)
cnt = 0
new_n <- ceiling(n * overSample)
while (length(observed) < n) {
 t_start = rbeta(new_n, alpha_s, beta_s)
 raw_duration = rbeta(new_n, alpha_d, beta_d)
 duration = raw_duration * (1 - t_start)
 t_end = t_start + duration
 t_sampled = runif(new_n) #works because beta is in [0,1], the default min and max
 valid = (t_sampled >= t_start) & (t_sampled <= t_end)
 new_observed = t_sampled[valid]
 observed = c(observed, new_observed)
 cnt = cnt + 1
 if (cnt >= 5) {
	 cat(paste0("Number simulated: ", length(observed), "\n"))
  stop("Could not simulate enough observed times. Adjust parameters.")
 }
}
	observed = sample(observed, size = n, replace = FALSE)
	return(minResponse + observed * (maxResponse - minResponse))
}

rCkN.BB = function(n, N, mu_O, sigma_O, mu_D, sigma_D, minResponse = 0, maxResponse = 365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, n=n, N=N, minResponse=minResponse, maxResponse=maxResponse)
	dataO = convert_beta_check((maxResponse - minResponse) / 2, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check((maxResponse - minResponse) / 2, mu_D, sigma_D, minResponse, maxResponse)

	alpha_s = dataO$alpha_s
	beta_s  = dataO$beta_s
	alpha_d = dataD$alpha_s
	beta_d  = dataD$beta_s

	# Total samples needed
	total = n * N

	t_start = rbeta(total, alpha_s, beta_s)
	raw_duration = rbeta(total, alpha_d, beta_d)

	t_start_mat = matrix(t_start, nrow = n)
	raw_duration_mat = matrix(raw_duration, nrow = n)

	duration = raw_duration_mat * (1 - t_start_mat)
	t_end = t_start_mat + duration

	# Find the maximum t_end in each row (each group of N)
	let = apply(t_end, 1, max)

	return(minResponse + let * (maxResponse - minResponse))
}

dCkN.BB = function(x, N, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, N=N, minResponse=minResponse, maxResponse=maxResponse)
	dataO = convert_beta_check(x, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check(x, mu_D, sigma_D, minResponse, maxResponse)
	alpha_s = dataO$alpha_s
	beta_s = dataO$beta_s
	alpha_d = dataD$alpha_s
	beta_d = dataD$beta_s
	return(dCkN.BB.s(dataO$scaled_x,N,alpha_s, beta_s, alpha_d, beta_d)/(maxResponse-minResponse))
}

pCkN.BB = function(q, N, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, N=N, minResponse=minResponse, maxResponse=maxResponse)
	Fx = pC.BB(q, mu_O, sigma_O, mu_D, sigma_D, minResponse, maxResponse)
	return(Fx^N)
}

qCkN.BB = function(p, N, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, N=N, minResponse=minResponse, maxResponse=maxResponse)
	dataO = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check((maxResponse-minResponse)/2, mu_D, sigma_D, minResponse, maxResponse)
	alpha_s = dataO$alpha_s
	beta_s = dataO$beta_s
	alpha_d = dataD$alpha_s
	beta_d = dataD$beta_s
	qfunc = make_quantile_function_from_pdf(ddist= function(x){ dCkN.BB.s(x,N,alpha_s,beta_s,alpha_d,beta_d) }, support = c(0, 1), n = 100)
	return(minResponse + qfunc(p) * (maxResponse - minResponse))
}

#' @importFrom stats integrate
dR.BB = function(x, N, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, N=N, minResponse=minResponse, maxResponse=maxResponse)
	dataO = convert_beta_check(x, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check(x, mu_D, sigma_D, minResponse, maxResponse)
	alpha_s = dataO$alpha_s
	beta_s = dataO$beta_s
	alpha_d = dataD$alpha_s
	beta_d = dataD$beta_s

	densities = sapply(dataO$scaled_x, function(r) {
				  integrate(function(y) dCkN.BB.s(r + y, N, alpha_s, beta_s, alpha_d, beta_d) * dOk1.BB.s(y,N,alpha_s, beta_s), lower = 0, upper = 1-r)$value
	})

	return(densities / (maxResponse-minResponse))
}

rR.BB = function(n, N, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, N=N, n=n, minResponse=minResponse, maxResponse=maxResponse)
	dataO = convert_beta_check(0, mu_O, sigma_O, minResponse, maxResponse) #scaled data not needed
	dataD = convert_beta_check(0, mu_D, sigma_D, minResponse, maxResponse) #scaled data not needed
	alpha_s = dataO$alpha_s
	beta_s = dataO$beta_s
	alpha_d = dataD$alpha_s
	beta_d = dataD$beta_s

	# Total samples needed
	total = n * N

	t_start = rbeta(total, alpha_s, beta_s)
	raw_duration = rbeta(total, alpha_d, beta_d)

	t_start_mat = matrix(t_start, nrow = n)
	raw_duration_mat = matrix(raw_duration, nrow = n)

	duration = raw_duration_mat * (1 - t_start_mat)
	t_end = t_start_mat + duration

	let = apply(t_end, 1, max)
	fst = apply(t_start_mat, 1, min)

	return((let - fst) * (maxResponse - minResponse))
}

#' @importFrom stats dbinom
dPNt.BB = function(x, t, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, res = 1000) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse)
	if(t <= minResponse || t >= maxResponse ) {
		stop("t must be between the minimum and maximum response values.")
	}
	if(res<100 || res>1000000) {
		stop("Please set res to a number between 100 and 1000000 inclusive. Smaller values are faster but less precise.")
	}

	if (any(x < 0 | x > 1)) {
 stop("Provide proportion of individuals in phenophase for input in dPNt.BB. All elements of x must be between 0 and 1 (exclusive).")
}

	N = res
	n = round(x * res)

	Pt = Pt.BB(t, mu_O, sigma_O, mu_D, sigma_D, minResponse, maxResponse)
	return( dbinom(x=n, size=N, prob=Pt) )
}

rPNt.BB = function(n, t, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, res=1000) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, n=n, minResponse=minResponse, maxResponse=maxResponse)
	if(t <= minResponse || t >= maxResponse ) {
		stop("t must be between the minimum and maximum response values.")
	}
	dataO = convert_beta_check(0, mu_O, sigma_O, minResponse, maxResponse) #scaled data not used
	dataD = convert_beta_check(0, mu_D, sigma_D, minResponse, maxResponse) #scaled data not used
	alpha_s = dataO$alpha_s
	beta_s = dataO$beta_s
	alpha_d = dataD$alpha_s
	beta_d = dataD$beta_s

	t = (t-minResponse)/(maxResponse-minResponse)

	total = n * res

	t_start = rbeta(total, alpha_s, beta_s)
	raw_duration = rbeta(total, alpha_d, beta_d)

	t_start_mat = matrix(t_start, nrow = n)
	raw_duration_mat = matrix(raw_duration, nrow = n)

	duration = raw_duration_mat * (1 - t_start_mat)
	t_end_mat = t_start_mat + duration

	hit = (t > t_start_mat) & (t < t_end_mat)

	return(apply(hit, 1, mean)*100)
}

convert_beta_check = function(x, mean, sigma, minResponse, maxResponse) {
	if(minResponse != 0) {
		stop("The minimum must be set to 0. Other values are not currently handled.")
	}
	if(maxResponse <= minResponse) {
		stop("The maximum needs to be greater than the minimum.")
	}
	if( mean < minResponse || mean > maxResponse) {
		stop(paste("The mean needs to be between ", minResponse, " and ", maxResponse))
	}
	if(sigma <=0 ) {
		stop(paste("The standard deviation, sigma, needs to be positive."))
	}

	x_scaled = (x - minResponse) / (maxResponse - minResponse)
	if(any(x_scaled < 0 | x_scaled > 1)) {
		stop("The input values need to be between min and max.")
	}

	sigma = sigma / (maxResponse - minResponse)
	if(sigma > 0.5) {
		stop(paste("The provided sigma is too large for the beta distribtuion. The maximum for your data would be ", (maxResponse-minResponse)*0.5, ". More typical values are less than ", (maxResponse-minResponse)*sqrt(1/12)))
	}


	mean = (mean-minResponse) / (maxResponse - minResponse)

	alpha_s = beta_alpha(mean = mean, sd=sigma)
	beta_s = beta_beta(mean = mean, sd=sigma)

	result = list(
		   scaled_x = x_scaled,
		   scaled_mean = mean,
		   scaled_sigma = sigma,
		   alpha_s = alpha_s,
		   beta_s = beta_s
	)
	return(result)
}

dD.BB.s = Vectorize(function(x, alpha_s, beta_s, alpha_d, beta_d) {
 integrand = function(o) {
  f_o = dbeta(o, alpha_s, beta_s)
  scaled_beta = (1 / (1 - o)) * dbeta(x / (1 - o), alpha_d, beta_d)
  return(f_o * scaled_beta)
 }
 return(integrate(integrand, lower = 1e-10, upper = 1 - x - 1e-10 , rel.tol=1e-6)$value)
})

pD.BB.s = Vectorize(function(q, alpha_s, beta_s, alpha_d, beta_d) {
 integrand = function(d) {
	dD.BB.s(d, alpha_s, beta_s, alpha_d, beta_d)
 }
 return(integrate(integrand, lower = 0, upper = q, rel.tol = 1e-6)$value)
})

dO.BB.s = function(x, alpha_s, beta_s) {
	return(dbeta_safe(x, alpha_s, beta_s) )
}

pO.BB.s = function(q, alpha_s, beta_s) {
	return(pbeta_safe(q, alpha_s, beta_s))
}

dOk1.BB.s = function(x, N, alpha_s, beta_s) {
	return(N * (1 - pO.BB.s(x, alpha_s, beta_s))^(N - 1) * dO.BB.s(x, alpha_s, beta_s))
}

#' @importFrom stats dbeta integrate 
dC.BB.s = function(x, alpha_s, beta_s, alpha_d, beta_d) {
	sapply(x, function(t) {
		    integrand = function(s) {
			    u = (t - s) / (1 - s)
			    f_ts = dbeta(s, alpha_s, beta_s)
			    f_u  = dbeta(u, alpha_d, beta_d) 
			    #f_u  = dbeta(u, alpha_d, beta_d) * u 
			    return(f_ts * f_u / (1 - s))
		    }
		    out = tryCatch(
				    integrate(integrand, lower = 0, upper = t, rel.tol = 1e-6)$value,
				    error = function(e) 1e-12
		    )
		    return(out)
	})
}

pC.BB.s = Vectorize(function(q, alpha_s, beta_s, alpha_d, beta_d) {
			  return(integrate(function(u) dC.BB.s(u, alpha_s, beta_s, alpha_d, beta_d), lower = 0, upper = q, rel.tol = 1e-6)$value)
	})

#' @importFrom stats dbeta pbeta integrate
Pt.BB.s = function(x, alpha_s, beta_s, alpha_d, beta_d) {
	sapply(x, function(t) {
		    integrand = function(s) {
			    f_ts = dbeta(s, alpha_s, beta_s)
			    scaled_tail = (t - s) / (1 - s)
			    tail_prob = 1 - pbeta(scaled_tail, alpha_d, beta_d)
			    f_ts * tail_prob
		    }
		    tryCatch(
				integrate(integrand, lower = 0, upper = t, rel.tol = 1e-6)$value,
				error = function(e) 1e-12
		    )
	})
}

Pt.nc.BB = function(mu_O, mu_D, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, mu_D=mu_D, minResponse=minResponse, maxResponse=maxResponse)
	mu_O = (mu_O-minResponse) / (maxResponse-minResponse)
	mu_D = (mu_D-minResponse) / (maxResponse-minResponse)
	nc = mu_D * (1 - mu_O) + 1e-12
	return(nc)
}

#' @importFrom stats integrate
pT.BB.s = function(q, alpha_s, beta_s, alpha_d, beta_d, nc) {
	probs = sapply(q, function(t) {
			    tryCatch(
					integrate(
						 function(u) Pt.BB.s(u, alpha_s, beta_s, alpha_d, beta_d),
						 lower = 0,
						 upper = t,
						 rel.tol = 1e-6
						 )$value,
					error = function(e) 1e-12
			    )
	})
	return(probs/nc)
}

dCkN.BB.s = function(x,N,alpha_s,beta_s,alpha_d,beta_d) {
	return(N * (pC.BB.s(x,alpha_s,beta_s,alpha_d,beta_d))^(N - 1) * dC.BB.s(x,alpha_s,beta_s,alpha_d,beta_d))
}

