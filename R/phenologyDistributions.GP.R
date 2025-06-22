dD.GP = function(x, mu_D) {
	parameter_checks(mu_D=mu_D)
	if(x == mu_D) { return(1) }
	return(0)
}

pD.GP = function(q, mu_D) {
	parameter_checks(mu_D=mu_D)
	if(q<mu_D) { return(0) }
	return(1)
}

qD.GP = function(p, mu_D) {
	parameter_checks(mu_D=mu_D)
	return(mu_D)
}

rD.GP = function(n, mu_D) {
	parameter_checks(mu_D=mu_D, n=n)
	return(rep(mu_D,n))
}

dO.GP = function(x, mu_O, sigma) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma)
	return(dnorm(x, mu_O, sigma) )
}

pO.GP = function(q, mu_O, sigma) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma)
	return(pnorm(q, mu_O, sigma) )
}

qO.GP = function(p, mu_O, sigma) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma)
	return(qnorm(p, mu_O, sigma) )
}

rO.GP = function(n, mu_O, sigma) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma, n=n)
	return(rnorm(n, mu_O, sigma) )
}

dOk1.GP = function(x, N, mu_O, sigma) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma, N=N)
	N * (1 - pO.GP(x, mu_O, sigma))^(N - 1) * dO.GP(x, mu_O, sigma)
}

qOk1.GP = function(p, N, mu_O, sigma) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma, N=N)
	q_fst = qnorm(1 - (1-p)^(1/N), mu_O, sigma)
	return(q_fst)
}

pOk1.GP <- function(q, N, mu_O, sigma) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma, N=N)
	Fx <- pnorm(q=q, mean = mu_O, sd = sigma)
	return(1 - (1 - Fx)^N)
}

rOk1.GP = function(n, N, mu_O, sigma) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma, N=N, n=n)
	return(replicate(n, min(rnorm(N, mean = mu_O, sd = sigma))))
}

dC.GP = function(x, mu_C, sigma) {
	parameter_checks(mu_C=mu_C, sigma_C=sigma)
	return(dnorm(x, mu_C, sigma))
}

pC.GP = function(q, mu_C, sigma) {
	parameter_checks(mu_C=mu_C, sigma_C=sigma)
	return(pnorm(q, mu_C, sigma))
}

qC.GP = function(p, mu_C, sigma) {
	parameter_checks(mu_C=mu_C, sigma_C=sigma)
	return(qnorm(p, mu_C, sigma))
}

rC.GP = function(n, mu_C, sigma) {
	parameter_checks(mu_C=mu_C, sigma_C=sigma, n=n)
	return(rnorm(n, mu_C, sigma))
}

dCkN.GP = function(x, N, mu_C, sigma) {
	parameter_checks(mu_C=mu_C, sigma_C=sigma, N=N)
	return(N * (pC.GP(x,mu_C, sigma))^(N - 1) * dC.GP(x,mu_C, sigma))
}

qCkN.GP = function(p, N, mu_C, sigma) {
	parameter_checks(mu_C=mu_C, sigma_C=sigma, N=N)
	return(qnorm(p^(1 / N), mean = mu_C, sd = sigma))
}

pCkN.GP = function(q, N, mu_C, sigma) {
	parameter_checks(mu_C=mu_C, sigma_C=sigma, N=N)
	Fx <- pnorm(q, mean = mu_C, sd = sigma)
	return(Fx^N)
}

rCkN.GP = function(n, N, mu_C, sigma) {
	parameter_checks(mu_C=mu_C, sigma_C=sigma, N=N, n=n)
	return(replicate(n, max(rnorm(N, mean = mu_C, sd = sigma))))
}

Pt.GP = function(x, mu_O, mu_C, sigma, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma, mu_C=mu_C, minResponse=minResponse, maxResponse=maxResponse)
	if(maxResponse<=minResponse) { stop("The maximum response time must be greater than the minimum response time.") }
	if(mu_C<=mu_O) { stop("The mean start time must be before the mean end time.") }
	if(mu_C > maxResponse) { stop("The mean end time must be earlier than the maximum response time.") }
	if(mu_O < minResponse) { stop("The mean start time must be later than the minimum response time.") }

	prob = pnorm(x,mu_O,sigma)*(1-pnorm(x,mu_C,sigma))
	return(prob)
}

Pt.nc.GP = function(mu_O, mu_C, sigma, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma, mu_C=mu_C, minResponse=minResponse, maxResponse=maxResponse)
	return(integrate(function(x) Pt.GP(x, mu_O, mu_C, sigma, minResponse, maxResponse), lower = minResponse, upper = maxResponse, rel.tol = 1e-8)$value)
}

dT.GP = function(x, mu_O, mu_C, sigma, minResponse=0, maxResponse=365, nc=NULL) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma, mu_C=mu_C, minResponse=minResponse, maxResponse=maxResponse)
	numerator = Pt.GP(x, mu_O, mu_C, sigma, minResponse, maxResponse)
	if(is.null(nc)) {
		nc = Pt.nc.GP(mu_O, mu_C, sigma, minResponse, maxResponse)
	}
	return(numerator / nc)
}

pT.GP = Vectorize(function(q, mu_O, mu_C, sigma, minResponse=0, maxResponse=365, nc=NULL) {
					  parameter_checks(mu_O=mu_O, sigma_O=sigma, mu_C=mu_C, minResponse=minResponse, maxResponse=maxResponse)
					  numerator <- integrate(function(x) Pt.GP(x, mu_O, mu_C, sigma, minResponse, maxResponse), lower = minResponse, upper = q, rel.tol = 1e-8)$value
					  if(is.null(nc)) {nc = Pt.nc.GP(mu_O, mu_C, sigma, minResponse, maxResponse)}
					  return(numerator / nc)
})

qT.GP = Vectorize(function(p, mu_O, mu_C, sigma, minResponse=0, maxResponse=365, nc=NULL) {
					  parameter_checks(mu_O=mu_O, sigma_O=sigma, mu_C=mu_C, minResponse=minResponse, maxResponse=maxResponse)
					  if(is.null(nc)) {nc = Pt.nc.GP(mu_O, mu_C, sigma, minResponse, maxResponse)}
					  lower <- mu_O - 5 * sigma  # reasonable lower bound
					  upper <- mu_C + 5 * sigma  # reasonable upper bound
					  root_fun <- function(t) { pT.GP(t, mu_O, mu_C, sigma, minResponse, maxResponse, nc) - p }
					  res <- tryCatch( { uniroot(root_fun, lower = lower, upper = upper, tol = 1e-8)$root },
						  error = function(e) {
							  warning("T quantile function failed. Returning NA.")
							  stop("T quantile function failed. Returning NA.")
							  NA
						  }
					  )
					  return(res)
})

rT.GP = function(n, mu_O, mu_C, sigma, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma, mu_C=mu_C, n=n, minResponse=minResponse, maxResponse=maxResponse)
	d = mu_C - mu_O
	if(d<=0) { stop("Mean onset must be less than mean cessation.") }
	onset = rO.GP(n=n,mu_O=mu_O,sigma=sigma)
	cessation = onset+d
	T = runif(n,onset, cessation)
	return(T)
}

dR.GP = function(x, N, mu_O, mu_C, sigma, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma, mu_C=mu_C, N=N, minResponse=minResponse, maxResponse=maxResponse)
	if (any(x < 0 | x > (maxResponse - minResponse))) {
		stop("Population phenophase range time x must be within [0, maxResponse - minResponse].")
	}
	return(sapply(x, function(x) {
					  integrate(function(y) dCkN.GP(x + y, N, mu_C, sigma) * dOk1.GP(y,N,mu_O, sigma), lower = minResponse, upper = maxResponse-x)$value
					  }))
}

rR.GP = function(n, N, mu_O, mu_C, sigma) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma, mu_C=mu_C, n=n, N=N, minResponse=minResponse, maxResponse=maxResponse)

	d <- mu_C - mu_O
	if(d<=0) { stop("Mean onset must be less than mean cessation.") }
	pop_matrix <- matrix(rO.GP(n * N, mu_O, sigma), nrow = n, ncol = N)
	rSim <- apply(pop_matrix, 1, function(row) (max(row) + d) - min(row))
	return(rSim)
}

dPNt.GP = function(x, t, mu_O, mu_C, sigma, minResponse=0, maxResponse=365, res = 1000) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma, mu_C=mu_C, minResponse=minResponse, maxResponse=maxResponse)
	if(res<100 || res>1000000) {
		stop("Please set res to a number between 100 and 1000000 inclusive. Smaller values are faster but less precise. The default gives precision to the thousandth's place.")
	}

	if (any(x < 0 | x > 1)) {
		stop("Provide proportions of individuals in the phenophase for x in dPNt.GP. All x values must be between 0 and 1 since they represent proportions.")
	}
	if(t<minResponse || t>maxResponse) {
		return(rep(0, length(x)))
	}

	N = res
	n = round(x * res)

	Pt = Pt.GP(t, mu_O, mu_C, sigma, minResponse, maxResponse)
	return( dbinom(x=n, size=N, prob=Pt))

}

rPNt.GP = function(n, t, mu_O, mu_C, sigma, minResponse=0, maxResponse=365, res=1000, asPercent=TRUE) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma, mu_C=mu_C, n=n, minResponse=minResponse, maxResponse=maxResponse)
	    if(t <= 0) {
        stop("t must be greater than 0.")
    }
	s = 1
	if(asPercent) {
		s = 100
	}
	total <- n * res

	t_start <- rnorm(total, mu_O, sigma)
	t_end  = rnorm(total, mu_C, sigma)
	t_start_mat <- matrix(t_start, nrow = n)
	t_end_mat <- matrix(t_end, nrow = n)

	hit = (t > t_start_mat) & (t < t_end_mat)

	return(apply(hit, 1, mean)*s)
}
