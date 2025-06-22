getProportionOverlap.OC.BB = function(mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	dataO = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check((maxResponse-minResponse)/2, mu_D, sigma_D, minResponse, maxResponse)
	alpha_s = dataO$alpha_s
	beta_s = dataO$beta_s
	alpha_d = dataD$alpha_s
	beta_d = dataD$beta_s
	return(getProportionOverlap.OC.BB.s(alpha_s, beta_s, alpha_d, beta_d))
}

getPeak.T.BB <- function(mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	dataO = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check((maxResponse-minResponse)/2, mu_D, sigma_D, minResponse, maxResponse)
	alpha_s = dataO$alpha_s
	beta_s = dataO$beta_s
	alpha_d = dataD$alpha_s
	beta_d = dataD$beta_s

	obj <- function(t) -Pt.BB.s(t, alpha_s, beta_s, alpha_d, beta_d)
	peak = optimize(obj, interval = c(1e-4, 1 - 1e-4))$minimum
	return(minResponse + peak * (maxResponse-minResponse))
}

E.Ok1.BB = function(N, mu_O, sigma_O, minResponse=0, maxResponse=365, threshApprox=NA) {
	dataO = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	alpha_s = dataO$alpha_s
	beta_s = dataO$beta_s
	vals = minResponse + E.Ok1.BB.s(N, alpha_s, beta_s) * (maxResponse - minResponse)

	if(!is.na(threshApprox) && threshApprox>0) {
		vals[is.na(vals)] = Inf
		vals[is.nan(vals)] = Inf
		vals.approx = E.Kth.approx(N=N, mu=mu_O, sigma=sigma_O, k=1)
		vals[abs(vals - vals.approx) > threshApprox] = vals.approx[abs(vals - vals.approx) > threshApprox]
		nRep = sum(abs(vals - vals.approx) > threshApprox)
		if(nRep>0) {
			warning(paste(nRep, " instances failed numerical integration and were replaced by an asymptotic approximation."))
		}
	}
	return(vals)
}

E.O.BB = function(mu_O, sigma_O, minResponse=0, maxResponse=365) {
	dataO = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	alpha_s = dataO$alpha_s
	beta_s = dataO$beta_s
	mu_O = beta_mean(alpha_s, beta_s)
	return(minResponse + mu_O * (maxResponse - minResponse))
}

E.T.BB = function(mu_O, mu_D, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=NA, mu_D=mu_D, sigma_D=NA, mu_C=NA, sigma_C=NA, N=NA, n=NA, minResponse=minResponse, maxResponse=maxResponse)
	muo = (mu_O-minResponse) / (maxResponse - minResponse)
	mud = (mu_D) / (maxResponse - minResponse)
	e = muo + mud * (1-muo)/2 #this should be checked! // is approximate
	return(minResponse + e * (maxResponse-minResponse)) 
}

E.C.BB = function(mu_O, mu_D, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=NA, mu_D=mu_D, sigma_D=NA, mu_C=NA, sigma_C=NA, N=NA, n=NA, minResponse=minResponse, maxResponse=maxResponse)
	muo = (mu_O-minResponse) / (maxResponse - minResponse)
	mud = (mu_D) / (maxResponse - minResponse)
	e = muo + mud * (1-muo) #this should be checked! // is approximate
	return(minResponse + e * (maxResponse - minResponse))
}

E.CkN.BB = function(N, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365,threshApprox=NA) {
	dataO = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check((maxResponse-minResponse)/2, mu_D, sigma_D, minResponse, maxResponse)
	alpha_s = dataO$alpha_s
	beta_s = dataO$beta_s
	alpha_d = dataD$alpha_s
	beta_d = dataD$beta_s
	vals = minResponse + E.CkN.BB.s(N, alpha_s, beta_s, alpha_d, beta_d) * (maxResponse - minResponse)

	if(!is.na(threshApprox) && threshApprox>0) {
		vals[is.na(vals)] = Inf
		vals[is.nan(vals)] = Inf
		mu_O_scaled = ( mu_O - minResponse ) / ( maxResponse - minResponse )
		mu_D_scaled = ( 1 - mu_O_scaled) * ( mu_D ) / ( maxResponse - minResponse )
		mu_C = minResponse + (mu_O_scaled + mu_D_scaled) * (maxResponse-minResponse)
		sigma_C = sd(rC.BB(n=1000, mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, minResponse=minResponse, maxResponse=maxResponse))
		vals.approx = E.Kth.approx(N=N, mu=mu_C, sigma=sigma_C, k=N)
		vals[abs(vals - vals.approx) > threshApprox] = vals.approx[abs(vals - vals.approx) > threshApprox]
		nRep = sum(abs(vals - vals.approx) > threshApprox)
		if(nRep>0) {
			warning(paste(nRep, " instances failed numerical integration and were replaced by an asymptotic approximation."))
		}
	}
	return(vals)
}

E.D.BB = function(mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, mu_C=NA, sigma_C=NA, N=NA, n=NA, minResponse=minResponse, maxResponse=maxResponse)
	muo = (mu_O-minResponse) / (maxResponse - minResponse)
	mud = (mu_D) / (maxResponse - minResponse)
	e = mud * (1-muo) #this should be checked!
	return(e * (maxResponse - minResponse))
}

SD.Ok1.BB = function(N, mu_O, sigma_O, minResponse=0, maxResponse=365) {
	dataO = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	alpha_s = dataO$alpha_s
	beta_s = dataO$beta_s
	E = E.Ok1.BB.s(N, alpha_s, beta_s)
	integrand <- function(x) { x * x * dOk1.BB.s(x, N, alpha_s, beta_s) } 
	E2 = integrate(integrand, 1e-4, 1-1e-4)$value 
	var = (E2 - E^2) * (maxResponse-minResponse)^2
	if(var < 0) { stop("Error calculating variance for first onset times") }
	return(sqrt(var))
}

SD.O.BB = function(sigma_O) {
	if(sigma_O<=0) {
		stop("sigma_O needs to be positive.")
	}
	return(sigma_O)
}

SD.T.BB = function(mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	dataO = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check((maxResponse-minResponse)/2, mu_D, sigma_D, minResponse, maxResponse)
	alpha_s = dataO$alpha_s
	beta_s = dataO$beta_s
	alpha_d = dataD$alpha_s
	beta_d = dataD$beta_s
	E = (E.T.BB(mu_O, mu_D, minResponse, maxResponse) - minResponse) / (maxResponse - minResponse)
	nc = Pt.nc.BB(mu_O, mu_D, minResponse, maxResponse)
	integrand <- function(x) { x * x * Pt.BB.s(x, alpha_s, beta_s, alpha_d, beta_d)/nc } 
	E2 = integrate(integrand, 0, 1)$value 
	var = (E2 - E^2) * (maxResponse-minResponse)^2
	if(var < 0) { stop("Error calculating variance for observed times") }
	return(sqrt(var))
}

SD.C.BB = function(mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	dataO = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check((maxResponse-minResponse)/2, mu_D, sigma_D, minResponse, maxResponse)
	alpha_s = dataO$alpha_s
	beta_s = dataO$beta_s
	alpha_d = dataD$alpha_s
	beta_d = dataD$beta_s

	E = ( E.C.BB(mu_O, mu_D, minResponse, maxResponse) - minResponse ) / ( maxResponse - minResponse )
	integrand <- function(x) {
		x*x * dC.BB.s(x, alpha_s=alpha_s, beta_s=beta_s, alpha_d=alpha_d, beta_d=beta_d)
	}
	E2  = integrate(integrand, 1e-4, 1-1e-4)$value
	var = (E2 - E^2) * (maxResponse-minResponse)^2
	if(var < 0) { stop("Error calculating variance for observed times") }
	return(sqrt(var))
}

SD.CkN.BB = function(N, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	dataO = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check((maxResponse-minResponse)/2, mu_D, sigma_D, minResponse, maxResponse)
	alpha_s = dataO$alpha_s
	beta_s = dataO$beta_s
	alpha_d = dataD$alpha_s
	beta_d = dataD$beta_s

	E = E.CkN.BB.s(N, alpha_s, beta_s, alpha_d, beta_d)
	integrand <- function(x) { x * x * dCkN.BB.s(x, N, alpha_s, beta_s, alpha_d, beta_d) } 
	E2 = integrate(integrand, 1e-4, 1-1e-4)$value 
	var = (E2 - E^2) * (maxResponse-minResponse)^2
	if(var < 0) { stop("Error calculating variance for first onset times") }
	return(sqrt(var))
}

SD.D.BB = function(mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365) {
	dataO = convert_beta_check((maxResponse-minResponse)/2, mu_O, sigma_O, minResponse, maxResponse)
	dataD = convert_beta_check((maxResponse-minResponse)/2, mu_D, sigma_D, minResponse, maxResponse)
	alpha_s = dataO$alpha_s
	beta_s = dataO$beta_s
	alpha_d = dataD$alpha_s
	beta_d = dataD$beta_s

	mean = E.D.BB(mu_O, sigma_O, mu_D, sigma_D, minResponse, maxResponse)  / (maxResponse - minResponse)

	second_moment_integrand <- function(d) d^2 * dD.BB.s(d, alpha_s, beta_s, alpha_d, beta_d)
	second_moment <- integrate(second_moment_integrand, lower = 0, upper = 1, rel.tol = 1e-6)$value

	var = (second_moment - mean^2 ) * ( maxResponse - minResponse )^2
	sd = sqrt(var)
	return(sd)
}

PI.Ok1.BB = function(N, mu_O, sigma_O, minResponse=0, maxResponse=365, alpha=0.05) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=NA, sigma_D=NA, mu_C=NA, sigma_C=NA, N=N, n=NA, minResponse=minResponse, maxResponse=maxResponse)
	lower = alpha/2
	upper = 1 - alpha/2
	return(qOk1.BB(c(lower,upper), N, mu_O, sigma_O))
}

PI.O.BB = function(mu_O, sigma_O, minResponse=0, maxResponse=365, alpha=0.05) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=NA, sigma_D=NA, mu_C=NA, sigma_C=NA, N=NA, n=NA, minResponse=minResponse, maxResponse=maxResponse)
	if(alpha<=0 || alpha>=1) {
		stop("alpha must be between 0 and 1 exclusive.")
	}
	lower = alpha/2
	upper = 1 - alpha/2
	return(qO.BB(c(lower,upper), mu_O, sigma_O))
}

PI.T.BB = function(mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, alpha=0.05) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, mu_C=NA, sigma_C=NA, N=NA, n=NA, minResponse=minResponse, maxResponse=maxResponse)
	if(alpha<=0 || alpha>=1) {
		stop("alpha must be between 0 and 1 exclusive.")
	}
	lower = alpha/2
	upper = 1 - alpha/2
	return(qT.BB(c(lower,upper), mu_O, sigma_O, mu_D, sigma_D))
}

PI.C.BB = function(mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, alpha=0.05) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, mu_C=NA, sigma_C=NA, N=NA, n=NA, minResponse=minResponse, maxResponse=maxResponse)
	if(alpha<=0 || alpha>=1) {
		stop("alpha must be between 0 and 1 exclusive.")
	}
	lower = alpha/2
	upper = 1 - alpha/2
	return(qC.BB(c(lower,upper), mu_O, sigma_O, mu_D, sigma_D))
}

PI.CkN.BB = function(N, mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, alpha=0.05) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, mu_C=NA, sigma_C=NA, N=N, n=NA, minResponse=minResponse, maxResponse=maxResponse)
	if(alpha<=0 || alpha>=1) {
		stop("alpha must be between 0 and 1 exclusive.")
	}
	lower = alpha/2
	upper = 1 - alpha/2
	return(qCkN.BB(c(lower,upper), N, mu_O, sigma_O, mu_D, sigma_D))
}

PI.D.BB = function(mu_O, sigma_O, mu_D, sigma_D, minResponse=0, maxResponse=365, alpha=0.05) {
	parameter_checks(mu_O=mu_O, sigma_O=sigma_O, mu_D=mu_D, sigma_D=sigma_D, mu_C=NA, sigma_C=NA, N=NA, n=NA, minResponse=minResponse, maxResponse=maxResponse)
	if(alpha<=0 || alpha>=1) {
		stop("alpha must be between 0 and 1 exclusive.")
	}
	lower = alpha/2
	upper = 1 - alpha/2
	return(qD.BB(c(lower,upper), mu_O, sigma_O, mu_D, sigma_D))
}


getMAP.T.BB = function(fileOrData, minResponse=0, maxResponse=365,minS=1, maxS=3000,  init_params = c(180,20,60,7), hyperparameters = c(100, 7, 60, 6, 24, 12, 24, 12), dataProvided=F) {

	if(length(hyperparameters) != 8) {
		print("Provide hyperparameters. Hyperparameters are the mean and sd for: mean onset, mean duration, sd onset, and sd duration. So, there should be 8 hyperparameters total. Defaults are provided.")
		return()
	}

	#convert raw initial parameter values to beta parameters
	init_params[1] = (init_params[1] - minResponse)
	init_params = init_params / (maxResponse - minResponse)
	init_params_beta = init_params
	init_params_beta[1] =  beta_alpha(init_params[1], init_params[2])
	init_params_beta[2] =  beta_beta(init_params[1], init_params[2])
	init_params_beta[3] =  beta_alpha(init_params[3], init_params[4])
	init_params_beta[4] =  beta_beta(init_params[3], init_params[4])

	if(!isBetaFeasible(init_params_beta[1],init_params_beta[2], minS, maxS) || !isBetaFeasible(init_params_beta[3],init_params_beta[4], minS, maxS)) {
		print(init_params)
		print(init_params_beta)
		result = list(error = T, error_m = "Infeasible beta parameters")
		return(result)
	}

	#gets and scales times
	if(dataProvided) {
		observed = (fileOrData - minResponse ) / (maxResponse - minResponse)
	}
	else {
		observed = getObservations(fileOrData,minResponse=minResponse,maxResponse=maxResponse)
	}

	#generate hyperparameters
	hMean_mO = (hyperparameters[1] - minResponse) / (maxResponse - minResponse)
	hSigma_mO = hyperparameters[2] / (maxResponse - minResponse)
	hMean_SDO = hyperparameters[5]  / (maxResponse - minResponse)
	hSigma_SDO = hyperparameters[6] / (maxResponse - minResponse)

	hMean_mD = hyperparameters[3] / (maxResponse - minResponse)
	hSigma_mD = hyperparameters[4] / (maxResponse - minResponse)
	hMean_SDD = hyperparameters[7] / (maxResponse - minResponse)
	hSigma_SDD = hyperparameters[8] / (maxResponse - minResponse)

	a_o_m = beta_alpha(hMean_mO, hSigma_mO)
	b_o_m = beta_beta(hMean_mO, hSigma_mO)
	a_d_m = beta_alpha(hMean_mD, hSigma_mD)
	b_d_m = beta_beta(hMean_mD, hSigma_mD)

	a_o_sd = beta_alpha(hMean_SDO, hSigma_SDO)
	b_o_sd = beta_beta(hMean_SDO, hSigma_SDO)
	a_d_sd = beta_alpha(hMean_SDD, hSigma_SDD)
	b_d_sd = beta_beta(hMean_SDD, hSigma_SDD)

	#check hyperparameter feasibility
	if(!isBetaFeasible(a_o_m, b_o_m,0.01,maxS) || !isBetaFeasible(a_d_m, b_d_m,0.01,maxS) || !isBetaFeasible(a_o_sd, b_o_sd,0.01,maxS) || !isBetaFeasible(a_d_sd, b_d_sd,0.01,maxS)) {
		print(paste("hm_mo: ", hMean_mO, " hsd_mo: ", hSigma_mO, " hm_sdo: ", hMean_SDO, " hsd_sdo: ", hSigma_SDO))
		print(c(a_o_m,b_o_m,a_o_sd,b_o_sd))
		print(paste("hm_md: ", hMean_mD, " hsd_md: ", hSigma_mD, " hm_sdd: ", hMean_SDD, " hsd_sdd: ", hSigma_SDD))
		print(c(a_d_m,b_d_m,a_d_sd,b_d_sd))
		result = list(error = T, error_m = "Infeasible beta hyperparameters")
		return(result)
	}

	cat("Starting MAP optimization. This may take a while.\n")

	result <- optim(
					par = init_params_beta,
					fn = neg_log_posterior,
					t_obs = observed,
					a_o_m = a_o_m, b_o_m = b_o_m, a_d_m = a_d_m, b_d_m = b_d_m,	#beta priors for the means
					a_o_sd = a_o_sd, b_o_sd = b_o_sd, a_d_sd = a_d_sd, b_d_sd = b_d_sd,	#beta priors for the SDs
					method = "L-BFGS-B",
					lower = c(minS,minS,minS,minS),
					upper = c(maxS,maxS,maxS,maxS)
	)

	cat("Finished MAP optimization.\n")

	par = result$par
	mO = minResponse + beta_mean(par[1],par[2]) * (maxResponse - minResponse)
	sdO = beta_sd(par[1],par[2]) * (maxResponse - minResponse)
	mD = beta_mean(par[3],par[4]) * (maxResponse - minResponse)
	sdD = beta_sd(par[3],par[4]) * (maxResponse - minResponse)
	result$par_orig_scale = c(mO, sdO, mD, sdD)
	result$error = F

	if(result$convergence != 0) {
		result$error = T
		result$error_m = "optimization failed to converged."
	}

	return(result)
}

getMLE.T.BB = function(fileOrData, minResponse=0, maxResponse=365, minS=1, maxS=3000, init_params = c(180, 20, 60, 7), dataProvided=F) {

	#convert raw initial parameter values to beta parameters
	init_params[1] = init_params[1] - minResponse
	init_params = (init_params) / (maxResponse - minResponse)
	init_params_beta = init_params;
	init_params_beta[1] =  beta_alpha(init_params[1], init_params[2])
	init_params_beta[2] =  beta_beta(init_params[1], init_params[2])
	init_params_beta[3] =  beta_alpha(init_params[3], init_params[4])
	init_params_beta[4] =  beta_beta(init_params[3], init_params[4])

	if(!isBetaFeasible(init_params_beta[1],init_params_beta[2], minS, maxS) || !isBetaFeasible(init_params_beta[3],init_params_beta[4], minS, maxS))
	{
		result = list(error = T, error_m = "Infeasible beta parameters")
		return(result)
	}

	#gets and scales times
	if(dataProvided) {
		observed = (fileOrData - minResponse ) / (maxResponse - minResponse)
	}
	else {
		observed = getObservations(fileOrData,minResponse=minResponse,maxResponse=maxResponse)
	}


	cat("Starting MLE optimization. This may take a while.\n")
	result <- optim(
					par = init_params_beta,
					fn = neg_loglik_observed,
					data = observed,
					method = "L-BFGS-B",
					lower = c(minS,minS,minS,minS),
					upper = c(maxS,maxS,maxS,maxS)
	)
	par = result$par
	mO = minResponse + beta_mean(par[1],par[2]) * (maxResponse - minResponse)
	sdO = beta_sd(par[1],par[2]) * (maxResponse - minResponse)
	mD = beta_mean(par[3],par[4]) * (maxResponse - minResponse)
	sdD = beta_sd(par[3],par[4]) * (maxResponse - minResponse)
	result$par_orig_scale = c(mO, sdO, mD, sdD)
	result$error = F

	if(result$convergence != 0) {
		result$error = T
	}
	cat("Finished MLE optimization.\n")

	return(result)
}


neg_loglik_observed <- function(param, data) {
	return(-loglik_observed(param,data) )
}

loglik_observed <- function(param, data) {
	a_s <- param[1]; b_s <- param[2]
	a_d <- param[3]; b_d <- param[4]
	if (any(param <= 0)) return(1e10)
	mu_s = a_s/(a_s+b_s)
	mu_d = a_d/(a_d+b_d)
	nc = mu_d * (1 - mu_s) + 1e-12

	dens = Pt.BB.s(data, a_s, b_s, a_d, b_d) / nc
	if (any(!is.finite(dens)) || any(dens <= 0)) stop("Integration error in likelihood calculation.")
	sum(log(dens + 1e-12))  # avoid log(0)
}

neg_log_posterior <- function(params, t_obs, a_o_m, b_o_m, a_d_m, b_d_m, a_o_sd, b_o_sd, a_d_sd, b_d_sd) {
	alpha_s <- params[1]
	beta_s <- params[2]
	alpha_d <- params[3]
	beta_d <- params[4]

	mu_o = beta_mean(alpha_s,beta_s)
	sd_o = beta_sd(alpha_s,beta_s)
	mu_d = beta_mean(alpha_d,beta_d)
	sd_d = beta_sd(alpha_d,beta_d)


	if (alpha_s <= 0 || beta_s <= 0 || sd_o <= 0 || sd_o > 0.25 ||
		alpha_d <= 0 || beta_d <= 0 || sd_d <= 0 || sd_d > 0.25) {
		print("Parameters are out of range during inference. Try setting priors")
		return(Inf)
	}

	# Log-likelihood
	log_lik <- loglik_observed(c(alpha_s,beta_s,alpha_d,beta_d),t_obs)

	# Log-prior
	log_prior <- dbeta(mu_o, a_o_m, b_o_m, log = TRUE) +
		dbeta(mu_d, a_d_m, b_d_m, log = TRUE) +
		dbeta(sd_o, a_o_sd, b_o_sd, log = TRUE) +
		dbeta(sd_d, a_d_sd, b_d_sd, log = TRUE)

	# negative log-posterior
	return(- (log_lik + log_prior))  
}

getProportionOverlap.OC.BB.s = function(alpha_s, beta_s, alpha_d, beta_d) {

	f <- function(x) dO.BB.s(x, alpha_s=alpha_s, beta_s=beta_s)
	g <- function(x) dC.BB.s(x, alpha_s=alpha_s, beta_s=beta_s, alpha_d=alpha_d, beta_d=beta_d)

	result <- try(integrate(function(x) pmin(f(x), g(x)), 0.00001, 0.9999), silent = TRUE)
	if (inherits(result, "try-error")) return(NA) else return(result$value)

}

E.CkN.BB.s = Vectorize(function(N, alpha_s, beta_s, alpha_d, beta_d) {
						   integrand <- function(x) { x * dCkN.BB.s(x, N=N, alpha_s=alpha_s, beta_s=beta_s, alpha_d=alpha_d, beta_d=beta_d) }
						   expected = integrate(integrand, 1e-4, 1 - 1e-4)$value
						   return(expected)
		} )

E.Ok1.BB.s = Vectorize(function(N, alpha_s, beta_s) {
						   integrand <- function(x) { x * dOk1.BB.s(x, N=N, alpha_s=alpha_s, beta_s=beta_s) }
						   expected = integrate(integrand, 1e-4, 1 - 1e-4)$value
						   return(expected)
		})

