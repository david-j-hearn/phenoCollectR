extract_draw  =  function(posterior_samples, draw_index, covariate_names, betaSuffix="_O", alphaSuffix="_O") {
	draw  =  posterior_samples[draw_index, ]

	# Extract alpha (intercept)
	alpha  =  as.numeric(draw[[paste0("alpha", alphaSuffix)]])

	# Extract beta vector (assumes beta[1], beta[2], ..., beta[K])
	beta_indices  =  paste0("beta" , betaSuffix, "[", seq_along(covariate_names), "]")
	beta  =  as.numeric(draw[, beta_indices])
	names(beta)  =  covariate_names

	#extract sigma
	sigma  =  as.numeric(draw[["sigma"]])

	list(alpha = alpha, beta = beta, sigma = sigma)
}

predict_response  =  function(covariate_df, posterior_draw) {
# return a vector of predictions, one for each provided set of covariate values, based on provided posterior draw
	# covariate_df: data.frame (n_samples x K)
	# posterior_draw: list or named vector with elements:
	#                 - alpha: scalar
	#                 - beta: named vector of coefficients

	# Ensure column names match
	stopifnot(all(names(posterior_draw$beta) %in% colnames(covariate_df)))

	# Subset and order covariates to match beta
	X  =  as.matrix(covariate_df[, names(posterior_draw$beta), drop = FALSE])

	# Compute linear predictor
	linear_pred  =  X %*% posterior_draw$beta + posterior_draw$alpha
	return(as.vector(linear_pred))
}

sample_posterior_predictions  =  function(posterior_samples, cov_samplesO, covariate_namesO, cov_samplesD, covariate_namesD, n_draws = 100, N = 1000) {
	total_draws  =  nrow(posterior_samples)

	# Random sample of draw indices
	sampled_indices  =  sample(total_draws, size = n_draws, replace = FALSE)

	# Predict for each sampled draw
	preds  =  sapply(sampled_indices, function(i) {
						drawO  =  extract_draw(posterior_samples=posterior_samples, draw_index=i, covariate_names=covariate_namesO, betaSuffix="_O", alphaSuffix="_O")
						drawD  =  extract_draw(posterior_samples=posterior_samples, draw_index=i, covariate_names=covariate_namesD, betaSuffix="_D", alphaSuffix="_D")

						# Predict using linear model
						Os = predict_response(cov_samplesO, drawO)
						Ds = predict_response(cov_samplesD, drawD)
						Cs = Os+Ds
						Ts = Os+Ds/2
						Ok1s = E.Kth.approx(N=N, mu=Os, sigma=drawO$sigma, k=1)
						CkNs = E.Kth.approx(N=N, mu=Cs, sigma=drawO$sigma, k=N) #same sigma for cessation and onset
						out = c(mOk1 = mean(Ok1s), mO = mean(Os), mT = mean(Ts), mC = mean(Cs), mCkN = mean(CkNs), S = drawO$sigma)
						return(out)
})

	# Return: matrix with rows = data points, columns = sampled draws
	preds
}


#' @export
#' @importFrom posterior as_draws_df
#' @importFrom copula pobs normalCopula fitCopula
#' @importFrom dplyr tibble bind_rows group_by filter ungroup n row_number group_modify
#' @importFrom mgcv gam
#' @importFrom stats predict
#' @importFrom splines bs
#' @importFrom cowplot plot_grid
makePosteriorPredictivePlot = function(stanResult, responseData, responseVariableName="DOY", targetCovariateName, onsetCovariateData, durationCovariateData, n_samples=100, n_draws=100, n_hist=10000, n_bin=100, resolution=50, slice = 0.25, N=1000, minResponse=0, maxResponse=365, adjustYLim=T, smooth=c("GAM", "LOESS", "SPLINE", "NONE"), keepSigma=T) {

	smooth = match.arg(smooth)

	if (!(targetCovariateName %in% colnames(onsetCovariateData))) {
		stop("The provided target covariate must be present in the onset covariate data.")
	}

	#get draws
	posterior_samples  =  as_draws_df(stanResult$result$draws())

	## DANIEL: Try to prevent warning messages:
	posterior_samples = as_tibble( posterior_samples )

	#onset
	covariate_namesO = colnames(onsetCovariateData)
	covariate_namesD = colnames(durationCovariateData)

	target_idx = which(colnames(onsetCovariateData) == targetCovariateName)

	nColO = ncol(onsetCovariateData)
	nColD = ncol(durationCovariateData)

	if((nColO == 1 || nColD == 1) && nColO != nColD) {
		stop("Not currently implemented: one covariate for either onset or duration, but not both one is not currently implemented.")
	}


	#Set up values to sample covariates using Gaussian copula to model correlation structure among covariates
	if(nColO > 1) {
		cat("Setting up copula for sampling covariates.\n")

		UO  =  pobs(as.matrix(onsetCovariateData))  # Pseudo-observations (uniform marginals)
		copO  =  normalCopula(dim = ncol(UO), dispstr = "un")
		fitO  =  fitCopula(copO, UO, method = "ml")

		UD  =  pobs(as.matrix(durationCovariateData))  # Pseudo-observations (uniform marginals)
		copD  =  normalCopula(dim = ncol(UD), dispstr = "un")
		fitD  =  fitCopula(copD, UD, method = "ml")
	}

	#needed for the time slice histograms
	minTarget = min(onsetCovariateData[[targetCovariateName]])
	maxTarget = max(onsetCovariateData[[targetCovariateName]])
	rangeTarget = maxTarget - minTarget

	#requires that the target be in the onset covariate data set
	target_vals  =  seq(minTarget, maxTarget, length.out = resolution)
	idx_s1 = round(slice * resolution)
	idx_s2 = round((1-slice) * resolution)

	#get positions where time slice plots will be made
	s1 = target_vals[idx_s1]
	s2 = target_vals[idx_s2]

	results  =  data.frame()
	for (x in target_vals) {
		cat(paste0("Processing ", x, "\n"))

		#get simulated covariate data
		if(nColO > 1) {
			cov_samplesO  =  sample_conditional_covariates(
														  x_target = x,
														  x_column = target_idx,
														  covars = onsetCovariateData,
														  copula_fit = fitO,
														  n_samples = n_samples
			)

			cov_samplesD  =  sample_conditional_covariates(
														  x_target = x,
														  x_column = target_idx,
														  covars = durationCovariateData,
														  copula_fit = fitD,
														  n_samples = n_samples
			)


			pred_means = sample_posterior_predictions(
													  posterior_samples=posterior_samples,
													  cov_samplesO = cov_samplesO,
													  covariate_namesO = covariate_namesO,
													  cov_samplesD = cov_samplesD,
													  covariate_namesD = covariate_namesD,
													  n_draws = n_draws,
													  N=N)

			pred_summary  =  tibble(
								   x_target = x,
								   series = c("Ok1", "O", "T", "C", "CkN", "Sigma"),
								   mean = apply(pred_means, 1, mean),
								   lower = apply(pred_means, 1, function(x) quantile(x, 0.025)),
								   upper = apply(pred_means, 1, function(x) quantile(x, 0.975))
			)
		}
		else {
			cov_samplesO = data.frame(targetCovariateName = c(x)) # a single draw when no other covariates present
			cov_samplesD = data.frame(targetCovariateName = c(x)) # a single draw when no other covariates present
			pred_means = sample_posterior_predictions(
													  posterior_samples=posterior_samples,
													  cov_samplesO = cov_samplesO,
													  covariate_namesO = covariate_namesO,
													  cov_samplesD = cov_samplesD,
													  covariate_namesD = covariate_namesD,
													  n_draws = n_draws,
													  N = N)
			pred_summary  =  tibble(
								   x_target = x,
								   series = c("Ok1", "O", "T", "C", "CkN", "Sigma"),
								   mean = apply(pred_means, 1, mean),
								   lower = apply(pred_means, 1, function(x) quantile(x, 0.025)),
								   upper = apply(pred_means, 1, function(x) quantile(x, 0.975))
			)
		}
		if(x == s1) {
			s1_summary = pred_summary
		}
		if(x == s2) {
			s2_summary = pred_summary
		}

		results  =  bind_rows(results, pred_summary)
	}

	minY = min(pred_summary$lower)
	maxY = max(pred_summary$upper)


	#make the time slice panels
	cat("Making the time slice panels.\n")

	mu_O1 = s1_summary[s1_summary$series == "O",]$mean
	mu_C1 = s1_summary[s1_summary$series == "C",]$mean
	sigma1= s1_summary[s1_summary$series == "Sigma",]$mean

	Ok1_s1 = rOk1(n=n_hist,N=N,mu_O=mu_O1,sigma_O=sigma1,minResponse=minResponse, maxResponse=maxResponse, type="GP")
	CkN_s1 = rCkN(n=n_hist,N=N,mu_O=mu_O1,sigma_O=sigma1, mu_D=mu_C1-mu_O1, minResponse=minResponse, maxResponse=maxResponse, type="GP")

	mu_O2 = s2_summary[s2_summary$series == "O",]$mean
	mu_C2 = s2_summary[s2_summary$series == "C",]$mean
	sigma2= s2_summary[s2_summary$series == "Sigma",]$mean

	Ok1_s2 = rOk1(n=n_hist,N=N,mu_O=mu_O2,sigma_O=sigma2,minResponse=minResponse, maxResponse=maxResponse, type="GP")
	CkN_s2 = rCkN(n=n_hist,N=N,mu_O=mu_O2,sigma_O=sigma2, mu_D=mu_C2-mu_O2, minResponse=minResponse, maxResponse=maxResponse, type="GP")

	xlim = c(min(Ok1_s1,Ok1_s2,CkN_s1,CkN_s2), max(Ok1_s1,Ok1_s2,CkN_s1,CkN_s2))

	pS1 = makeSimulatedAndTheoreticalOverlayGraph(mu_O=mu_O1, mu_C=mu_C1, sigma=sigma1, n_hist=n_hist, N=N, minResponse=minResponse, maxResponse=maxResponse, responseVariableName=NULL, nBin=n_bin, xlim=xlim)
	pS2 = makeSimulatedAndTheoreticalOverlayGraph(mu_O=mu_O2, mu_C=mu_C2, sigma=sigma2, n_hist=n_hist, N=N, minResponse=minResponse, maxResponse=maxResponse, responseVariableName="DOY at Dotted Slice",nBin=n_bin, xlim=xlim)

	#get rid of the sigma data for graphing
	if(!keepSigma) {
	results  =  results[results$series != "Sigma", ]
	}

	#The last data point appears to be highly biased (error in sampling from copula?), so graphing all but the last point
	filtered_results  =  results %>%
		group_by(series) %>%
		filter(row_number() < n()) %>%  # keep all but last row per group
		ungroup()

	if(is.na(smooth)) {
		stop("Missing the smoothing type.")
	}
	else if(smooth == "GAM") {
		if(resolution <=10) {
			stop("Please increase resolution to greater than 10 before attempting to use GAM to smooth the posterior predictive plot.")
		}
		filtered_results  =  filtered_results %>%
			group_by(series) %>%
			group_modify(~ {
									tibble(
									x_target = .x$x_target,
									series = .x$series,
									mean  = predict(gam(mean  ~ s(x_target), data = .x, method = "REML")),
									lower = predict(gam(lower ~ s(x_target), data = .x, method = "REML")),
									upper = predict(gam(upper ~ s(x_target), data = .x, method = "REML"))
							 )
			}) %>%
			ungroup()
	}
	else if(smooth == "LOESS") {
		filtered_results  =  filtered_results %>%
			group_by(series) %>%
			group_modify(~ {
									tibble(
									x_target = .x$x_target,
									series = .x$series,
									mean  = predict(loess(mean  ~ x_target, data = .x, span = 0.3)),
									lower = predict(loess(lower ~ x_target, data = .x, span = 0.3)),
									upper = predict(loess(upper ~ x_target, data = .x, span = 0.3))
							 )
			}) %>%
			ungroup()
	}
	else if(smooth == "SPLINE") {
		filtered_results  =  filtered_results %>%
			group_by(series) %>%
			group_modify(~ {
									tibble(
									x_target = .x$x_target,
									series = .x$series,
									mean  = predict(lm(mean  ~ bs(x_target, df = 5), data = .x)),
									lower = predict(lm(lower ~ bs(x_target, df = 5), data = .x)),
									upper = predict(lm(upper ~ bs(x_target, df = 5), data = .x))
							 )
			}) %>%
			ungroup()
	}


	#set up the observed data for plotting
	observed_data = data.frame(
							   x = onsetCovariateData[[targetCovariateName]],
							   y = responseData
	)
	#create the posterior predictive plot
	pPPP = posterior_predictive_graphic(observed_data, filtered_results,targetCovariateName,responseVariableName, timeSlice1 = s1, timeSlice2 = s2, adjustYLim=adjustYLim)

	slices = plot_grid(pS1, pS2, nrow = 2, ncol = 1, labels = c("B", "C"), label_x = 0, label_y = 0.1)
	figure = plot_grid(pPPP, slices, nrow = 1, ncol = 2, labels = c("A"), label_x = 0, label_y = 0.1, rel_widths = c(2, 1))

	return(figure)
}

#' @importFrom ggplot2 ggplot geom_ribbon aes geom_line geom_point scale_color_manual scale_fill_manual geom_vline theme_minimal labs theme element_rect ylim coord_cartesian
posterior_predictive_graphic = function(observed_data, filtered_results,targetCovariateName,responseVariableName,legend=TRUE, timeSlice1, timeSlice2, adjustYLim = T, minResponse=0, maxResponse=365) {

	if(adjustYLim) {
		ylim= c(min(filtered_results$lower), max(filtered_results$upper))
	}
	pPPP = ggplot() +
		# Posterior predictive ribbons
		geom_ribbon(
					data = filtered_results,
					aes(x = x_target, ymin = lower, ymax = upper, fill = series),
					alpha = 0.2,
					color = NA
					) +
		# Posterior predictive mean lines
		geom_line(
				  data = filtered_results,
				  aes(x = x_target, y = mean, color = series),
				  size = 0.5
				  ) +
		# Observed data points
		geom_point(
				   data = observed_data,
				   aes(x = x, y = y),
				   color = "black",
				   alpha = 0.5,
				   size = 1.5
				   ) +
		# Custom color and fill mapping
		scale_color_manual(
						   values = c(
									  Ok1 = "yellow",
									  O   = "red",
									  T   = "purple",
									  C   = "blue",
									  CkN = "cyan"
						   )
						   ) +
		scale_fill_manual(
						  values = c(
									 Ok1 = "yellow",
									 O   = "red",
									 T   = "purple",
									 C   = "blue",
									 CkN = "cyan"
						  )
						  ) +
		geom_vline(xintercept = c(timeSlice1, timeSlice2), linetype = "dashed", color = "black") +

		theme_minimal() +
		labs(
			 x = targetCovariateName,
			 y = responseVariableName,
			 color = "Series",
			 fill = "Series"
			 ) +
		theme(
			  legend.position = c(0.1, 0.95),  # slightly inset from top-left corner
			  legend.justification = c("left", "top"),
			  #legend.background = element_rect(fill = alpha("white", 0.6), color = NA)  # semi-transparent bg
			  #legend.background = element_rect(fill = "white", color = "black", linewidth=0.5)  # opaque bg, black outline around legend
			  legend.background = element_rect(fill = "white", color=NA )  # opaque bg
		)
		if(!legend) {
			pPPP = pPPP + theme(legend.position = "none")
		}
		if(adjustYLim) {
			print(ylim)
			pPPP = pPPP + ylim(ylim[1], ylim[2])
		}
		else {
			pPPP = pPPP + coord_cartesian(ylim = c(minResponse, maxResponse))
		}
		return(pPPP)
}

#' @importFrom ggplot2 ggplot geom_histogram aes scale_fill_manual scale_color_manual stat_function xlim geom_vline labs theme_minimal geom_area geom_line theme element_rect
#' @importFrom dplyr tibble
#' @importFrom tidyr pivot_longer
makeSimulatedAndTheoreticalOverlayGraph= function(mu_O, mu_C, sigma, n_hist, N, minResponse, maxResponse, responseVariableName="DOY", xlim=c(0,365), legend=FALSE, nBin=50, includeHistograms=FALSE) {
	# Step 1 & 2: Generate samples and put in a tidy df
	mu_D = mu_C-mu_O

	Ok1_s = rOk1(n=n_hist,N=N,mu_O=mu_O,sigma_O=sigma,minResponse=minResponse, maxResponse=maxResponse, type="GP")
	O_s = rO(n=n_hist, mu_O=mu_O, sigma_O=sigma, minResponse=minResponse, maxResponse=maxResponse, type="GP")
	T_s = rT(n=n_hist, mu_O=mu_O, sigma_O=sigma, mu_D=mu_D, minResponse=minResponse, maxResponse=maxResponse, type="GP")
	C_s = rC(n=n_hist, mu_O=mu_O, sigma_O=sigma, mu_D=mu_D, minResponse=minResponse, maxResponse=maxResponse, type="GP")
	CkN_s = rCkN(n=n_hist,N=N,mu_O=mu_O,sigma_O=sigma, mu_D=mu_D, minResponse=minResponse, maxResponse=maxResponse, type="GP")

	mOk1_s = mean(Ok1_s)
	mO_s = mean(O_s)
	mT_s = mean(T_s)
	mC_s = mean(C_s)
	mCkN_s = mean(CkN_s)

	xIs = c(mOk1_s, mO_s, mT_s, mC_s, mCkN_s)

	if(includeHistograms) {
		df  =  data.frame(
						 value = c(Ok1_s, O_s, T_s, C_s, CkN_s),
						 dist = factor(rep(c("Ok1", "O", "T", "C", "CkN"), each = n_hist))
		)

		# Step 3 & 4: Plot histograms and overlay theoretical density
		p = ggplot(df, aes(x = value, fill = dist, color = dist)) +
			geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.3, bins = nBin) +
			scale_fill_manual(
							  values = c(
										 "Ok1" = "yellow",
										 "O"   = "red",
										 "T"   = "purple",
										 "C"   = "blue",
										 "CkN" = "cyan"
							  )
							  ) +
			scale_color_manual(
							   values = c(
										  "Ok1" = "yellow",
										  "O"   = "red",
										  "T"   = "purple",
										  "C"   = "blue",
										  "CkN" = "cyan"
							   )
							   ) +
			# Add theoretical density lines per distribution
			stat_function(fun = dOk1, args = list(N=N,mu_O=mu_O, sigma_O=sigma, minResponse=minResponse, maxResponse=maxResponse, type="GP"), aes(color = "Ok1"), size = 1) +
			stat_function(fun = dO, args = list(mu_O=mu_O, sigma_O=sigma, minResponse=minResponse, maxResponse=maxResponse, type="GP"), aes(color = "O"), size = 1) +
			stat_function(fun = dT, args = list(mu_O=mu_O, sigma_O=sigma, mu_D=mu_D, minResponse=minResponse, maxResponse=maxResponse, type="GP"), aes(color = "T"), size = 1) +
			stat_function(fun = dC, args = list(mu_O=mu_O, sigma_O=sigma, mu_D=mu_D, minResponse=minResponse, maxResponse=maxResponse, type="GP"), aes(color = "C"), size = 1) +
			stat_function(fun = dCkN, args = list(N=N,mu_O=mu_O,sigma_O=sigma, mu_D=mu_D, minResponse=minResponse, maxResponse=maxResponse, type="GP"), aes(color = "CkN"), size = 1) +
			xlim(xlim[1], xlim[2])+

			geom_vline(xintercept = xIs, color = "black") +

			labs(title = NULL,
				 x = responseVariableName, y = "Density") +
			theme_minimal()
	}
	else {

		cat("Making theoretical shaded plots.\n")
		custom_colors  =  c(
						   Ok1_d = "yellow",
						   O_d   = "red",
						   T_d   = "purple",
						   C_d   = "blue",
						   CkN_d = "cyan"
		)

		x = seq(xlim[1], xlim[2], length.out=1000)

		Ok1_d = dOk1(x,N=N,mu_O=mu_O,sigma_O=sigma,minResponse=minResponse, maxResponse=maxResponse, type="GP")
		O_d = dO(x, mu_O=mu_O, sigma_O=sigma, minResponse=minResponse, maxResponse=maxResponse, type="GP")
		T_d = dT(x, mu_O=mu_O, sigma_O=sigma, mu_D=mu_D, minResponse=minResponse, maxResponse=maxResponse, type="GP")
		C_d = dC(x, mu_O=mu_O, sigma_O=sigma, mu_D=mu_D, minResponse=minResponse, maxResponse=maxResponse, type="GP")
		CkN_d = dCkN(x,N=N,mu_O=mu_O,sigma_O=sigma, mu_D=mu_D, minResponse=minResponse, maxResponse=maxResponse, type="GP")

		df  =  tibble(x = x, Ok1_d,O_d,T_d,C_d,CkN_d) %>%
			pivot_longer(cols = c(Ok1_d,O_d,T_d,C_d,CkN_d), names_to = "curve", values_to = "y")

		p = ggplot(df, aes(x = x, y = y, fill = curve)) +
			geom_area(alpha = 0.4, position = "identity") +
			geom_line(aes(color = curve), size = 0.5) +
			scale_fill_manual(values = custom_colors) +
			geom_vline(xintercept = xIs, color = "black") +
			scale_color_manual(values = custom_colors) +
			labs(title = NULL,
				 x = responseVariableName, y = "Density", fill = "Curve", color = "Curve") +
			theme_minimal() +
			theme(
				  legend.background = element_rect(fill = "white", color = "black")
			)


	}
	if(!legend) {
		p = p+	theme(legend.position = "none")
	}
	return(p)
}
