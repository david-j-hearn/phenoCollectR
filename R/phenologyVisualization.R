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

getWeightedColors = function(stageColors, simulatedIndividuals) {
  n = length(simulatedIndividuals)
  nStages = length(simulatedIndividuals[[1]]$stageCounts)
  for(i in 1:n) {
        simulatedIndividuals[[i]]$stageCounts = tabulate(simulatedIndividuals[[i]]$phenologicalUnits, nbins=nStages)
  }
  rgbCols = col2rgb(stageColors)
  stageCountMatrix <- t(
    sapply(simulatedIndividuals, `[[`, "stageCounts")
  )

  rgbIndividuals <- (stageCountMatrix %*% t(rgbCols)) / rowSums(stageCountMatrix)
  cols <- rgb(rgbIndividuals[,1], rgbIndividuals[,2], rgbIndividuals[,3], maxColorValue = 255)
  return(cols)
}

#' Make a plot of simulated multistage data
#'
#' @description Make a plot of simulated multistage data
#'
#' @param simulatedData The output object from the function simulateMultistageData. (default: NULL)
#' @param targetCovariateIndex The index number of the covariate to be plotted along the x-axis. (default: 1)
#' @param stageColors A vector of colors, one per stage. (default: NULL)
#' @param drawIndividuals Flag to draw horizontal line segments representing individuals through time. (default: FALSE)
#' @param drawObserved Flag to draw observed collection times. (default: TRUE)
#' @param drawLatentOnset Flag to draw the true, but latent (unobserved) individual onset times. (default: TRUE)
#' @param drawTrueModel Flag to draw the true mean onset model as a function of the target covariate. (default: TRUE)
#' @param drawInferredModel Flag to draw the inferred mean onset model as a function of the target covariate BASED ON UNOBSERVED LATENT ONSET TIMES, not inferred from observed data. Can be useful to see how far off the fit to the true data is from the true model. (default: FALSE)
#' @param shadeStage Flag to shade the stages with transparent shade-specific color. (default: TRUE)
#' @param observedCex Set the plot point size of the observed collection times. (default: 0.75)
#' @param onsetCex Set the plot point size of the latent onset times. (default: 0.05)
#' @param observedPch Set the point display type of the observed collection times. (default: 4, an 'x')
#' @param onsetPch Set the plot point display type of the latent onset times. (default: 17, a triangle)
#' @param minResponse The minimum response time for setting plot bounds. (default: 0)
#' @param maxResponse The maximum response time for setting plot bounds. (default: 365)
#'
#' @return A plot illustrating the simulated multistage data for one target covariate along the x-axis.
#' @export
#'
#' @examples
#' \donttest{
#' ##Set basic simulation parameters
#' numberStages = 3
#' numberCovariates = 1
#' sampleSize = 50
#'
#' ##Simulate data
#' simulatedData = simulateMultistageData(n=sampleSize, nStages=numberStages, nCovariates=numberCovariates)
#'
#' ##Plot simulated data with the first covariate (the only covariate in this case) along the x-axis.
#' #    set colors for the stages
#' stageColors = viridisLite::viridis(numberStages)
#' # stageColors = rainbow(numberStages)
#' #    create the plot
#' plotMultistageSimulation(simulatedData=simulatedData, targetCovariateIndex=1, stageColors=stageColors)
#' }
plotMultistageSimulation = function(simulatedData=NULL, targetCovariateIndex=1, stageColors=NULL, drawIndividuals=FALSE, drawObserved=TRUE, drawLatentOnset=TRUE, drawTrueModel=TRUE, drawInferredModel=FALSE, observedCex=0.75, onsetCex=0.05, observedPch=4, onsetPch=17, shadeStage=TRUE, minResponse=0, maxResponse=365, main=NULL, includeOverlap=FALSE) {


  if(includeOverlap) {
    if(is.null(simulatedData$simulatedData)) {
      stop("Please provide the output from simulateMultistageOverlapData when includeOverlap is set to TRUE. Otherwise, provide the output from simulateMultistageData when not.")
    }
    nStages = simulatedData$simulatedData$nStages
    if(is.null(stageColors) || length(stageColors)<nStages) {
	    stageColors = rainbow(nStages)
    }

    weightedStageColors = getWeightedColors(stageColors=stageColors, simulatedIndividuals=simulatedData$simulatedIndividuals)
    simulatedData = simulatedData$simulatedData
  }
  else {
    nStages = simulatedData$nStages
    if(is.null(stageColors) || length(stageColors)<nStages) {
	    stageColors = rainbow(nStages)
    }
  }

	if(!drawLatentOnset && !drawObserved) {
		stop("Nothing to plot in plotMultistageSimulation. Quitting.")
	}

	#Extract information
	if(is.null(simulatedData)) { stop("As input, provide the simulated data from the function simulateMultistageData.") }
	if(is.null(targetCovariateIndex) || targetCovariateIndex<1) { stop("As input, provide the covariate number you wish to plot. Covariates are indexed from 1 to the number of covariates.") }
	if(is.null(stageColors)) { stop("Provide a vector with a named color for each stage in your multistage model.") }

	numberStages = simulatedData$nStages
    	if(is.null(stageColors) || length(stageColors)<nStages) {
	   	stageColors = rainbow(nStages)
    	}
	numberCovariates = simulatedData$nCovariates
	trueSlopes = rep(0,numberStages-1)		#one slope per stage for the target covariate, not including onset of stage 1 (at 0)
	trueIntercepts = rep(0,numberStages-1)		#one intercept per stage for the target covariate
	if(targetCovariateIndex>numberCovariates) { stop("Please provide the index number of the covariate you wish to plot between 1 and the number of covariates, inclusive.") }
	if(length(stageColors) < numberStages) { 
		stop("Please provide a color for each stage.") 
	}

	#One covariate, so intercept (with anchor adjustment) and slope stay the same
	if(numberCovariates==1) { 
		trueSlopes[1] = simulatedData$stage2OnsetCovariateSlopes[targetCovariateIndex]
		trueIntercepts[1] = simulatedData$stage2OnsetMean - sum(simulatedData$stage2OnsetCovariateSlopes * simulatedData$covariateMeans) #sum not neede here, for for consistency
		cumOnsetSum = simulatedData$stage2OnsetMean
		if(numberStages>2)
		{
		for(i in 1:(numberStages-2)) {
			#trueSlopes[i+1] = simulatedData$stageDurationCovariateSlopes[i,targetCovariateIndex]
			trueSlopes[i+1] = trueSlopes[i] + simulatedData$stageDurationCovariateSlopes[i,targetCovariateIndex]
			trueIntercepts[i+1] = cumOnsetSum + simulatedData$stageDurationMeans[i] - sum(simulatedData$stageDurationCovariateSlopes[i,] * simulatedData$covariateMeans)
			cumOnsetSum = cumOnsetSum + simulatedData$stageDurationMeans[i]
		}
		}
	}
	#Get marginal model parameters based on correlation structure of the covariates
	else { 
		modelStage2 = true_marginal_line(alpha = simulatedData$stage2OnsetMean - sum(simulatedData$stage2OnsetCovariateSlopes * simulatedData$covariateMeans), 
						 beta = simulatedData$stage2OnsetCovariateSlopes, 
						 mu = simulatedData$covariateMeans, 
						 j = targetCovariateIndex, 
						 Sigma = simulatedData$Sigma, 
						 R = simulatedData$R, 
						 sd_x = simulatedData$covariateSDs) 
		trueSlopes[1] = modelStage2$slope
		trueIntercepts[1] = modelStage2$intercept
		cumOnsetSum = simulatedData$stage2OnsetMean
		if(numberStages>2) {
		for(i in 1:(numberStages-2)) {
			model = true_marginal_line(alpha = simulatedData$stageDurationMeans[i] - sum(simulatedData$stageDurationCovariateSlopes[i,] * simulatedData$covariateMeans),
						   beta = simulatedData$stageDurationCovariateSlopes[i,],
						   mu = simulatedData$covariateMeans, 
						   j = targetCovariateIndex, 
						   Sigma = simulatedData$Sigma, 
						   R = simulatedData$R, 
						   sd_x = simulatedData$covariateSDs) 
			#trueSlopes[i+1] = model$slope
			trueSlopes[i+1] = model$slope + trueSlopes[i]
			trueIntercepts[i+1] = model$intercept + cumOnsetSum
			cumOnsetSum = cumOnsetSum + simulatedData$stageDurationMeans[i]
		}
		}
	}

	colInc = 1
	#if(nonCyclical) {
		#if(length(stageColors) < numberStages+1) { 
			#stop("For non-cyclical stage data, please provide an additional color.") 
		#}
		#colInc = 1
	#}

	#Begin the plotting
	x1 = min(simulatedData$outputData[,targetCovariateIndex])-50
	x2 = max(simulatedData$outputData[,targetCovariateIndex])+50
  if(includeOverlap) {
		plot(simulatedData$outputData[,targetCovariateIndex], simulatedData$outputData$sampledTime, col=weightedStageColors, pch=observedPch, cex=observedCex, ylim=c(minResponse,maxResponse), main=main, xlab="Target covariate value", ylab="Time of observation")
  }
	for(i in 1:(numberStages-1)) {
		if(i == 1) {
			if(!drawLatentOnset && !includeOverlap) {
				plot(simulatedData$outputData[,targetCovariateIndex], simulatedData$outputData$sampledTime, col=stageColors[simulatedData$outputData$sampledStage], ylim=c(minResponse,maxResponse), pch=observedPch, cex=observedCex, main=main, xlab="Target covariate value", ylab="Time of observation")

			}
			if(drawLatentOnset && !includeOverlap) {
				plot(simulatedData$outputData[,targetCovariateIndex],simulatedData$outputData[,numberCovariates+1],col=stageColors[1+colInc],ylim=c(minResponse,maxResponse),pch=onsetPch,cex=observedCex, main=main, xlab="Target covariate value", ylab="Time of observation")
			}
			if(shadeStage) {
				y1 = c(0,0)
				y2 = trueIntercepts[i] + trueSlopes[i] * c(x1,x2)
				y3 = trueIntercepts[i+1] + trueSlopes[i+1] * c(x1,x2)
				#Under stage 1
				#if(nonCyclical) {
					tc = adjustcolor(stageColors[1], alpha.f = 0.1)
				#}
				#else {
					#tc = adjustcolor(stageColors[numberStages], alpha.f = 0.1)
				#}
				polygon( x = c(x1, x2, x2, x1), 
					y = c(y1[1], y1[2], y2[2], y2[1]), 
					col = tc,
					border = NA
				)
				if(i==numberStages-1) {		#case when 2 stages
					y1 = trueIntercepts[i] + trueSlopes[i] * c(x1,x2)
					y2 = c(maxResponse,maxResponse)
					tc = adjustcolor(stageColors[i+colInc], alpha.f = 0.1)
					polygon( x = c(x1, x2, x2, x1), 
						y = c(y1[1], y1[2], y2[2], y2[1]), 
						col = tc,
						border = NA
					)
				}
				##Above stage 1 onset
				#tc = adjustcolor(stageColors[1+colInc], alpha.f = 0.1)
				#polygon( x = c(x1, x2, x2, x1), 
					#y = c(y2[1], y2[2], y3[2], y3[1]), 
					#col = tc,
					#border = NA
				#)
			}
		}
		else {
			if(drawLatentOnset && !includeOverlap) {
				points(simulatedData$outputData[,targetCovariateIndex],simulatedData$outputData[,numberCovariates+i],col=stageColors[i+colInc],pch=onsetPch,cex=observedCex)
			}
			if(shadeStage) {
				if(i==numberStages-1) {
					y1 = trueIntercepts[i] + trueSlopes[i] * c(x1,x2)
					y2 = c(maxResponse,maxResponse)
					tc = adjustcolor(stageColors[i+colInc], alpha.f = 0.1)
					polygon( x = c(x1, x2, x2, x1), 
						y = c(y1[1], y1[2], y2[2], y2[1]), 
						col = tc,
						border = NA
					)
				}
				#else {
					y1 = trueIntercepts[i-1] + trueSlopes[i-1] * c(x1,x2)
					y2 = trueIntercepts[i] + trueSlopes[i] * c(x1,x2)
					#y1 = trueIntercepts[i] + trueSlopes[i] * c(x1,x2)
					#y2 = trueIntercepts[i+1] + trueSlopes[i+1] * c(x1,x2)
					tc = adjustcolor(stageColors[i], alpha.f = 0.1)
					polygon( x = c(x1, x2, x2, x1), 
						y = c(y1[1], y1[2], y2[2], y2[1]), 
						col = tc,
						border = NA
					)
				#}
			}
		}
		if(drawInferredModel) {
			fit = lm(simulatedData$outputData[,numberCovariates+i] ~ simulatedData$outputData[,targetCovariateIndex])
			abline(fit,col=stageColors[i+colInc], lty="dashed")
		}
		if(drawTrueModel) {
			print(paste("a= ", trueIntercepts[i], ", b= ", trueSlopes[i]))
			abline(a=trueIntercepts[i], b=trueSlopes[i], col=stageColors[i+colInc])
		}
	}
	if(drawLatentOnset && drawObserved && !includeOverlap) {
		points(simulatedData$outputData[,targetCovariateIndex], simulatedData$outputData$sampledTime, col=stageColors[simulatedData$outputData$sampledStage], pch=observedPch, cex=observedCex)
	}
	if(drawIndividuals) {
		tc = adjustcolor("gray25", alpha.f = 0.1)
		segments(x0=simulatedData$outputData[,targetCovariateIndex],x1=simulatedData$outputData[,targetCovariateIndex],y0=minResponse,y1=maxResponse,col=tc)
	}

	return(list(trueIntercepts = trueIntercepts, trueSlopes = trueSlopes))
}

#' Make a posterior predictive plot from a "full" presence-only dataset analysis focused on one covariate and marginalized across the other covariates
#'
#' @description Make a graphic of the posterior predictive values of the extremes, onset, observed, and cessation values as well as sub-plots indication directions of change for onset, duration, and observed collection times. 
#'
#' @param stanResult  The output from the function 'runStanPhenology' with a 'full' model
#' @param responseData A vector of the response data
#' @param responseVariableName The name of the response variable (default: DOY)
#' @param targetCovariateName A string with the name of the target covariate to be graphed on the x-axis of the output plot
#' @param onsetCovariateData A data frame with the onset covariate data. Can be obtained from the function 'preparePhenologyData'.
#' @param durationCovariateData A data frame with the duration covariate data. Can be obtained from the function 'preparePhenologyData'.
#' @param n_samples The number of covariate values to simulate for each value of the target covariate. Higher values give more accurate and more precise estimates. (default: 100)
#' @param n_draws The number of draws from the posterior distribution. Higher values give more accurate and more precise estimates. (default: 100)
#' @param n_hist  The sample size of values to be used to approximate the onset, duration, and observed distributions. (default: 10000)
#' @param n_bin The number of bins to be used for the histograms. (default: 100)
#' @param resolution The number of equi-distant values of the target covariate to be used for graphing. Higher values give smoother curves, but take longer to render. (default: 50)
#' @param slice The relative position where time slices will be made (also at 1 - slice). Plots of phenological distributions are made at the time slices. (default: 0.25)
#' @param N The population size used to estimate phenological extremes. (default: 1000)
#' @param minResponse The minimum possible value of the response variable. (default: 0)
#' @param maxResponse The maximum possible value of the response variable. (default: 365)
#' @param adjustYLim Boolean flag specifying whether to set the limits of the y-axis adaptively. Use default. (default: TRUE)
#' @param smooth Define the type of smoothing, if any, to be used in the output plot. Possibilities are GAM, LOESS, SPLINE, and NONE. (default: GAM)
#' @param keepSigma Boolean flag indicating whether the sigma parameter should be included in the graphical output. (default: FALSE)
#' @param makeHistograms Boolean flag indicating whether histograms (vs. densities) at the time slices should be made. (default: FALSE)
#'
#' @return A ggplot2 graphics object
#' @export
#' @importFrom posterior as_draws_df
#' @importFrom copula pobs normalCopula fitCopula
#' @importFrom dplyr tibble bind_rows group_by filter ungroup n row_number group_modify %>% 
#' @importFrom mgcv gam
#' @importFrom stats predict
#' @importFrom splines bs
#' @importFrom cowplot plot_grid
#' @importFrom tibble as_tibble
#' @importFrom stats quantile sd rnorm 
#'
#' @examples
#' \donttest{
#' ##get the file name with data for the blood root plant
#' file  =  getDatasetPath("Sanguinaria_canadensis")
#' ## See documentation for more species:
#' help(getDatasetPath)
#' ##define the covariate names - remove up to all but 1
#' vars = c("Latitude", "Year", "Elevation", "AnnualMonthlyAverageTemp", "SpringMonthlyAverageTemp"
#'          , "FirstQuarterMonthlyAverageTemp")
#' ##get the phenology data
#' data  =  preparePhenologyData(dataFile=file, responseVariableName="DOY"
#'                               , onsetCovariateNames=vars, durationCovariateNames=vars
#'                               , taxonName="Sanguinaria_canadensis", removeOutliers=TRUE)
#' ##run the Stan sampler
#' stanResult  =  runStanPhenology(type="full", responseData = data$responseData
#'                                 , onsetCovariateData = data$onsetCovariateData
#'                                 , durationCovariateData = data$durationCovariateData
#'                                 , partitionDataForPriors = TRUE)
#' ##summarize the Stan run
#' stanSummary  =  summarizePhenologyResults(stanRunResult = stanResult
#'                                           , taxonName = "Sanguinaria_canadensis"
#'                                           , standardLinearModel = TRUE)
#' ##make posterior predictive graph
#' pp  =  makePosteriorPredictivePlot(stanResult = stanResult, responseData = data$responseData
#'                                    , targetCovariateName = "SpringMonthlyAverageTemp"
#'                                    , onsetCovariateData = data$onsetCovariateData
#'                                    , durationCovariateData = data$durationCovariateData)
#' ##display the posterior predictive graph
#' print(pp)
#' }
makePosteriorPredictivePlot = function(stanResult, responseData, responseVariableName="DOY", targetCovariateName, onsetCovariateData, durationCovariateData, n_samples=100, n_draws=100, n_hist=10000, n_bin=100, resolution=50, slice = 0.25, N=1000, minResponse=0, maxResponse=365, adjustYLim=TRUE, smooth=c("GAM", "LOESS", "SPLINE", "NONE"), keepSigma=FALSE, makeHistograms=FALSE) {

	smooth = match.arg(smooth)

  if(length(smooth)>1) {
    smooth = "GAM"
  }

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
			s1_raw = pred_means
		}
		if(x == s2) {
			s2_summary = pred_summary
			s2_raw = pred_means
		}

		results  =  bind_rows(results, pred_summary)
	}

	minY = min(pred_summary$lower)
	maxY = max(pred_summary$upper)

	if(!makeHistograms) {

		O1s0 = s1_raw["mO",]
		D1s0 = s1_raw["mC",]- s1_raw["mO",]
		T1s0 = s1_raw["mT",]
		S1s0 = s1_raw["S",]

		#print("O1s")
		#print(O1s)

		O2s0 = s2_raw["mO",]
		D2s0 = s2_raw["mC",]- s2_raw["mO",]
		T2s0 = s2_raw["mT",]
		S2s0 = s2_raw["S",]

		cat("Simulating shift distributions.\n")
		n = 100

		cat("\tOnset\n")
		O1list = list()
		O2list = list()
		for( i in 1:length(O1s0))
		{
			O1list[[i]] = rO(n=n, mu_O=O1s0[i], sigma_O=S1s0[i], type="GP")
			O2list[[i]] = rO(n=n, mu_O=O2s0[i], sigma_O=S2s0[i], type="GP")
		}
		O1s = unlist(O1list)
		O2s = unlist(O2list)

		cat("\tDuration\n")
		D1list = list()
		D2list = list()
		sd1 = sd(D1s0)
		sd2 = sd(D2s0)
		for( i in 1:length(D1s0))
		{
			#includes intrinsic variation in D, which is not part of the GP model
			D1list[[i]] = rD.GP(n=n, mu_D=D1s0[i]) + rnorm(n=n, 0, sd1)
			D2list[[i]] = rD.GP(n=n, mu_D=D2s0[i]) + rnorm(n=n, 0, sd2)
		}
		D1s = unlist(D1list)
		D2s = unlist(D2list)

		cat("\tObserved\n")
		T1list = list()
		T2list = list()
		for( i in 1:length(T1s0))
		{
			T1list[[i]] = rT(n=n, mu_O=O1s0[i], sigma_O=S1s0[i], mu_D=D1s0[i], type="GP")
			T2list[[i]] = rT(n=n, mu_O=O2s0[i], sigma_O=S2s0[i], mu_D=D2s0[i], type="GP")
		}
		T1s = unlist(T1list)
		T2s = unlist(T2list)

		cat("Making shift panels.\n")
		pO = makeShiftPanel(O1s, O2s,"red", "Mean Onset (O)")
		pD = makeShiftPanel(D1s, D2s, "gray", "Mean Duration (D)")
		pT = makeShiftPanel(T1s, T2s, "purple", "Mean Observed (T)")

		slices = cowplot::plot_grid(pO, pD, pT, nrow = 3, ncol = 1, labels = c("B", "C", "D"), label_x = 0, label_y = 0.1)
		#print(pO)
		#dev.new()
		#print(slices)
	}

	#make histograms of phenological distributions at two different time slices
	else {

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
		slices = cowplot::plot_grid(pS1, pS2, nrow = 2, ncol = 1, labels = c("B", "C"), label_x = 0, label_y = 0.1)
	}


	#get rid of the sigma data for graphing
	if(!keepSigma) {
		results  =  results[results$series != "Sigma", ]
	}

	#The last data point appears to be highly biased (error in sampling from copula?), so graphing all but the last point
	## DANIEL: Issue with object "series" called below without prior definition:
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

	#slices = plot_grid(pS1, pS2, nrow = 2, ncol = 1, labels = c("B", "C"), label_x = 0, label_y = 0.1)
	figure = plot_grid(pPPP, slices, nrow = 1, ncol = 2, labels = c("A"), label_x = 0, label_y = 0.1, rel_widths = c(2, 1))

	return(figure)
}

#' @importFrom ggplot2 ggplot geom_ribbon aes geom_line geom_point scale_color_manual scale_fill_manual geom_vline theme_minimal labs theme element_rect ylim coord_cartesian
posterior_predictive_graphic = function(observed_data, filtered_results,targetCovariateName,responseVariableName,legend=TRUE, timeSlice1, timeSlice2, adjustYLim = TRUE, minResponse=0, maxResponse=365) {

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

#' @importFrom ggplot2 ggplot aes  geom_vline labs theme_minimal geom_area theme annotate
#' @importFrom grid arrow unit
#' @importFrom dplyr tibble group_by bind_rows summarize %>%
#' @importFrom stats density
makeShiftPanel = function(sample1, sample2, col, xlab) {

	df <- bind_rows(
			tibble(value = sample1, timepoint = "T1"),
			tibble(value = sample2, timepoint = "T2")
	)

	# Compute sample means
	means <- df %>%
		group_by(timepoint) %>%
		summarize(mean_value = mean(value))

	x_range <- range(df$value)
	buffer <- 0.25 * diff(x_range)
	x_seq <- seq(x_range[1] - buffer, x_range[2] + buffer, length.out = 512)

	density_T1 <- density(df$value[df$timepoint == "T1"], from = min(x_seq), to = max(x_seq), n = 1000)
	density_T2 <- density(df$value[df$timepoint == "T2"], from = min(x_seq), to = max(x_seq), n = 1000)

	# Convert to data frames
	dens_T1_df <- tibble(x = density_T1$x, y = density_T1$y, timepoint = "T1")
	dens_T2_df <- tibble(x = density_T2$x, y = density_T2$y, timepoint = "T2")

	densities_df <- bind_rows(dens_T1_df, dens_T2_df)

	ymax <- max(densities_df$y)
	ymid <- ymax / 2

	ggplot(densities_df, aes(x = x, y = y)) +
		#geom_area(fill = col, alpha = 0.1, position = "identity") +
		geom_area(data = dens_T1_df, aes(x = x, y = y),
			  fill = col, alpha = 0.3, color = "black") +
		geom_area(data = dens_T2_df, aes(x = x, y = y),
			  fill = col, alpha = 0.3, color = "black") +
		geom_vline(data = means, aes(xintercept = mean_value), linetype = "dashed", color = "black") +
		annotate("segment",
			 x = means$mean_value[means$timepoint == "T1"],
			 xend = means$mean_value[means$timepoint == "T2"],
			 y = ymid, yend = ymid,
			 arrow = arrow(length = unit(0.25, "cm"), type="closed"),
			 color = "black") +
		labs( x = xlab, y = "Density") +
		theme_minimal() +
		theme(legend.position = "none")

}


#' @importFrom ggplot2 ggplot geom_histogram aes scale_fill_manual scale_color_manual stat_function xlim geom_vline labs theme_minimal geom_area geom_line theme element_rect
#' @importFrom dplyr tibble %>%
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

#' Make a posterior predictive plot from a "multistage-full" dataset analysis focused on one covariate
#' @description Make a graphic of the posterior predictive values of the onsets of multiple stages including 95% credible intervals marginalized for a single covariate displayed along the x-axis.
#'
#' @param stanResult The output from the function 'runStanPhenology' with a 'multistge-full' model
#' @param responseData A vector of the response data
#' @param stageData A vector of the stage for each item in the responseData vector
#' @param responseVariableName The name of the response variable (default: DOY)
#' @param targetCovariateName A string with the name of the target covariate to be graphed on the x-axis of the output plot
#' @param covariateData A data frame with the onset covariate data. Columns are labeled with covariate names and rows have the covariate data. Rows should correspond to rows in the responseData vector.
#' @param nReps The number of replicates to be used during the Monte Carlo simulation for each value of the target covariate. (default: 100)
#' @param nXs The number of intervals to be used for the target covariate. (default: 100)
#' @param nStages The number of stages. 
#' @param slice he relative position where time slices will be made (also at 1 - slice). Plots of phenological distributions are made at the time slices. (default: 0.25)
#' @param minResponse Must be 0. (default: 0)
#' @param maxResponse The maxiumum possible value of the observed collection times. (default: 365)
#' @return a ggplot2 graphics object
#' @export
#' @importFrom posterior as_draws_df
#' @importFrom tibble as_tibble
#' @importFrom stats quantile sd rnorm
#' @examples
#' \donttest{
#' ##Run a full multistage analysis with SIMULATION PARAMETERS themselves simulated and SIMULATED COVARIATES with MULTIPLE STAGES followed by a posterior predictive plot of the results
#' #
#' library(dplyr)
#' library(posterior)
#'  #Set the parameters
#'  nStages = 4
#'  nCovariates = 3
#'  nXs = 101
#'  nReps = 100
#'  n = 500
#'  minResponse = 0             #The minimum observed time should always be at least 0
#'  maxResponse = 365
#'  simulatedData = simulateMultistageData(n=n, nStages=nStages,nCovariates=nCovariates)
#'
#'  #Plotting
#'      #Set which covariate is the x-axis in the plot
#'  targetCovariateIndex = 1
#'
#'      #Set colors
#'  #stageColors = viridisLite::viridis(nStages+1)
#'  stageColors = RColorBrewer::brewer.pal(nStages+1, "Set1")
#'
#'  trueModels = plotMultistageSimulation(simulatedData=simulatedData,
#'                   targetCovariateIndex=targetCovariateIndex,
#'                   stageColors=stageColors,
#'                   drawLatentOnset=FALSE,
#'                   drawObserved=TRUE,
#'                   drawTrueModel=FALSE,
#'                   drawInferredModel=FALSE,
#'                   shadeStage=TRUE,
#'                   minResponse=minResponse, maxResponse=maxResponse)
#'  dev.new()
#'  #Run inference with Stan
#'  stanResult =runStanPhenology(
#'      type="multistage-full",                             #model with many stages and covariates
#'      responseData=simulatedData$outputData$sampledTime,  #observed collection times
#'      stage=simulatedData$outputData$sampledStage,        #observed stage at collection time
#'      nStages=nStages,                                            #number of stages
#'      nOnsetCovariates=nCovariates,                       #number of onset covariates (same as duration)
#'      nDurationCovariates=nCovariates,                     #number of duration covariates (same as onset)
#'      onsetCovariateData=simulatedData$X,                 #covariate data (same as for duration)
#'      durationCovariateData=simulatedData$X,              #covariate data (same as for onset)
#'	nReps=nReps,					    #MC replicates for posterior predictive distribution
#'	nXs=nXs,					    #Number of intervals on x-axis of posterior predictive plot
#'      maxDiv=4000                                         #should be set to 0, but in case of the stray divergence, set high for this example
#'  )
#' ###Make the posterior predictive plot
#' p = makeMultistagePosteriorPredictivePlot(stanResult=stanResult, 
#'			responseData=simulatedData$outputData$sampledTime, 
#'			targetCovariateName="cov1", 
#'			covariateData=simulatedData$X, 
#'			stageData=simulatedData$outputData$sampledStage, 
#'			nReps=nReps,	#this is the default for runStanPhenology
#'			nXs=nXs,	#this is the default for runStanPhenology
#'			smooth="GAM",
#'			nStages=nStages+1)
#'  
#'  for(i in 1:nStages) {
#'  	p = p + geom_abline(intercept = trueModels$trueIntercepts[i], slope = trueModels$trueSlopes[i], color = "black", linewidth = 0.5, linetype="dashed")
#'  }
#'  #Extract basic summary data
#'  probs = c(2.5, 97.5)
#'  measures=c("mean", "median", "sd", "mad")
#'
#'  summary = print(
#'    stanResult$result$summary(
#'      variables = c(
#'        "sigma","anchor_d", "beta_d", "alpha_d", "anchor_o", "beta_o", "alpha_o"
#'       ),
#'      quantiles = ~ quantile2(., probs=probs/100),
#'      measures
#'    ),
#'    n = Inf
#'  ) %>% as.data.frame()
#'
#'  #Set up vectors to store mean durations (means) and mean onsets (means_o) with low and high bounds of credible intervals.
#'  means = rep(0,nStages+1)
#'  mean_low = rep(0,nStages)
#'  mean_high = rep(0,nStages)
#'
#'  #Extract slopes and intercepts (anchors) and overlay the inferred lines onto earlier plot
#'  for(j in 1:(nStages+1)) {
#'          beta = rep(0,nCovariates)
#'          beta_low = rep(0,nCovariates)
#'          beta_high = rep(0,nCovariates)
#'          varO = paste0("anchor_o[", j, "]")
#'          if(j<=nStages) {
#'                  varD = paste0("anchor_d[", j, "]")
#'                  means[j] = summary[summary$variable==varD, "mean"]
#'                  mean_low[j] = summary[summary$variable==varD, paste0("q",probs[1])]
#'                  mean_high[j] = summary[summary$variable==varD, paste0("q",probs[2])]
#'          }
#'          alpha = summary[summary$variable==varO, "mean"]
#'          alpha_low = summary[summary$variable==varO, paste0("q",probs[1])]
#'          alpha_high = summary[summary$variable==varO, paste0("q",probs[2])]
#'
#'         for(i in 1:nCovariates) {
#'                  varO = paste0("beta_o[", j, ",", i, "]")
#'                  beta[i] = summary[summary$variable==varO, "mean"]
#'                  beta_low[i] = summary[summary$variable==varO, paste0("q",probs[1])]
#'                  beta_high[i] = summary[summary$variable==varO, paste0("q",probs[2])]
#'          }
#'         model = phenoCollectR:::true_marginal_line(alpha=alpha, beta=beta, mu=simulatedData$covariateMeans, Sigma=simulatedData$Sigma, j=targetCovariateIndex)
#'         model_low = phenoCollectR:::true_marginal_line(alpha=alpha_low, beta=beta_low, mu=simulatedData$covariateMeans, Sigma=simulatedData$Sigma, j=targetCovariateIndex)
#'         model_high = phenoCollectR:::true_marginal_line(alpha=alpha_high, beta=beta_high, mu=simulatedData$covariateMeans, Sigma=simulatedData$Sigma, j=targetCovariateIndex)
#'  	p = p + geom_abline(intercept = model$intercept, slope = model$slope, color = "red", linewidth = 0.5, linetype="dashed")
#'  }
#' print(p)
#' }
makeMultistagePosteriorPredictivePlot = function(stanResult, responseData, responseVariableName="DOY", targetCovariateName, covariateData, stageData, nReps=100, nXs=101, nStages, slice = 0.25, minResponse=0, maxResponse=365, y_pred=FALSE, nChains=4, nIts=1000, smooth=c("GAM", "LOESS", "SPLINE", "NONE")) {

	smooth = match.arg(smooth)

  if(length(smooth)>1) {    #shouldn't be reached...
    smooth = "GAM"
  }

	stageNames = paste0("stage", 1:(nStages))
	onsetMeans = as.data.frame( matrix(NA_real_, nrow = nXs-1, ncol = (nStages+1)))		#copula / smoothing doesn't do well with last x coordinate
	onsetQ2.5 = as.data.frame( matrix(NA_real_, nrow = nXs-1, ncol = (nStages+1)))		#copula / smoothing doesn't do well with last x coordinate
	onsetQ97.5 = as.data.frame( matrix(NA_real_, nrow = nXs-1, ncol = (nStages+1)))		#copula / smoothing doesn't do well with last x coordinate
	names(onsetMeans) <- c("x",stageNames)
	names(onsetQ2.5) <- c("x",stageNames)
	names(onsetQ97.5) <- c("x",stageNames)
	#Storage for output rows are stages, columns are x values, replicates are slices
	onsetsAll = numeric(nReps)
	stageMeans = numeric(nStages)
	stageQ2.5 = numeric(nStages)
	stageQ97.5 = numeric(nStages)
	sigma = 0

	#set the covariate information
	nCovariates = ncol(covariateData)
	Sigma = cov(covariateData)
	covariateMeans = colMeans(covariateData)
	target_idx = which(colnames(covariateData) == targetCovariateName)

	#get the minimum and maximum target covariate values and set the x increments
	minXt = min(covariateData[targetCovariateName])
	maxXt = max(covariateData[targetCovariateName])
	rangeXt = maxXt - minXt
	Xts = seq(from = minXt, to = maxXt, length.out=nXs)

	#set the original sample size
	N = nrow(covariateData)
	sqrtN = sqrt(N)

	#get draws
	posterior_samples  =  as_draws_df(stanResult$result$draws())
	posterior_samples = as_tibble( posterior_samples ) ## DANIEL: Try to prevent warning messages

	if(y_pred){

    cat("Using posterior predictive values generated by Stan.\n")

		totIts = nReps * nChains * nIts
		stageVals = numeric(totIts)
		for(i in 1:(nXs-1)) {
			low = (i-1)*nReps+1
			high = i*nReps
			for(s in 1:nStages) {
				predInds = paste0("y_pred[",low:high , ",", s , "]")
				vals = as.numeric(as.matrix(posterior_samples[, predInds]))
				stageMeans[s] = mean(vals)
				stageQ2.5[s] = quantile(vals, 0.025)
				stageQ97.5[s] = quantile(vals, 0.975)
			}
			onsetMeans[i,] = as.numeric(c(Xts[i],stageMeans))
			onsetQ2.5[i,] = as.numeric(c(Xts[i],stageQ2.5))
			onsetQ97.5[i,] = as.numeric(c(Xts[i],stageQ97.5))
		}
	if(is.na(smooth)) {   #shouldn't be reached
		stop("Missing the smoothing type.")
	}
	else if(smooth == "GAM") {
		if(nXs <=10) {
			stop("Please increase x-axis resolution (nXs) to greater than 10 before attempting to use GAM to smooth the posterior predictive plot.")
		}
    #print(stageNames)
		for(name in stageNames) {
        onsetMeans[[name]] = smooth_safe(onsetMeans[[name]], onsetMeans$x)
        onsetQ2.5[[name]] = smooth_safe(onsetQ2.5[[name]], onsetQ2.5$x)
        onsetQ97.5[[name]] = smooth_safe(onsetQ97.5[[name]], onsetQ97.5$x)
		    #onsetMeans[[name]]  = predict(gam(onsetMeans[[name]]  ~ s(onsetMeans$x), method = "REML"))
		    #onsetQ2.5[[name]]  = predict(gam(onsetQ2.5[[name]]  ~ s(onsetQ2.5$x), method = "REML"))
		    #onsetQ97.5[[name]]  = predict(gam(onsetQ97.5[[name]]  ~ s(onsetQ97.5$x), method = "REML"))
		}
	}
	else if(smooth == "LOESS") {
		for(name in stageNames) {
						    onsetMeans[[name]]  = predict(loess(onsetMeans[[name]]  ~ onsetMeans$x, span = 0.3))
						    onsetQ2.5[[name]] = predict(loess(onsetQ2.5[[name]]  ~ onsetQ2.5$x, span = 0.3))
						    onsetQ97.5[[name]] = predict(loess(onsetQ97.5[[name]] ~ onsetQ97.5$x, span = 0.3))
    }
	}
	else if(smooth == "SPLINE") {
		for(name in stageNames) {
						    onsetMeans[[name]]  = predict(lm(onsetMeans[[name]]  ~ bs(onsetMeans$x, df = 5)))
						    onsetQ2.5[[name]] = predict(lm(onsetQ2.5[[name]] ~ bs(onsetQ2.5$x, df = 5)))
						    onsetQ97.5[[name]] = predict(lm(onsetQ97.5[[name]] ~ bs(onsetQ97.5$x, df = 5)))
    }
	}
		p = make_multistage_PPD_plot(onsetMeans, onsetQ2.5, onsetQ97.5, rainbow(nStages), responseData, covariateData[targetCovariateName], stageData, minResponse, maxResponse)
		return(p)
	}

	else {

		if(minResponse!=0) {
			stop("The minimum possible observed time (minResponse) must be kept at the default of 0.")
		}

    cat("Generating posterior predictive values. This may take some time.\n")

		#Draw nReps * nXs theta from Stan posterior draws - this can be done without knowledge of x, can be done outside of loop
		# 	Get sampled indices
		indices = sample(x=1:nrow(posterior_samples), size=nReps*nXs, replace=TRUE)
		# 	Extract alpha (intercept)
		alphaInds = paste0("alpha_o[", 1:nStages, "]")
		# 	Extract beta vector (assumes beta[stage,covariate])
		betaInds  =  beta = paste0("beta_o[", rep(1:nStages, each = nCovariates), ",", rep(1:nCovariates, times = nCovariates), "]")
		# 	Named sigma variable
		sigmaInds = "sigma"
		#	Extract samples
		#		this includes alphas for each stage, betas for each stage and each covariate, and sigma, which is the same for all stages
		theta = posterior_samples[indices, c(betaInds,alphaInds,sigmaInds)]
		print(head(theta))

		#Sample nReps * nXs * nStages onset times from normal distribution - the normal deviates do not need knowledge of x, mu, or sd, can be done outside of loop
		normalDeviatesO = rnorm(nReps*nXs*nStages)
		#normalDeviatesSE = rnorm(nReps*nXs)

		#sqrtNReps = sqrt(nReps)

		#for each x value of a target covariate Ct iterated across the range of Ct 
		cnt = 1
		cntX = 1
		for(x in Xts) {
			print(x)

			#Calculate the mean onsets of all stages for all replicates using the marginal model parameters and sampled X
			#	get nReps samples for each x
			#print("going through stages")
			for(s in 1:nStages) {
				cntR = 1
				#print("going through replicates")
				for(n in ((cntX-1)*nReps + 1):(cntX*nReps)) {
					beta = as.numeric(theta[,paste0("beta_o[", s, ",", 1:nCovariates, "]")][n,])
					alpha = as.numeric(theta[,paste0("alpha_o[", s, "]")][n,])
					sigma = as.numeric(theta[,c("sigma")][n,])
					m = true_marginal_line(alpha=alpha, beta=beta, mu=covariateMeans, j=target_idx, Sigma=Sigma)
					#population level, not generated per individual
					#print("")
					#print(paste0("alpha_o[", s, "]"))
					#print(alpha)
					#print(paste0("beta_o[", s, ",", 1:nCovariates, "]"))
					#print(beta)
					#print("marginal model")
					#print(m$intercept)
					#print(m$slope)
					#print("x")
					#print(x)
					#print("rep")
					#print(n)
					#print("stage")
					#print(s)
					#print("sigma")
					#print(sigma)
					#print(normalDeviatesO[cnt])
					if(s == 1) {
						onsetsAll[cntR] = m$intercept + m$slope * x
					}
					else {
						onsetsAll[cntR] = m$intercept + m$slope * x + sigma * normalDeviatesO[cnt]
					}
					cntR = cntR+1
					cnt = cnt+1
				}
				#Add error to the mean based on sample size and standard deviation (replaces sample from covariates and bases on mean marginal expectation)
				#stageMeans[s] = mean(onsetsAll) + sd(onsetsAll) / sqrtNReps
				#Calculate stage-specific means across replicates - use rowMeans, as C vectorized and faster
				#print("Setting stage means.")
				stageMeans[s] = mean(onsetsAll)
				#Calculate 2.5% and 97.5% quantiles bases on the replicates
				stageQ2.5[s] = quantile(onsetsAll, 0.025)
				stageQ97.5[s] = quantile(onsetsAll, 0.975)
			}
			#print("Setting onset means.")
			#print(c(x,stageMeans))
			onsetMeans[cntX,] = as.numeric(c(x,stageMeans))
			onsetQ2.5[cntX,] = as.numeric(c(x,stageQ2.5))
			onsetQ97.5[cntX,] = as.numeric(c(x,stageQ97.5))
			#print("Done setting onset means.")
			cntX = cntX + 1
		}
		print("Done with loops, now making plot of the following data")
		#print(onsetMeans)
	}
	if(is.na(smooth)) {
		stop("Missing the smoothing type.")
	}
	else if(smooth == "GAM") {
		if(nXs <=10) {
			stop("Please increase x-axis resolution (nXs) to greater than 10 before attempting to use GAM to smooth the posterior predictive plot.")
		}
		for(name in stageNames) {
        onsetMeans[[name]] = smooth_safe(onsetMeans[[name]], onsetMeans$x)
        onsetQ2.5[[name]] = smooth_safe(onsetQ2.5[[name]], onsetQ2.5$x)
        onsetQ97.5[[name]] = smooth_safe(onsetQ97.5[[name]], onsetQ97.5$x)
		    #onsetMeans[[name]]  = predict(gam(onsetMeans[[name]]  ~ s(onsetMeans$x), method = "REML"))
		    #onsetQ2.5[[name]]  = predict(gam(onsetQ2.5[[name]]  ~ s(onsetQ2.5$x), method = "REML"))
		    #onsetQ97.5[[name]]  = predict(gam(onsetQ97.5[[name]]  ~ s(onsetQ97.5$x), method = "REML"))
		}
	}
	else if(smooth == "LOESS") {
		for(name in stageNames) {
						    onsetMeans[[name]]  = predict(loess(onsetMeans[[name]]  ~ onsetMeans$x, span = 0.3))
						    onsetQ2.5[[name]] = predict(loess(onsetQ2.5[[name]]  ~ onsetQ2.5$x, span = 0.3))
						    onsetQ97.5[[name]] = predict(loess(onsetQ97.5[[name]] ~ onsetQ97.5$x, span = 0.3))
      }
	}
	else if(smooth == "SPLINE") {
		for(name in stageNames) {
						    onsetMeans[[name]]  = predict(lm(onsetMeans[[name]]  ~ bs(onsetMeans$x, df = 5)))
						    onsetQ2.5[[name]] = predict(lm(onsetQ2.5[[name]] ~ bs(onsetQ2.5$x, df = 5)))
						    onsetQ97.5[[name]] = predict(lm(onsetQ97.5[[name]] ~ bs(onsetQ97.5$x, df = 5)))
      }
	}
	p = make_multistage_PPD_plot(onsetMeans, onsetQ2.5, onsetQ97.5, rainbow(nStages), responseData, covariateData[targetCovariateName], stageData, minResponse, maxResponse)
	return(p)
}


smooth_safe <- function(y, x) {
#print("smooth safe")
  if (sd(y) == 0) {
    return(rep(y[1], length(y)))
  } else {
    m <- mgcv::gam(y ~ s(x, bs = "cs", k = 10), method = "REML")
    return(predict(m))
  }
}

#' @importFrom dplyr tibble left_join %>%
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2  ggplot aes geom_ribbon geom_line geom_point coord_cartesian labs theme_minimal
#' @importFrom mgcv  gam 
make_multistage_PPD_plot = function(onsetMeans, onsetQ2.5, onsetQ97.5, stageCols, responseData, covariateData, stageData, minResponse, maxResponse) {
	library(ggplot2)

  warning("User-supplied stage colors (stageCols) are currently ignored and default RColorBrewer Set1 palette is used.")

	# Pivot to long format
	df_long <- onsetMeans %>%
		pivot_longer(-x, names_to = "stage", values_to = "mean") %>%
		left_join(
			  onsetQ2.5 %>%
				  pivot_longer(-x, names_to = "stage", values_to = "q2.5"),
			  by = c("x", "stage")
			  ) %>%
		left_join(
			  onsetQ97.5 %>%
				  pivot_longer(-x, names_to = "stage", values_to = "q97.5"),
			  by = c("x", "stage")
		)

    
#df_smooth <- df_long %>%
  #group_by(stage) %>%
  #arrange(x) %>%
  #mutate(
    #mean_s  = predict(gam(mean  ~ s(x, k = 10)), newdata = data.frame(x = x)),
    #q2.5_s  = predict(gam(q2.5 ~ s(x, k = 10)), newdata = data.frame(x = x)),
    #q97.5_s = predict(gam(q97.5 ~ s(x, k = 10)), newdata = data.frame(x = x))
  #)
#
		#print("Reorganized data for plotting")
		#print(df_long)

		df_data = data.frame(x=covariateData, observed=responseData, stage = stageData)
		names(df_data)[1] = "x"
		df_data$stage = paste0("stage", df_data$stage)

		#print("Raw data points")
		#print(df_data)

		# Plot
		p = ggplot(df_long, aes(x = x, y = mean, color = stage, fill = stage)) +
			geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.2, color = NA) +
    ## GAM smoothing of mean curve
      #geom_smooth(aes(y = mean), method = "gam", formula = y ~ s(x), se = FALSE, size = 1.2) +
			geom_line(linewidth = 1) +
			# Raw data points
			geom_point(data = df_data, aes(x = x, y = observed, color = stage), alpha = 0.7) +
			coord_cartesian(ylim = c(minResponse, maxResponse)) +
			labs(y = "Onset", x = "x") +
			scale_color_brewer(palette = "Set1") +
			scale_fill_brewer(palette = "Set1") +
			theme_minimal()
  #p = ggplot(df_smooth, aes(x = x, color = stage, fill = stage)) +
    #geom_ribbon(aes(ymin = q2.5_s, ymax = q97.5_s), alpha = 0.2, color = NA) +
    #geom_line(aes(y = mean_s), size = 1) +
			### Raw data points
		#geom_point( data = df_data, aes(x = x, y = observed, color = stage), alpha = 0.7) +
    #coord_cartesian(ylim = c(minResponse, maxResponse)) +
    #labs(y = "Onset", x = "x") +
    #scale_color_brewer(palette = "Set1") +
    #scale_fill_brewer(palette = "Set1") +
    #theme_minimal()

  return(p)
}
