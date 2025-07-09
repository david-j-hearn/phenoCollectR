
#' Title
#'
#' @param mu_O 
#' @param mu_D 
#' @param minResponse 
#' @param maxResponse 
#' @param type 
#'
#' @returns
#' @export
#'
#' @examples
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


#' Title
#'
#' @param N 
#' @param mu_O 
#' @param sigma_O 
#' @param mu_D 
#' @param sigma_D 
#' @param minResponse 
#' @param maxResponse 
#' @param type 
#' @param threshApprox 
#'
#' @returns
#' @export
#'
#' @examples
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


#' Title
#'
#' @param mu_O 
#' @param sigma_O 
#' @param mu_D 
#' @param sigma_D 
#' @param minResponse 
#' @param maxResponse 
#' @param type 
#'
#' @returns
#' @export
#'
#' @examples
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

#' Title
#'
#' @param mu_O 
#' @param sigma_O 
#' @param minResponse 
#' @param maxResponse 
#' @param type 
#'
#' @returns
#' @export
#'
#' @examples
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


#' Title
#'
#' @param N 
#' @param mu_O 
#' @param sigma_O 
#' @param minResponse 
#' @param maxResponse 
#' @param type 
#' @param threshApprox 
#'
#' @returns
#' @export
#'
#' @examples
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


#' Title
#'
#' @param mu_O 
#' @param mu_D 
#' @param minResponse 
#' @param maxResponse 
#' @param type 
#'
#' @returns
#' @export
#'
#' @examples
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


#' Title
#'
#' @param mu_O 
#' @param sigma_O 
#' @param mu_D 
#' @param sigma_D 
#' @param minResponse 
#' @param maxResponse 
#' @param alpha 
#' @param type 
#'
#' @returns
#' @export
#'
#' @examples
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

#' Title
#'
#' @param N 
#' @param mu_O 
#' @param sigma_O 
#' @param mu_D 
#' @param sigma_D 
#' @param minResponse 
#' @param maxResponse 
#' @param alpha 
#' @param type 
#'
#' @returns
#' @export
#'
#' @examples
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

#' Title
#'
#' @param mu_O 
#' @param sigma_O 
#' @param mu_D 
#' @param sigma_D 
#' @param minResponse 
#' @param maxResponse 
#' @param alpha 
#' @param type 
#'
#' @returns
#' @export
#'
#' @examples
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

#' Title
#'
#' @param mu_O 
#' @param sigma_O 
#' @param minResponse 
#' @param maxResponse 
#' @param alpha 
#' @param type 
#'
#' @returns
#' @export
#'
#' @examples
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

#' Title
#'
#' @param N 
#' @param mu_O 
#' @param sigma_O 
#' @param minResponse 
#' @param maxResponse 
#' @param alpha 
#' @param type 
#'
#' @returns
#' @export
#'
#' @examples
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

#' Title
#'
#' @param mu_O 
#' @param sigma_O 
#' @param mu_D 
#' @param sigma_D 
#' @param minResponse 
#' @param maxResponse 
#' @param alpha 
#' @param type 
#'
#' @returns
#' @export
#'
#' @examples
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

#' Title
#'
#' @param mu_O 
#' @param sigma_O 
#' @param mu_D 
#' @param sigma_D 
#' @param minResponse 
#' @param maxResponse 
#' @param type 
#'
#' @returns
#' @export
#'
#' @examples
SD.C =  function(mu_O=NA, sigma_O, mu_D=NA, sigma_D=NA, minResponse=0, maxResponse=365, type=c("GP","BB")) {
	type = match.arg(type)
	if(type=="BB") {
		if(is.na(sigma_D) || is.na(mu_O) || is.na(mu_C)) {
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

#' Title
#'
#' @param N 
#' @param mu_O 
#' @param sigma_O 
#' @param mu_D 
#' @param sigma_D 
#' @param minResponse 
#' @param maxResponse 
#' @param type 
#' @param intFailLow 
#' @param intFailHigh 
#'
#' @returns
#' @export
#'
#' @examples
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

#' Title
#'
#' @param mu_O 
#' @param sigma_O 
#' @param mu_D 
#' @param sigma_D 
#' @param minResponse 
#' @param maxResponse 
#' @param type 
#'
#' @returns
#' @export
#'
#' @examples
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

#' Title
#'
#' @param sigma_O 
#' @param type 
#'
#' @returns
#' @export
#'
#' @examples
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

#' Title
#'
#' @param N 
#' @param mu_O 
#' @param sigma_O 
#' @param minResponse 
#' @param maxResponse 
#' @param intFailLow 
#' @param intFailHigh 
#' @param type 
#'
#' @returns
#' @export
#'
#' @examples
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

#' Title
#'
#' @param mu_O 
#' @param sigma_O 
#' @param mu_D 
#' @param sigma_D 
#' @param minResponse 
#' @param maxResponse 
#' @param type 
#'
#' @returns
#' @export
#'
#' @examples
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

#' Title
#'
#' @param responseData 
#' @param minResponse 
#' @param maxResponse 
#' @param minS 
#' @param maxS 
#' @param init_params 
#' @param hyperparameters 
#' @param type 
#'
#' @returns
#' @export
#'
#' @examples
getMAP = function(responseData, minResponse=0, maxResponse=365,minS=1, maxS=3000,  init_params = c(180,20,60,7), hyperparameters = c(100, 7, 60, 6, 24, 12, 24, 12), type="BB") {
	if(type == "BB") {
		getMAP.T.BB(fileOrData=responseData, minResponse=minResponse, maxResponse=maxResponse, minS=minS, maxS=maxS,  init_params = init_params, hyperparameters = hyperparameters, dataProvided=TRUE)
	}
	else {
		stop(paste("Only model type beta onset, beta duration (BB) is implemented."))
	}
}

#' Title
#'
#' @param responseData 
#' @param minResponse 
#' @param maxResponse 
#' @param minS 
#' @param maxS 
#' @param init_params 
#' @param type 
#'
#' @returns
#' @export
#'
#' @examples
getMLE = function(responseData, minResponse=0, maxResponse=365, minS=1, maxS=3000, init_params = c(180, 20, 60, 7), type="BB") {
	if(type == "BB") {
		getMLE.T.BB(fileOrData=responseData, minResponse=minResponse, maxResponse=maxResponse, minS=minS, maxS=maxS, init_params = init_params, dataProvided=TRUE)
	}
	else {
		stop(paste("Only model type beta onset, beta duration (BB) is implemented."))
	}
}

#' Title
#'
#' @param mu_O 
#' @param sigma_O 
#' @param mu_D 
#' @param sigma_D 
#' @param minResponse 
#' @param maxResponse 
#' @param type 
#'
#' @returns
#' @export
#'
#' @examples
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

#' Title
#'
#' @param N 
#' @param mu_O 
#' @param sigma_O 
#' @param mu_D 
#' @param sigma_D 
#' @param minResponse 
#' @param maxResponse 
#' @param type 
#' @param precision 
#' @param includePlot 
#'
#' @returns
#' @export
#' @importFrom ForestFit fitWeibull
#'
#' @examples
fitWeibullExtremes = function(N, mu_O, sigma_O, mu_D, sigma_D=NA, minResponse=0, maxResponse=365, type=c("GP","BB"), precision=10000, includePlot=FALSE) {
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

		Ok1 = rOk1.BB(n=precision, N=N, mu_O=mu_O, sigma=sigma_O, minResponse=minResponse, maxResponse=maxResponse)
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

		hist(Ok1, breaks=x, col=rgb(1,1,0,0.5), xlim=c(xlim[1]-inc, xlim[2]+inc), probability=TRUE)
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

#' Title
#'
#' @param responseData 
#' @param onsetCovariateData 
#' @param durationCovariateData 
#'
#' @returns
#' @export
#'
#' @examples
runStandardLinearModel = function(responseData, onsetCovariateData, durationCovariateData) {
        # Combine response and predictors into one data frame
        merged = merge_df_by_column(onsetCovariateData, durationCovariateData)
        df <- data.frame(DOY = responseData, merged)

        # Dynamically construct formula
        formula <- as.formula(paste("DOY ~", paste(colnames(merged), collapse = " + ")))

        # Fit the model
        fit <- lm(formula, data = df)

        print(summary(fit))
		return(fit)
}

#' Title
#'
#' @param type 
#' @param responseData 
#' @param hyperparams_noCovariates 
#' @param onsetCovariateData 
#' @param durationCovariateData 
#' @param onsetHyperBeta 
#' @param onsetHyperAnchor 
#' @param durationHyperBeta 
#' @param durationHyperAnchor 
#' @param cessationHyperAnchor 
#' @param sigmaHyper 
#' @param minResponse 
#' @param maxResponse 
#' @param maxDiv 
#' @param setStringent 
#' @param runMAP 
#' @param processExtremes 
#' @param N 
#' @param partitionDataForPriors 
#' @param maximizeSampleSize 
#' @param threshApprox 
#' @param byPassChecks 
#' @param priorLevel
#' @param ... 
#'
#' @returns
#' @export
#'
#' @examples
runStanPhenology = function(type=c("intercept-only","full"), responseData, hyperparams_noCovariates=NA, onsetCovariateData=NA, durationCovariateData=NA, onsetHyperBeta=NA, onsetHyperAnchor=NA, durationHyperBeta=NA, durationHyperAnchor=NA, cessationHyperAnchor=NA, sigmaHyper=NA, minResponse=0, maxResponse=365, maxDiv=0, setStringent=TRUE, runMAP=FALSE, processExtremes=TRUE, N=500, partitionDataForPriors=TRUE, maximizeSampleSize=FALSE, threshApprox=NA, byPassChecks=FALSE,priorLevel=2, ...) {

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

	if(!(type=="intercept-only" || type=="full" ) ) {
		cat(paste("Unsupported type: ", type, "\nType should be 'intercept-only' or 'full'.\n"))
		stop("Unsupported type error.")
	}

	if(partitionDataForPriors) {
		cat("The data will be partitioned into two sets with 30% and 70% of the data. \n\n30% will be used to carry out a preliminary analysis using quantiles to estimate the prior distribution hyperparameter values. \n\n70% of the data will be used to carry out a Stan Bayesian analysis to obtain the posterior distributions of parameters.\n\nIf other hyperparameter information was provided as input, it will be ignored. \n\nThe calculated values based on quantiles are approximate; you may need to use other sources of data to get better estimates of prior hyperparameter values, especially if the Stan run results in divergences or other poor diagnostics.\n\n")

		prop = 0.3
		scale = (maxResponse-minResponse) / (365 - 0)

		if(type=="intercept-only") {
			if(partitionDataForPriors) {
			partition = partitionResponseData(responseData = responseData, prop = prop)
			responseData = partition$dataForInference
			hyperparams_noCovariates = getHyperparametersViaQuantiles(responseDataForPrior = partition$dataForPrior, scale = scale)
			}
		}
		else if(type=="full") {
			if(partitionDataForPriors) {
			#partition data
			partition = partitionResponseCovariateData(responseData=responseData, onsetCovariateData=onsetCovariateData, durationCovariateData=durationCovariateData, prop=prop)

			#get the data for inference
			#if(length(partition$responseDataForInference)<150 && maximizeSampleSize) {
			if(maximizeSampleSize) {
				   #don't reduce the amount of data - this is statistically invalid because the same data are used to estimate the prior hyperparameters and carry out the full Bayesian analysis, but when the sample size is too small, the error on the estimates will be unacceptably large unless a fuller data set is used.
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
			cessationHyperAnchor = prior$cessationHyperAnchor
			sigmaHyper = prior$sigmaHyper
		}
		}
	}

	if(type == "intercept-only") {
		cat("No covariates will be included in this analysis.\n")
		if(sum(is.na(hyperparams_noCovariates)) || length(hyperparams_noCovariates) != 6) stop("Expecting six hyperparameter values (mean and sd for mean onset, mean and sd for mean duration, mean and sd for sigma. Or, if you want hyperparameter values to be estimated for you, set 'partitionDataForPriors' to TRUE.")
		runStan.NoCovariates.T.GP(fileOrData=responseData, minResponse=minResponse, maxResponse=maxResponse, hyperparameters = hyperparams_noCovariates, dataProvided=TRUE, runMAP=runMAP, setStringent=setStringent, maxDiv=maxDiv, processExtremes=processExtremes, N=N, threshApprox=threshApprox, ...)
	}
	else if(type == "full" ) {
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
			if(!is.vector(cessationHyperAnchor)) {
				cat("Cessation anchor mean and standard deviation hyperparameter values are needed for the full model and should be provided as a two-element vector.\n")
				stop("Please provide appropriate inputs")
			}
			cat("Calling specialized runStan functions.\n")
		runStan.WithCovariates.T.GP(response=responseData, minResponse=minResponse, maxResponse=maxResponse, onsetCovariateData=onsetCovariateData, durationCovariateData=durationCovariateData, onsetHyperBeta=onsetHyperBeta, onsetHyperAnchor=onsetHyperAnchor, durationHyperBeta=durationHyperBeta, durationHyperAnchor=durationHyperAnchor, cessationHyperAnchor=cessationHyperAnchor, sigmaHyper=sigmaHyper, setStringent=setStringent, dataProvided=TRUE, priorLevel=priorLevel)
		}
	}
}

