parameter_checks = function(mu_O=NA, sigma_O=NA, mu_D=NA, sigma_D=NA, mu_C=NA, sigma_C=NA, N=NA, n=NA, minResponse=NA, maxResponse=NA) {
	flag = FALSE
	if(is.na(minResponse) || is.na(maxResponse)) {
		flag = TRUE
	}

	if(!flag && minResponse<0) {
		stop("The minimum response value should be non-negative. Please shift your data so that all response values are positive.")
	}
	if(!flag && minResponse>=maxResponse) {
		stop("The minimum response value must be less than the maximum response value.")
	}
	if(!is.na(n)) {
		if(n <= 0) {
			stop("n must be larger than 0.")
		}
	}
	if(!is.na(N)) {
		if(N <= 0) {
			stop("N must be larger than 0.")
		}
	}
	if(!is.na(mu_C)) {
		if(!flag && (mu_C <= minResponse || mu_C >= maxResponse)) {
			stop(paste0("mu_C must be between ", minResponse, " and ", maxResponse))
		}
	}
	if(!is.na(mu_C) && !is.na(mu_O)) {
		if(mu_C <= mu_O) {
			stop("mu_C must be larger than mu_O")
		}
	}
	if(!is.na(mu_O) ) {
		if(mu_O<=0) {
			stop("mu_O must be larger than 0.")
		}
		if(!flag && (mu_O <= minResponse || mu_O >= maxResponse)) {
			stop(paste0("mu_O must be between ", minResponse, " and ", maxResponse))
		}
	}
	if(!is.na(mu_D) ) {
		if(mu_D<=0) {
			stop("mu_D must be larger than 0.")
		}
		if(!flag && (mu_D <= minResponse || mu_D >= maxResponse)) {
			stop(paste0("mu_D must be between ", minResponse, " and ", maxResponse))
		}
	}
	if(!is.na(sigma_O)) {
		if(sigma_O<=0) {
			stop("sigma_O must be larger than 0.")
		}
		if(!flag && (sigma_O > maxResponse-minResponse)) {
			stop("sigma_O is too large for the range of possible values.")
		}
	}
	if(!is.na(sigma_D)) {
		if(sigma_D<=0) {
			stop("sigma_D must be larger than 0.")
		}
		if(!flag && (sigma_D > maxResponse-minResponse)) {
			stop("sigma_D is too large for the range of possible values.")
		}
	}
	if(!is.na(sigma_C)) {
		if(sigma_C<=0) {
			stop("sigma_C must be larger than 0.")
		}
		if(!flag && (sigma_C > maxResponse-minResponse)) {
			stop("sigma_C is too large for the range of possible values.")
		}
	}
}

merge_df_by_column = function(df1, df2) {
	# Identify common and unique columns
	common_cols  =  intersect(names(df1), names(df2))
	unique_df2  =  setdiff(names(df2), common_cols)

	for (col in common_cols) {
		if (!all(df1[[col]] == df2[[col]])) {
			stop(paste("Column", col, "differs between data frames"))
		}
	}

	# Combine: df1 + only the unique columns from df2
	merged_df  =  cbind(df1, df2[unique_df2])
	return(merged_df)
}

generate_summary_table  =  function(data, columns, output_file = "summary_table.csv", taxonName, originalFile, originalN, rawN) {
	# Check if all specified columns exist in the data frame
	if (!all(columns %in% colnames(data))) {
		stop("One or more specified columns are not present in the data frame.")
	}

	# Subset the data frame to the selected columns
	selected_data  =  data[ , columns, drop = FALSE]

	# Compute summary statistics for each column
	summary_stats  =  data.frame(
								Variable = columns,
								N = sapply(selected_data, function(x) sum(!is.na(x))),
								Min = sapply(selected_data, min, na.rm = TRUE),
								Max = sapply(selected_data, max, na.rm = TRUE),
								Mean = sapply(selected_data, mean, na.rm = TRUE),
								SD = sapply(selected_data, sd, na.rm = TRUE),
								stringsAsFactors = FALSE
	)

	# Add additional descriptive data
	attr(summary_stats, "TaxonName")  =  taxonName
	attr(summary_stats, "OriginalDataFile")  =  originalFile

	total_row  =  data.frame(
							Variable = "PostCleanN",
							N = nrow(data),
							Min = NA, Max = NA, Mean = NA, SD = NA,
							stringsAsFactors = FALSE
	)
	summary_stats  =  rbind(summary_stats, total_row)

	originalN  =  data.frame(
							Variable = "ReproductiveN",
							N = originalN,
							Min = NA, Max = NA, Mean = NA, SD = NA,
							stringsAsFactors = FALSE
	)
	summary_stats  =  rbind(summary_stats, originalN)

	taxonName  =  data.frame(
							Variable = "TaxonName",
							N = taxonName,
							Min = NA, Max = NA, Mean = NA, SD = NA,
							stringsAsFactors = FALSE
	)
	summary_stats  =  rbind(summary_stats, taxonName)

	originalFile  =  data.frame(
							   Variable = "DataFile",
							   N = originalFile,
							   Min = NA, Max = NA, Mean = NA, SD = NA,
							   stringsAsFactors = FALSE
	)
	summary_stats  =  rbind(summary_stats, originalFile)

	rawN  =  data.frame(
					   Variable = "AllSpecimenN",
					   N = rawN,
					   Min = NA, Max = NA, Mean = NA, SD = NA,
					   stringsAsFactors = FALSE
	)
	summary_stats  =  rbind(summary_stats, rawN)


	# Write summary table to CSV file
	write.csv(summary_stats, file = output_file, row.names = FALSE)

	return(summary_stats)
}

remove_outliers_lm = function(data, response, predictors,
#mahalanobis: Outlier Analysis 2017 Aggarwal (section 2.3.4)
#residuals: Generalized Linear Model Diagnostics Using the Deviance and Single Case Deletions
#Cook 1977
							  mahal_p = 0.975,
							  residual_cutoff = 3,
							  cooks_cutoff = NULL) {
	# Extract model frame
	formula = as.formula(paste(response, "~", paste(predictors, collapse = "+")))
	model_data = model.frame(formula, data)

	y = model_data[[response]]
	X = model.matrix(formula, data = model_data)[, -1, drop = FALSE] # drop intercept

	n = nrow(X)
	p = ncol(X)

	## 1. Multivariate outliers in predictors (Mahalanobis)
	md = mahalanobis(X, colMeans(X), cov(X))
	md_cutoff = qchisq(mahal_p, df = p)
	keep_md = md <= md_cutoff

	## Filter data before linear modeling
	data_md_clean = model_data[keep_md, ]
	lm_clean = lm(formula, data = data_md_clean)

	## 2. Outliers in response: studentized residuals and Cook’s distance
	rstud = rstudent(lm_clean)
	cooks = cooks.distance(lm_clean)

	if (is.null(cooks_cutoff)) {
		cooks_cutoff = 4 / nrow(data_md_clean)
	}

	keep_resid = abs(rstud) <= residual_cutoff
	keep_cooks = cooks <= cooks_cutoff
	keep_response = keep_resid & keep_cooks

	## Final filtered data
	cleaned_data = data_md_clean[keep_response, ]

	## Final model
	final_model = lm(formula, data = cleaned_data)

	return(list(
				cleaned_data = cleaned_data,
				final_model = final_model,
				removed_rows = setdiff(rownames(model_data), rownames(cleaned_data)),
				diagnostics = list(
								   mahalanobis = md,
								   mahalanobis_cutoff = md_cutoff,
								   studentized_residuals = rstud,
								   cooks_distance = cooks,
								   cooks_cutoff = cooks_cutoff
				)
				))
}

scaled_dbeta = function(y, shape1, shape2, minResponse = 0, maxResponse = 1) {
	x = (y - minResponse) / (maxResponse - minResponse)
	dbeta(x, shape1, shape2) / (maxResponse - minResponse)
}

#' @export
preparePhenologyData = function(dataFile, responseVariableName, onsetCovariateNames, durationCovariateNames, removeDuplicateRows=TRUE, removeOutliers=FALSE, removeIncomplete=TRUE, dataSummaryDirectory=NA, taxonName=NA, origN=NA) {


	if(!file.exists(dataFile)) {
		stop(paste("Your phenology data file, ", dataFile, " could not be found. Please check the file name, path, and current working directory."))
	}

	data = read.table(dataFile, header=TRUE, sep='\t')
	dataOrig = data

	if(removeIncomplete) {
		data = data[complete.cases(data), ]
	}

	allCovariates = Reduce(union, list(responseVariableName,onsetCovariateNames,durationCovariateNames))

	if (!all(allCovariates %in% names(data))) {
		stop("Some variable names are not present in the data file. Check variable names, or alter them by removing spaces and special characters in the names.")
	}

	data = data[, allCovariates]

	if(removeDuplicateRows) {
		data = removeDuplicateRows(data)
	}

	if(removeOutliers) {
		noOutliers = remove_outliers_lm(data=data, response=responseVariableName, predictors = allCovariates)
		data = noOutliers$cleaned_data
	}

	if(!is.na(dataSummaryDirectory) && !is.na(taxonName) && !is.na(origN)) {
		if (!dir.exists(dataSummaryDirectory)) {
			dir.create(dataSummaryDirectory, recursive = TRUE, showWarnings = FALSE)
		}

		output_file = paste0(dataSummaryDirectory, "/dataSummary.", taxonName, ".txt")
		summaryTable = generate_summary_table(data=data, columns=allCovariates, output_file = output_file, taxonName=taxonName, originalFile=dataFile, originalN=nrow(dataOrig), rawN=origN)
	}

	responseData = data[[responseVariableName]]
	onsetCovariates = data[, onsetCovariateNames, drop = FALSE]
	durationCovariates = data[, durationCovariateNames, drop = FALSE]

	out = list(
			   originalData = dataOrig,
			   responseData = responseData,
			   onsetCovariateData = onsetCovariates,
			   durationCovariateData = durationCovariates
	)
	return(out)
}

#' @export
#' @importFrom quantreg rq
getHyperparametersViaQuantileRegression = function(responseDataForPrior,
												   onsetCovariateDataForPrior,
												   durationCovariateDataForPrior,
												   lowerQuantile = 0.1,
												   upperQuantile = 0.9,
												   confidence = 2,
												   scale = 1
												   ) {

	if (!identical(names(onsetCovariateDataForPrior), names(durationCovariateDataForPrior))) {
		stop("To set hyperparameter values automatically, the covariates for onset must be the same as the covariates for duration.")
	}

	X_onset = as.data.frame(onsetCovariateDataForPrior)
	df_onset = cbind(y = responseDataForPrior, X_onset)

	formula_onset = as.formula(paste("y ~", paste(names(X_onset), collapse = " + ")))
	fitO = rq(formula_onset, data = df_onset, tau = lowerQuantile)
	coeffsO = coef(fitO)
	if (any(is.na(coeffsO))) stop("NA in onset quantile regression coefficients.")
	alphaO = coeffsO["(Intercept)"]
	betaO = coeffsO[names(coeffsO) != "(Intercept)"]

	summary_fitO = summary(fitO, se = "boot")
	SEO = summary_fitO$coefficients[, "Std. Error"]
	SEO_b = SEO[names(SEO) != "(Intercept)"]
	SEO_a = SEO["(Intercept)"]

	x_means_onset = colMeans(X_onset)

	H_M_Anchor_meanOnset = alphaO + sum(betaO * x_means_onset)
	H_SD_Anchor_meanOnset = 14
	onsetHyperAnchor = c(mean = H_M_Anchor_meanOnset, sd = H_SD_Anchor_meanOnset)

	onsetHyperBeta = data.frame(mean = betaO, sd = confidence * SEO_b, row.names = names(betaO))

	## --- Duration Regression ---
	X_duration = as.data.frame(durationCovariateDataForPrior)
	df_duration = cbind(y = responseDataForPrior, X_duration)
	formula_duration = as.formula(paste("y ~", paste(names(X_duration), collapse = " + ")))
	fitC = rq(formula_duration, data = df_duration, tau = upperQuantile)
	coeffsC = coef(fitC)
	if (any(is.na(coeffsC))) stop("NA in duration quantile regression coefficients.")
	alphaC = coeffsC["(Intercept)"]
	betaC = coeffsC[names(coeffsC) != "(Intercept)"]

	summary_fitC = summary(fitC, se = "boot")
	SEC = summary_fitC$coefficients[, "Std. Error"]
	SEC_b = SEC[names(SEC) != "(Intercept)"]
	SEC_a = SEC["(Intercept)"]

	x_means_duration = colMeans(X_duration)
	H_M_Anchor_meanDuration = alphaC + sum(betaC * x_means_duration) - H_M_Anchor_meanOnset
	H_SD_Anchor_meanDuration = 14*scale
	durationHyperAnchor = c(mean = H_M_Anchor_meanDuration, sd = H_SD_Anchor_meanDuration)

	H_M_Anchor_meanCessation = H_M_Anchor_meanDuration + H_M_Anchor_meanOnset
	H_SD_Anchor_meanCessation = 14*scale
	cessationHyperAnchor = c(mean = H_M_Anchor_meanCessation, sd = H_SD_Anchor_meanCessation)

	durationHyperBeta = data.frame(mean = betaC - betaO, sd = confidence * (SEC_b + SEO_b), row.names = names(betaO))

	sigmaHyper = c(7*scale,3.5*scale)

	out = list(
			   onsetHyperAnchor = onsetHyperAnchor,
			   onsetHyperBeta = onsetHyperBeta,
			   durationHyperAnchor = durationHyperAnchor,
			   durationHyperBeta = durationHyperBeta,
			   cessationHyperAnchor = cessationHyperAnchor,
			   sigmaHyper = sigmaHyper
	)

	return(out)
}

#' @export
getHyperparametersViaQuantiles = function(responseDataForPrior, scale=1, lowerQuantile=0.1, upperQuantile=0.9) {

	if (lowerQuantile >= upperQuantile) stop("Lower quantile must be less than upper quantile.")

	#need to get hyperparams_noCovariates=NA, a vector of 6 values:
	#  mean of mean onset taken at the 10% quantile (q10%)
	H_M_MO = quantile(responseDataForPrior, probs = c(lowerQuantile))
	#  sd of mean onset estimated based on prior experience
	H_SD_MO = 14 * scale  #uncertainty of +/- two weeks, scaled to match user-input response range
	#  mean of mean duration inferred as the difference between the q90% and q10%
	H_M_MD = quantile(responseDataForPrior, probs = c(upperQuantile)) - H_M_MO
	#  sd of mean duration estimated based on prior experience
	H_SD_MD = 14 * scale  #uncertainty of +/- two weeks, scaled to match user-input response range
	#  mean of sigma estimated based on prior experience
	H_M_S = 7 * scale#plants start flowering within +/- a week of each other in their population, scaled to match user-input response range
	#  sd fo sigma estimated based on prior experience, scaled to match user-input response range
	H_SD_S = 3.5 * scale

	return(c(H_M_MO, H_SD_MO, H_M_MD, H_SD_MD, H_M_S, H_SD_S))
}

#' @export
partitionResponseCovariateData = function(responseData, onsetCovariateData, durationCovariateData, prop=0.3) {

	#perhaps set seed...

	nRes = length(responseData)

	if(nrow(onsetCovariateData) != nRes ) {
		stop("The number of onset covariate data points does not match the number of response data points.")
	}
	if(nrow(durationCovariateData) != nRes ) {
		stop("The number of duration covariate data points does not match the number of response data points.")
	}

	train_indices = sample(nRes, size = floor(prop * nRes), replace=FALSE)

	responseDataForInference = responseData[-train_indices]
	onsetCovariateDataForInference = onsetCovariateData[-train_indices,,drop=FALSE]
	durationCovariateDataForInference = durationCovariateData[-train_indices,,drop=FALSE]

	if(!is.data.frame(onsetCovariateData)) {
		stop("Partitioning didn't START WITH a data frame")
	}

	responseDataForPrior = responseData[train_indices]
	onsetCovariateDataForPrior = onsetCovariateData[train_indices,,drop=FALSE]
	durationCovariateDataForPrior = durationCovariateData[train_indices,,drop=FALSE]

	if(!is.data.frame(onsetCovariateDataForPrior)) {
		print(onsetCovariateDataForPrior)
		if(is.data.frame(onsetCovariateData)) { print("started with a data frame in onsetCovariateData") }
		stop("Partitioning didn't create a data frame. ug")
	}

	out = list(
			   responseDataForInference = responseDataForInference,
			   onsetCovariateDataForInference = onsetCovariateDataForInference,
			   durationCovariateDataForInference = durationCovariateDataForInference,

			   responseDataForPrior = responseDataForPrior,
			   onsetCovariateDataForPrior = onsetCovariateDataForPrior,
			   durationCovariateDataForPrior = durationCovariateDataForPrior
	)


	return(out)
}

#' @export
partitionResponseData = function(responseData, prop=0.3) {
	nRes = length(responseData)
	train_indices = sample(nRes, size = floor(prop * nRes), replace=FALSE)

	dataForPrior = responseData[train_indices]
	dataForInference = responseData[-train_indices]

	out = list(
			   dataForPrior = dataForPrior,
			   dataForInference = dataForInference
	)
	return(out)
}

removeDuplicateRows = function(df) {
	df_unique = df[!duplicated(df), ]
	return(df_unique)
}

#' @export
E.Kth.approx = function(N, mu, sigma, k) {
	if(k<0 || k>N || N<0) {
		stop("k must be between 1 and N inclusive, and N must be a positive integer.")
	}
	approx = mu + sigma * qnorm((k - pi/8)/(N - pi/4 + 1))
	return(approx)
}


#Obtains the unscaled hyperparameters from a file associated with a provided file name
getHyperparameters = function(hyperparameterFile) {
	hyperparameters = read.table(hyperparameterFile, header=TRUE, sep='\t')
	return(hyperparameters)
}

#Obtains min max scaled covariate data from a file associated with a provided file name
getCovariates = function(covariatesFile) {
	covariates = read.table(covariatesFile, header=TRUE, sep='\t')
	return(processCovariates(covariates))
}

processCovariates = function(covariates) {
	mins = apply(covariates,2,min)
	maxs = apply(covariates,2,max)
	means = apply(covariates,2,mean)
	scaledMeans = (means - mins) / (maxs - mins)
	SDs = apply(covariates,2,sd)
	scaledSDs = SDs / (maxs - mins)

	if(any(maxs-mins==0)) {
		stop("Constant covariate data is not allowed.")
	}

	scaledCovariates = as.data.frame(lapply(names(covariates), function(col_name) {
												(covariates[[col_name]] - mins[col_name]) / (maxs[col_name] - mins[col_name])
				}))

	colnames(scaledCovariates) = names(covariates)


	output = list(
				  mins = mins,
				  maxs = maxs,
				  means = means,
				  scaledMeans = scaledMeans,
				  SDs = SDs,
				  scaledSD = scaledSDs,
				  covariates = covariates,
				  scaledCovariates = scaledCovariates,
				  K = ncol(covariates)
	)

	return(output)
}

make_quantile_function_from_pdf = function(ddist, support = c(-10, 10), n = 1000) {
	# Grid of x values
	x = seq(support[1], support[2], length.out = n)
	dx = diff(x)[1]

	# Evaluate PDF
	pdf_vals = ddist(x)

	# Normalize the PDF to ensure it integrates to 1
	pdf_vals = pdf_vals / sum(pdf_vals * dx)

	# Compute the CDF
	cdf_vals = cumsum(pdf_vals) * dx

	# Construct the quantile function (inverse CDF)
	qfun = approxfun(cdf_vals, x, method = "linear", ties = "ordered")

	# Return a named function that accepts probability p
	quantile_func = function(p) {
		if (any(p < 0 | p > 1)) stop("Probabilities must be between 0 and 1.")
		qfun(p)
	}

	return(quantile_func)
}

isBetaFeasible = function(alpha, beta, minS=1, maxS=3000) {
	if(alpha<minS || alpha>maxS || beta<minS || beta>maxS) {
		return(FALSE)
	}
	return(TRUE)
}

getObservations = function(responseFile,minResponse=0,maxResponse=365) {
	observed = read.table(responseFile, header=TRUE)
	if(ncol(observed)>1) {
		print("The response file should contain only one column with the time of the observed response for each specimen. Each row represents a different specimen. There should be a header in the file.")
		return();
	}

	#scale data between max and min
	observed = observed[[1]]
	observed = (observed - minResponse) / (maxResponse - minResponse)
	return(observed)
}

beta_mean = function(alpha, beta) {
	if(any(alpha<0) || any(beta<0)) {
		stop("alpha and beta shape parameters for the beta distribution must be positive.")
	}
	return(alpha/(alpha+beta))
}

beta_sd = function(alpha, beta) {
	if(any(alpha<0) || any(beta<0)) {
		stop("alpha and beta shape parameters for the beta distribution must be positive.")
	}
	var = alpha*beta / ((alpha+beta)*(alpha+beta)*(alpha+beta+1))
	sd = sqrt(var)
	return(sd)
}

beta_alpha = function(mean, sd) {
	if (any(mean <= 0) || any(mean >= 1)) stop("Beta mean must be in (0, 1).")
	if(any(sd <=0)) stop("Standard deviations must be positive.")
	((1 - mean) / sd^2 - 1 / mean) * mean^2
}

beta_beta = function(mean, sd) {
	if (any(mean <= 0) || any(mean >= 1)) stop("Beta mean must be in (0, 1).")
	if(any(sd <=0)) stop("Standard deviations must be positive.")
	((1 - mean) / sd^2 - 1 / mean) * mean * (1 - mean)
}

# --- 1. Safe beta helpers ---
dbeta_safe = function(x, shape1, shape2) {
	ifelse(x <= 0 | x >= 1, 0, dbeta(x, shape1, shape2))
}

pbeta_safe = function(q, shape1, shape2) {
	ifelse(q <= 0, 0,
		   ifelse(q >= 1, 1, pbeta(q, shape1, shape2)))
}

#' @export
summaryStanDiagnostics = function(stanRunResult, NTotal, taxonName, diagnostics=c("num_divergent", "num_max_treedepth", "ebfmi")) {

	print(NTotal)
	print(taxonName)
	print(stanRunResult$data$N)

	summary = stanRunResult$result$diagnostic_summary() %>% as.data.frame()

	n = nrow(summary)
	c = ncol(summary)

	summary = data.frame(
						 taxon.name = taxonName,
						 n.total = NTotal,
						 n.inference = stanRunResult$data$N,
						 as.data.frame(t(colMeans(summary)))
	)

	return(summary)
}

#' @export
#' @importFrom posterior as_draws_df quantile2
summarizePhenologyResults = function(stanRunResult, taxonName, measures=c("mean", "median", "sd", "mad"), quantiles=c(2.5, 97.5), convergence=c("rhat","ess_bulk", "ess_tail"), standardLinearModel=NA, processExtremes=TRUE, N=500) {
	if(processExtremes) {
		variables=c("anchor_O", "anchor_D", "anchor_C", "anchor_T", "anchorOk1", "anchorCkN", "alpha_O", "alpha_D", "alpha_C", "alpha_T", "alphaOk1", "alphaCkN", "beta_O", "beta_D", "beta_C", "beta_T", "betaOk1", "betaCkN", "sigma")
	}
	else {
		variables=c("anchor_O", "anchor_D", "anchor_C", "anchor_T", "alpha_O", "alpha_D", "alpha_C", "alpha_T", "beta_O", "beta_D", "beta_C", "beta_T", "sigma")
	}

	probs = quantiles/100

	onsetCovariateData = stanRunResult$onsetCovariateData
	durationCovariateData = stanRunResult$durationCovariateData

	onsetCovariateNames = colnames(onsetCovariateData)
	durationCovariateNames = colnames(durationCovariateData)
	unionCovariateNames = union(onsetCovariateNames,durationCovariateNames)

	idx_O  =  match(onsetCovariateNames, unionCovariateNames)
	idx_D  =  match(durationCovariateNames, unionCovariateNames)

	K_O = length(onsetCovariateNames)
	K_D = length(durationCovariateNames)
	K_union = length(unionCovariateNames)

	summary = stanRunResult$result$summary(variables = variables, measures, quantiles = ~ quantile2(., probs=probs), convergence) %>% as.data.frame()

	posterior  =  stanRunResult$result$draws()
	posterior_df  =  as_draws_df(posterior)
	mean(posterior_df$beta < 0)

	taxon.name =rep(taxonName, nrow(summary))
	type = rep(NA, nrow(summary))
	covariate = rep(NA, nrow(summary))
	posterior.prob.neg.slope = rep(NA, nrow(summary))
	standard.linear.model = rep(NA, nrow(summary))
	in.95CI = rep(NA, nrow(summary))

	if(class(standardLinearModel)=="lm") {
		coefs  =  coef(standardLinearModel)
	}

	summary = data.frame(
						 taxon.name = taxon.name,
						 summary[1],
						 type = type,
						 covariate = covariate,
						 posterior.prob.neg.slope = posterior.prob.neg.slope,
						 standard.linear.model = standard.linear.model,
						 in.95CI = in.95CI,
						 summary[2:ncol(summary)]
	)


	for(i in 1:K_O) {
		varO = paste0("beta_O[", i, "]")
		summary[summary$variable==varO, "type"] = "onset (O)"
		summary[summary$variable==varO, "covariate"] = onsetCovariateNames[i]
		summary[summary$variable == varO, "posterior.prob.neg.slope"] = mean(posterior_df[[varO]] < 0)
		summary[summary$variable==varO, "variable"] = "slope"
		if(processExtremes) {
			varOk1 = paste0("betaOk1[", i, "]")
			summary[summary$variable==varOk1, "type"] = paste0("first onset (Ok1): N=", N)
			summary[summary$variable==varOk1, "covariate"] = onsetCovariateNames[i]
			summary[summary$variable == varOk1, "posterior.prob.neg.slope"] = mean(posterior_df[[varOk1]] < 0)
			summary[summary$variable==varOk1, "variable"] = "slope"
		}
	}
	for(i in 1:K_D) {
		varD = paste0("beta_D[", i, "]")
		summary[summary$variable==varD, "type"] = "duration (D)"
		summary[summary$variable==varD, "covariate"] = durationCovariateNames[i]
		summary[summary$variable == varD, "posterior.prob.neg.slope"] = mean(posterior_df[[varD]] < 0)
		summary[summary$variable==varD, "variable"] = "slope"
	}
	for(i in 1:K_union) {
		varC = paste0("beta_C[", i, "]")
		varT = paste0("beta_T[", i, "]")
		summary[summary$variable==varC, "type"] = "cessation (C)"
		summary[summary$variable==varC, "covariate"] = unionCovariateNames[i]
		summary[summary$variable == varC, "posterior.prob.neg.slope"] = mean(posterior_df[[varC]] < 0)
		summary[summary$variable==varT, "type"] = "observed times (T)"
		summary[summary$variable==varT, "covariate"] = unionCovariateNames[i]
		summary[summary$variable == varT, "posterior.prob.neg.slope"] = mean(posterior_df[[varT]] < 0)
		if(class(standardLinearModel)=="lm") {
			summary[summary$variable == varT, "standard.linear.model"] = coefs[[unionCovariateNames[i]]]
			summary[summary$variable == varT, "in.95CI"] = (summary[summary$variable==varT,]$q2.5<coefs[[unionCovariateNames[i]]] & summary[summary$variable==varT,]$q97.5>coefs[[unionCovariateNames[i]]])
		}

		summary[summary$variable==varC, "variable"] = "slope"
		summary[summary$variable==varT, "variable"] = "slope"

		if(processExtremes) {
			varCkN = paste0("betaCkN[", i, "]")
			summary[summary$variable==varCkN, "type"] = paste0("last cessation (CkN): N=", N)
			summary[summary$variable==varCkN, "covariate"] = unionCovariateNames[i]
			summary[summary$variable == varCkN, "posterior.prob.neg.slope"] = mean(posterior_df[[varCkN]] < 0)
			summary[summary$variable==varCkN, "variable"] = "slope"
		}
	}

	summary[summary$variable=="anchor_O", "type"] = "onset (O)"
	summary[summary$variable=="anchor_D", "type"] = "duration (D)"
	summary[summary$variable=="anchor_C", "type"] = "cessation (C)"
	summary[summary$variable=="anchor_T", "type"] = "observed (T)"

	summary[summary$variable=="anchor_O", "variable"] = "anchor"
	summary[summary$variable=="anchor_D", "variable"] = "anchor"
	summary[summary$variable=="anchor_C", "variable"] = "anchor"
	summary[summary$variable=="anchor_T", "variable"] = "anchor"

	if(processExtremes) {
		summary[summary$variable=="anchorOk1", "type"] = paste0("first onset (Ok1): N=", N)
		summary[summary$variable=="anchorCkN", "type"] = paste0("last cessation (CkN): N=", N)
		summary[summary$variable=="anchorOk1", "variable"] = "anchor"
		summary[summary$variable=="anchorCkN", "variable"] = "anchor"
	}

	summary[summary$variable=="alpha_O", "type"] = "onset (O)"
	summary[summary$variable=="alpha_D", "type"] = "duration (D)"
	summary[summary$variable=="alpha_C", "type"] = "cessation (C)"
	summary[summary$variable=="alpha_T", "type"] = "observed (T)"
	if(class(standardLinearModel)=="lm") {
		summary[summary$variable == "alpha_T", "standard.linear.model"] = coefs[1]
			summary[summary$variable == "alpha_T", "in.95CI"] = (summary[summary$variable=="alpha_T",]$q2.5<coefs[1] & summary[summary$variable=="alpha_T",]$q97.5>coefs[1])
	}

	summary[summary$variable=="alpha_D", "variable"] = "intercept"
	summary[summary$variable=="alpha_C", "variable"] = "intercept"
	summary[summary$variable=="alpha_O", "variable"] = "intercept"
	summary[summary$variable=="alpha_T", "variable"] = "intercept"

	if(processExtremes) {
		summary[summary$variable=="alphaOk1", "type"] = paste0("first onset (Ok1): N=", N)
		summary[summary$variable=="alphaCkN", "type"] = paste0("last cessation (CkN): N=", N)
		summary[summary$variable=="alphaOk1", "variable"] = "intercept"
		summary[summary$variable=="alphaCkN", "variable"] = "intercept"
	}

	return(summary)
}

