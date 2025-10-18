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

#' @importFrom utils write.csv
#' @importFrom stats sd
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

#' @importFrom stats model.frame mahalanobis cov qchisq lm rstudent cooks.distance model.matrix dist
remove_outliers_lm = function(data, response, predictors,
#mahalanobis: Outlier Analysis 2017 Aggarwal (section 2.3.4)
#residuals: Generalized Linear Model Diagnostics Using the Deviance and Single Case Deletions
#Cook 1977
							  mahal_p = 0.975,
							  residual_cutoff = 3,
							  cooks_cutoff = NULL) {
	# Extract model frame
	# formula = as.formula(paste(response, "~", paste(predictors, collapse = "+")))
  ## DANIEL: Fix warning with having the response among predictors:
  if( response %in% predictors ){
    predictors <- predictors[!predictors == response]
  }
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


#' Prepare data for analysis
#'
#' This function prepares the phenological data for Bayesian estimation analysis using "runStanPhenology". The main inputs are the path to the data file and variable names. Note that all variable names should match the names of the columns in file.
#' @details
#' The file whose name is provided as input must be a text-formatted, tab-delimited spreadsheet. Each row corresponds to a single specimen, and the columns provide the specimen data. Any number of columns is fine, but the input names of the response variable, the onset covariates, and the duration covariates must match the corresponding header names in the spreadsheet. 
#' @param dataFile A path to a file with the data. See 'Details' for the correct configuration of this file.
#' @param responseVariableName The name of the response variable. Must match the name of a column in 'dataFile'.
#' @param onsetCovariateNames A vector with the names of the onset covariates. Must match name(s) of column(s) in 'dataFile'.
#' @param durationCovariateNames A vector with the names of the duration covariates. Must match name(s) of column(s) in 'dataFile'.
#' @param removeDuplicateRows If duplicated rows should be removed (default: TRUE).
#' @param removeOutliers  If outliers should be removed (default: FALSE). Outliers will be removed using three methods: Mahalanobis distance from: Outlier Analysis 2017 Aggarwal (section 2.3.4); residuals beyond a threshold; Cooks distance (Cook, 1977). These will be run at: Mahalanobis p = 0.975, residual cutoff = 3, Cook's cutoff = 4 / number of items in dataset. If you prefer to set different criteria, you can call the internal function that removes outliers: phenoCollectR:::remove_outliers_lm, and set these values manually. This function provides multiple outputs. With the resulting cleaned dataset that is the output, you can call the current function, and set this removeOutliers to FALSE.
#' @param removeIncomplete If rows with missing data should be removed (default: TRUE).
#' @param dataSummaryDirectory Optionally provide the name of a directory where a summary of the data will be saved in a tab-delimited file named dataSummary.<taxonName>.txt
#' @param taxonName Optionally, the name of the taxon to be stored in the data summary file.
#' @param origN Optionally, the initial number of specimens before any filtering of specimens took place (Use only if creating optional data summary file). 
#'
#' @return Returns a list object with elements necessary to use the "runStanPhenology" function. Elements include the original data (originalData) and the cleaned response data (responseData), onset covariate data (onsetCovariateData), and duration covariate data (durationCovariateData). The output data remain unscaled. The runStanPhenology function will scale data appropraitely. 
#' @importFrom utils read.table
#' @importFrom stats complete.cases
#' @export
#'
#' @examples
#' \donttest{
#' ##get the file name with data for the blood root plant
#' file  =  getDatasetPath("Sanguinaria_canadensis")
#' ## See documentation for more species:
#' help(getDatasetPath)
#' ##define the covariate names - remove up to all but 1
#' vars = c("Latitude", "Year", "Elevation", "AnnualMonthlyAverageTemp", "SpringMonthlyAverageTemp"
#'           , "FirstQuarterMonthlyAverageTemp")
#' ##get the phenology data
#' data  =  preparePhenologyData(dataFile=file, responseVariableName="DOY", onsetCovariateNames=vars
#'                               , durationCovariateNames=vars, taxonName="Sanguinaria_canadensis"
#'                               , removeOutliers=TRUE)
#' ##run the Stan sampler
#' stanResult  =  runStanPhenology(type="full", responseData = data$responseData
#'                                 , onsetCovariateData = data$onsetCovariateData
#'                                 , durationCovariateData = data$durationCovariateData
#'                                 , partitionDataForPriors = TRUE)
#' ##summarize the Stan run
#' stanSummary  =  summarizePhenologyResults(stanRunResult = stanResult
#'                                           , taxonName = "Sanguinaria_canadensis"
#'                                           ,standardLinearModel = TRUE)
#' ##make posterior predictive graph
#' pp  =  makePosteriorPredictivePlot(stanResult = stanResult, responseData = data$responseData
#'                                    , targetCovariateName = "SpringMonthlyAverageTemp"
#'                                    , onsetCovariateData = data$onsetCovariateData
#'                                    , durationCovariateData = data$durationCovariateData)
#' ##display the posterior predictive graph		
#' print(pp)		
#' }
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


#' Hyperparameter estimates using quantile regression
#' 
#' This function runs a quantile regression analysis and provides the quantile estimates of slope coefficients for onset and cessation. The quantile regression approach for phenology was first implemented by Park et al. (2024) https://doi.org/10.1111/ecog.06961
#'
#' Slope estimates for coefficients tend to be accurate at the 10% (for onset) and 90% (for cessation) quantiles. However, estimates of the onset and duration values themselves are prone to systematic errors due to fixing aribrary quantile levels. 
#'
#' @param responseDataForPrior A vector of response data (e.g., day of year of collection values)
#' @param onsetCovariateDataForPrior  A data frame with the onset covariate data. This can be obtained from the preparePhenologyData function.
#' @param durationCovariateDataForPrior  A data frame with the duration covariate data. This can be obtained from the preparePhenologyData function.
#' @param lowerQuantile The quantile used to estimate the onset model
#' @param upperQuantile The quantile used to estimate the cessation model
#' @param confidence Set the number of standard deviations to use for the prior. (default: 2)
#' @param scale Set the scaling of hardcoded values for the standard deviation of the mean duration (two weeks), the standard deviation of the mean onset (two weeks), and the sigma mean (one week) and sigma standard deviation (half a week) values. For example, setting the scale to 0.5 would make the standard deviation of the prior distribution for mean onset equal to one week rather than two (default: 1)
#' @details 
#' This approach to define prior hyperparameters depends on arbitrary quantile cutoffs and, although it may work well in practice for simulated botanical collection data, it is not currently theoretically justified, nor is it known how well it works for empirical data or non-botanical data. Use quantile regression estimates as a starting place, or if you already have domain-specific knownledge about your covariates, use this knowledge to set your prior hyperparameters.
#'
#' Note that this function is currently only implemented when the onset and duration models have the same covariates.
#'
#' This procedure is automatically applied when the runStanPhenology function parameter partitionDataForPriors is set to TRUE.
#' 
#' @return A list with a two element vector (mean and SD) for the mean onset (onsetHyperAnchor), a data frame for the mean slope coefficient (first column) and SD (second column) for each of the onset coefficients (onsetHyperBeta), a two element vector (mean and SD) for the mean duration (durationHyperAnchor), a data frame for the mean slope coefficient (first column) and SD (second column) for each of the duration coefficients (durationHyperBeta), a two element vector (mean and SD) for the mean cessation (cessationHyperAnchor), and a two element vector (mean and SD) for the sigma phenology parameter (sigmaHyper).
#' @export
#' @importFrom quantreg rq
#' @importFrom stats as.formula coef quantile
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
	H_SD_Anchor_meanOnset = 14*scale
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


#' Quantile estimates for hyperparameters of intercept-only (no covariates) models
#' 
#' This function uses quantiles to estimate the mean onset, mean cessation, and mean duration (mean onset - mean duration). The quantile approach for phenology was first implemented by Park et al. (2024) https://doi.org/10.1111/ecog.06961
#'
#' Estimates of the mean onset and mean duration values are prone to systematic errors due to fixing aribrary quantile levels. Park et al. (2024) set quantiles to 10% (onset) and 90% (cessation).
#'
#' @details 
#' This approach to define prior hyperparameters depends on arbitrary quantile cutoffs and, although it may work well in practice for simulated botanical collection data, it is not currently theoretically justified, nor is it known how well it works for empirical data or non-botanical data. Use quantile estimates as a starting place, or if you already have domain-specific knownledge, use this knowledge to set your prior hyperparameters.
#'
#' The returned data structures can be used directly to provide the prior hyperparameters to the runStanPhenology function.
#' 
#' This procedure is automatically applied when the runStanPhenology function parameter partitionDataForPriors is set to TRUE.
#' 
#' @param responseDataForPrior A vector of response data (e.g., day of year of collection values)
#' @param lowerQuantile The quantile used to estimate the onset model
#' @param upperQuantile The quantile used to estimate the cessation model
#' @param scale Set the scaling of hardcoded values for the standard deviation of the mean duration (two weeks), the standard deviation of the mean onset (two weeks), and the sigma mean (one week) and sigma standard deviation (half a week) values. For example, setting the scale to 0.5 would make the standard deviation of the prior distribution for mean onset equal to one week rather than two (default: 1)
#'
#' @return A vector with six hyperparameter values in the following order: mean of mean onset, standard deviation (SD) of mean onset, mean of mean duration, SD of mean duration, mean of sigma, SD of sigma.
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

	return(c(mean_mu_O = H_M_MO, SD_mu_O = H_SD_MO, mean_mu_D=H_M_MD, SD_mu_D=H_SD_MD, mean_sigma=H_M_S, SD_sigma=H_SD_S))
}


#' Partition response and covariate data
#'
#' Randomly partitions response and covariate data into two non-intersecting sets whose sizes are determined by the input proportion. This may be useful if you would like to use a subset of your data for one inference task (e.g., hyperparameter estimates) and the other subset for a separate inference task (e.g., Bayesian inference).
#'
#' @details
#' 
#' This procedure is automatically applied when the runStanPhenology function parameter partitionDataForPriors is set to TRUE at the default proportion of 0.3.
#' @param responseData A vector of response data (e.g., day of year of collection values)
#' @param onsetCovariateData  A data frame with the onset covariate data. This can be obtained from the preparePhenologyData function.
#' @param durationCovariateData  A data frame with the duration covariate data. This can be obtained from the preparePhenologyData function.
#' @param prop The proportion of the data to be used for the first set of the partition. (default: 0.3)
#'
#' @return A list of six elements. The first three represent the first subset for response data (responseDataForInference), onset covariate data (onsetCovariateDataForInference) and duration covariate data (durationCovariateDataForInference), and the last three correspondingly represent the second subset.
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


#' Partition data for models with no covariates (intercept-only)
#' Randomly partitions response data into two non-intersecting sets whose sizes are determined by the input proportion. This may be useful if you would like to use a subset of your data for one inference task (e.g., hyperparameter estimates) and the other subset for a separate inference task (e.g., Bayesian inference).
#'
#' @param responseData A vector of response data (e.g., day of year of collection values)
#' @param prop The proportion of the data to be used for the first set of the partition. (default: 0.3)
#'
#' @return A list with two vectors corresponding to the two sets of the partition.
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


#' Approximate the kth largest value in a normally-distributed finite population
#'
#' @description Uses an asymptotic approximation of order statistics for individuals in a finite population (Elfving, 1947; Cody, 1993) when the underlying population is normally distributed. For "exact", but slower, calculations, use the E.Ok1 function for first onset and the E.CkN function for last cessation. 
#'
#' Input values can be vectors.
#' @param N The population size
#' @param mu The mean of the population
#' @param sigma The standard deviation of the population
#' @param k The rank order of the individual. k=1 gives the first, and k=N gives the last.
#'
#' @return A vector of estimates of the expected value of the input kth ranked individual in the population.
#' @export
#' @importFrom stats qnorm
#'
#' @examples
#' \donttest{
#' #Define population sizes
#' N = c(1000,10000,100000)
#' #Define the mean values of the populations
#' mu = rep(100, length(N))
#' #Define the standard deviation of the populations
#' sigma = rep(7, length(N))
#' #Set the rank order
#' k = rep(1, length(N))
#' #Calculate the expected value of the individuals at the kth rank
#' eV = E.Kth.approx(N=N, mu=mu, sigma=sigma, k=k)
#' }
E.Kth.approx = function(N, mu, sigma, k) {
	if(any(k<0 | k>N | N<0)) {
		stop("k must be between 1 and N inclusive, and N must be a positive integer.")
	}
	approx = mu + sigma * qnorm((k - pi/8)/(N - pi/4 + 1))
	return(approx)
}


#Obtains the unscaled hyperparameters from a file associated with a provided file name
#' @importFrom utils read.table
getHyperparameters = function(hyperparameterFile) {
	hyperparameters = read.table(hyperparameterFile, header=TRUE, sep='\t')
	return(hyperparameters)
}

#Obtains min max scaled covariate data from a file associated with a provided file name
#' @importFrom utils read.table
getCovariates = function(covariatesFile) {
	covariates = read.table(covariatesFile, header=TRUE, sep='\t')
	return(processCovariates(covariates))
}

#' @importFrom stats sd
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

#' @importFrom stats approxfun
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

#' @importFrom utils read.table
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
#' @importFrom stats dbeta
dbeta_safe = function(x, shape1, shape2) {
	ifelse(x <= 0 | x >= 1, 0, dbeta(x, shape1, shape2))
}

#' @importFrom stats pbeta
pbeta_safe = function(q, shape1, shape2) {
	ifelse(q <= 0, 0,
		   ifelse(q >= 1, 1, pbeta(q, shape1, shape2)))
}


#' Summary diagnostics of a Stan run
#' 
#' @description 
#' Provides a data frame with diagnostic information for a Stan run.
#'
#' @param stanRunResult The output from the runStanPhenology function
#' @param NTotal The sample size (i.e., number of specimens used in analysis)
#' @param taxonName The taxon name
#' @param diagnostics A vector of named diagnostic quantities (default: c("num_divergent", "num_max_treedepth", "ebfmi"))
#'
#' @return A list with the taxon name, sample size total, sample size used by Stan, and a data frame with the Stan diagnostics summary table
#' @export
#' @importFrom dplyr %>%
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

#' Fit a standard linear model to the observed collection times, T
#'
#' Fits a standard linear model to the observed collection times based on the covariates used to model the onset and duration. 
#' 
#' Data objects for the input to this function can be obtained from the preparePhenologyData function and are provided with the output of the Stan run. 
#'
#' This function is run automatically when the standardLinearModel parameter of the summarizePhenologyResults function is set to TRUE, with corresponding results provided in the output data frame columns 'standard.linear.model'  and 'in.95CI'.
#' @param responseData A vector of the response data.
#' @param onsetCovariateData A data frame with each named column providing a covariate for the onset model. 
#' @param durationCovariateData A data frame with each named column providing a covariate for the duration model.
#'
#' @return The fit linear model object that is the output of the R lm function.
#' @export
runStandardLinearModel = function(responseData, onsetCovariateData, durationCovariateData) {
        # Combine response and predictors into one data frame
        merged = merge_df_by_column(onsetCovariateData, durationCovariateData)
        df <- data.frame(response = responseData, merged)

        # Dynamically construct formula
        formula <- as.formula(paste("response ~", paste(colnames(merged), collapse = " + ")))

        # Fit the model
        fit <- lm(formula, data = df)

        print(summary(fit))
        return(fit)
}


#' Summary of Stan parameter estimates and associated diagnostics for a phenology model with covariates
#' 
#' Provides a summary as a data frame of the results of a Stan analysis of phenology. The model is assumed to include covariates with corresponding slope coefficients, mean onset, mean duration values. 
#'
#' @param stanRunResult The output from runStanPhenology 
#' @param taxonName The name of the taxon being analyzed
#' @param measures The metrics to calculate for each parameter. (default: c("mean", "median", "sd", "mad"))
#' @param quantiles The quantiles to calculate for each metric. (default: c(2.5, 97.5))
#' @param convergence Stan sample diagnostic values to calculate for each metric. (default: c("rhat","ess_bulk", "ess_tail"))
#' @param standardLinearModel Optional Boolean indicating whether to fit a standard linear regression model through the observed collection times. Only used when the onset and duration covariates are the same. (default: FALSE)
#' @param processExtremes Boolean indicating whether to process phenological extremes (default: TRUE)
#' @param N The population size used when processing extremes. Not used when processExtremes = FALSE. (default: 500)
#'
#' @return A data frame with the the following columns: taxon name, variable name, the type of phenological parameter being summarized (e.g., O onset), the name of the covariate, the marginal posterior probability that the slope coefficient is negative, the estimate for the observed collection time standard linear model, whether the estimate based on the standard linear model is in the 95%CI of the Bayesian estimate, metrics, diagnostics.
#' @export
#' @importFrom posterior as_draws_df quantile2
#' @importFrom dplyr %>%
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
#' data  =  preparePhenologyData(dataFile=file, responseVariableName="DOY", onsetCovariateNames=vars
#'                              , durationCovariateNames=vars, taxonName="Sanguinaria_canadensis"
#'                              , removeOutliers=TRUE)
#' ##run the Stan sampler
#' stanResult  =  runStanPhenology(type="full", responseData = data$responseData
#'                                 , onsetCovariateData = data$onsetCovariateData
#'                                 , durationCovariateData = data$durationCovariateData
#'                                 , partitionDataForPriors = TRUE)
#' ##summarize the Stan run
#' stanSummary  =  summarizePhenologyResults(stanRunResult = stanResult
#'                                           , taxonName = "Sanguinaria_canadensis"
#'                                           ,standardLinearModel = TRUE)
#' stanSummary
#' }
summarizePhenologyResults = function(stanRunResult, taxonName, measures=c("mean", "median", "sd", "mad"), quantiles=c(2.5, 97.5), convergence=c("rhat","ess_bulk", "ess_tail"), standardLinearModel=FALSE, processExtremes=TRUE, N=500) {
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

	if(standardLinearModel) {
		standardLinearModel = runStandardLinearModel(stanRunResult$responseData, stanRunResult$onsetCovariateData, stanRunResult$durationCovariateData) 
	}

	
	# if(class(standardLinearModel)=="lm") {
	# DANIEL: Changing to inherits because it is safer:
	if( inherits(x = standardLinearModel, what = "lm") ) {
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
		#if(class(standardLinearModel)=="lm") {
		# DANIEL: Changing to inherits because it is safer:
		if( inherits(x = standardLinearModel, what = "lm") ) {
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
	# if(class(standardLinearModel)=="lm") {
	# DANIEL: Changing to inherits because it is safer:
	if( inherits(x = standardLinearModel, what = "lm") ) {
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

