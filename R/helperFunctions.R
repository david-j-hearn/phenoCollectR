
rdirichlet <- function(n, alpha) {
	K <- length(alpha)
	G <- matrix(rgamma(n * K, shape = alpha, rate = 1), nrow = n, ncol = K, byrow = TRUE)
	G / rowSums(G)
}

softplus <- function(x) ifelse(x > 30, x, log1p(exp(x)))

true_marginal_line <- function(alpha, beta, mu, j,
			       Sigma = NULL,
			       R = NULL,
			       sd_x = NULL) {

	C <- length(beta)
	stopifnot(length(mu) == C)
	stopifnot(j >= 1 && j <= C)

	# --- Build Sigma if needed ---
	if (is.null(Sigma)) {

		# Must have R and sd_x
		if (is.null(R) || is.null(sd_x)) {
			stop("Provide either Sigma, OR (R and sd_x).")
		}

		stopifnot(all(dim(R) == c(C, C)))
		stopifnot(length(sd_x) == C)

		# Convert correlation + SDs to covariance
		Sigma <- diag(sd_x) %*% R %*% diag(sd_x)

	} else {
		stopifnot(all(dim(Sigma) == c(C, C)))
	}

	# --- Partition indices ---
	idx_other <- setdiff(seq_len(C), j)

	beta_j    <- beta[j]
	beta_rest <- beta[idx_other]

	mu_j      <- mu[j]
	mu_rest   <- mu[idx_other]

	# --- Covariance blocks ---
	Sigma_jj      <- Sigma[j, j]             # scalar
	Sigma_rest_j  <- Sigma[idx_other, j]     # (C-1) x 1 vector

	# --- Marginal line parameters ---
	slope <- beta_j + as.numeric(t(beta_rest) %*% Sigma_rest_j / Sigma_jj)

	intercept <- alpha + as.numeric(
					t(beta_rest) %*% (mu_rest - Sigma_rest_j * mu_j / Sigma_jj)
	)

	list(intercept = intercept, slope = slope, Sigma = Sigma)
}


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
		if(is.data.frame(onsetCovariateData)) { print("started with a data frame in onsetCovariateData. This should not be set if using automatic priors.") }
		stop("Partitioning didn't create a data frame.")
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

	#if(any(maxs-mins==0)) {
	#stop("Constant covariate data is not allowed.")
	#}

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
#' stanSummary  =  summarizePhenologyResults(stanRunResult = stanResult,
#'                                           taxonName = "Sanguinaria_canadensis"
#'                                           )
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
		summary[summary$variable==varT, "type"] = "observed times (T; presence only)"
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
	summary[summary$variable=="anchor_T", "type"] = "observed (T; presence only)"

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
	summary[summary$variable=="alpha_T", "type"] = "observed (T; presence only)"
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

checkInput = function(type=c("intercept-only","full","multistage-full"), responseData=NULL, onsetCovariateData=NULL, durationCovariateData=NULL, stage=NULL, nStages=NULL, nOnsetCovariates=NULL, nDurationCovariates=NULL, minResponse, maxResponse, maxDiv, N, processExtremes) {

	type = match.arg(type)
	N = length(responseData)
	N_C = nrow(onsetCovariateData)
	K = ncol(onsetCovariateData)
	N_O = nrow(onsetCovariateData)
	N_D = nrow(durationCovariateData)

	if(!(type=="intercept-only" || type=="full" || type=="multistage-full" )) {
		cat(paste("Unsupported type: ", type, "\nType should be 'intercept-only' or 'full' or 'multistage-full'.\n"))
		stop("Unsupported type error.")
	}

	if(processExtremes && N<0) {
		stop("A positive population size, N, must be supplied when processing extremes.")
	}

	if(maxDiv > 0) {
		warning("The maximum number of divergences tolerated is greater than 0. This can result in biased estimates.")
	}

	if(minResponse != 0){
		stop("Under current implementations, the minimum response time (minResponse) must be set to 0, which is the default.")
	}

	if(maxResponse <= minResponse) {
		stop("The maximum possible response time must be larger than the minimum response time of 0.")
	}

	if(!is.vector(responseData)) {
		stop("Expecting a vector of real numeric values of the collection times (e.g., day of year (DOY) of when specimens were collected).")
	}

	if(length(responseData)<10) {
		warning("Ten or fewer data items is a very small sample size and will likely result in inaccurate inferences and a high divergence rate during Bayesian inference. The sample size should be at least 60, especially when covariates are used.")
	}

	if(type=="multistage-full") {
		#check that stage values match
		N_S = length(stage)
		S = nStages
		if(N != N_S) {
			stop("Make sure that there is a stage for each observation, and vice versa.")

			#return(list(error=TRUE, error_m="When running the multistage model, please include a vector of stages of the individuals that were observed. The vector of stages should be the same length as the vector of response data. For each individual sampled, there should be a response time and there should be a stage associated with the time. If the time is before the fist phenophase of the time period, the stage is the last stage, if during the first phenophase, the stage is 1, and so forth."))
		}
		if(nStages < nlevels(factor(stage))-1) {
			stop("The number of input stages (nStages) is less than the number of different stages in 'stage'. Data can be missing for a stage, but there cannot be more actual stages than reported stages.")
		}
		if(max(stage)-1 > nStages) {
			stop("The input maximum stage is more than the reported number of stages (nStages).")
		}
		if(length(responseData) != length(stage)) {
			stop("The stage data vector must be the same length as the response data vector so that each time in the response data corresponds to a stage. Stages should be coded with integer values starting with 1 and increasing consecutively.")
		}

		if(!identical(onsetCovariateData, durationCovariateData)) {
			stop("The covariate data for duration models must be the same as the covariate data for the onset model when applying type 'multistage-full'.")
		}

		if(nOnsetCovariates != ncol(onsetCovariateData)) {
			stop("Input number of onset covariates does not match.")
		}

		if(nDurationCovariates != ncol(durationCovariateData)) {
			stop("Input number of duration covariates does not match.")
		}
	}

	if(type=="full" || type=="multistage-full") {
		if(sum(is.na(onsetCovariateData)) || !is.data.frame(onsetCovariateData) ) {
			cat("Please remove all NA values, and provide a data frame of the onset model covariate (predictor variable) data. \nThese might be temperature or precipitation data at the specimen collection sites, for example.\nEach column of the data frame should be named with the predictor variable name (e.g. 'meanAnnualTemperature').\nThe order of the rows should correspond to the order of the elements in the response variable data vector (i.e., the first element in the response vector corresponds with the first row in the onset covariate data frame, the second element with the second row, and so forth).\nPlease see documentation for additional information and examples.")
			stop("Please provide appropriate inputs")
		}
		if(sum(is.na(durationCovariateData)) || !is.data.frame(durationCovariateData) ) {
			cat("Please remove all NA values, and provide a data frame of the duration model covariate (predictor variable) data. \nThese might be temperature or precipitation data at the specimen collection sites, for example.\nEach column of the data frame should be named with the predictor variable name (e.g. 'meanAnnualTemperature').\nThe order of the rows should correspond to the order of the elements in the response variable data vector (i.e., the first element in the response vector corresponds with the first row in the duration covariate data frame, the second element with the second row, and so forth).\nPlease see documentation for additional information and examples.")
			stop("Please provide appropriate inputs")
		}
		cat("Checking data compatibility.\n")
		#check data dimensions
		if(N != N_O || N != N_D) {
			return(list(error=TRUE, error_m="Be sure the sample size is the the same for the response (observed) values, for the covariate values for the onset, and for the covariate values for the duration. Rows should be parallel. Row 1 in the response data vector should correspond to the same individual in row 1 of each covariate data frame. If files names are input, corresponding files should be text formatted with tab-separated columns which should have headers. If data frames are input, these should have column names."))
		}

		#check data dimensions
		if(N != N_C) {
			stop("Make sure that each covariate value is present for each observation, and vice versa.")
			#return(list(error=TRUE, error_m="Be sure the sample size is the the same for the response (observed) values and for the covariate values. Rows should be parallel. Row 1 in the response data vector should correspond to the same individual in row 1 of the covariate data frame. The input data frames should have column names."))
		}
	}
	return(TRUE)
}

checkPriors = function(type=c("intercept-only","full","multistage-full"), nStages, responseData=NULL, hyperparams_noCovariates=NULL, onsetCovariateData=NULL, durationCovariateData=NULL, onsetHyperBeta=NULL, onsetHyperBetaMean=NULL, onsetHyperBetaSD=NULL, onsetHyperAnchor=NULL, durationHyperBeta=NULL, durationHyperBetaMean=NULL, durationHyperBetaSD=NULL, durationHyperAnchor=NULL, sigmaHyper=NULL, partitionDataForPriors=FALSE, minResponse, maxResponse) {

	if(partitionDataForPriors) {
		cat("The data will be partitioned into two sets with 30% and 70% of the data. \n\n30% will be used to carry out a preliminary analysis using quantiles to estimate the prior distribution hyperparameter values. \n\n70% of the data will be used to carry out a Stan Bayesian analysis to obtain the posterior distributions of parameters.\n\nIf other hyperparameter information was provided as input, it will be ignored. \n\nThe calculated values based on quantiles are approximate; you may need to use other sources of data to get better estimates of prior hyperparameter values, especially if the Stan run results in divergences or other poor diagnostics.\n\n")
		prop = 0.3
		scale = (maxResponse-minResponse) / (365 - 0) 	#setting scale relative to parameters for annual variation


		if(type=="multistage-full") {
			stop("Automatic partitioning of data for prior hyperparameter estimation is not available for multistage analysis.")
		}

		if(type=="intercept-only") {
			warning("Automated hyperparameters are not recommended for data without covariates. Estimates of duration are likely to be inaccurate.")
			partition = partitionResponseData(responseData = responseData, prop = prop)
			responseData = partition$dataForInference
			hyperparams_noCovariates = getHyperparametersViaQuantiles(responseDataForPrior = partition$dataForPrior, scale = scale)
			return( list(
				     responseData=responseData,
				     hyperparams_noCovariates=hyperparams_noCovariates, 
				     onsetCovariateData=onsetCovariateData,
				     durationCovariateData=durationCovariateData,
				     onsetHyperBeta=onsetHyperBeta,
				     onsetHyperBetaMean=onsetHyperBetaMean,
				     onsetHyperBetaSD=onsetHyperBetaSD,
				     onsetHyperAnchor=onsetHyperAnchor,
				     durationHyperBeta=durationHyperBeta,
				     durationHyperBetaMean=durationHyperBetaMean,
				     durationHyperBetaSD=durationHyperBetaSD,
				     durationHyperAnchor=durationHyperAnchor,
				     sigmaHyper=sigmaHyper
	)
			)
		}
		else if(type=="full") {
			#partition data
			partition = partitionResponseCovariateData(responseData=responseData, onsetCovariateData=onsetCovariateData, durationCovariateData=durationCovariateData, prop=prop)

			responseData = partition$responseDataForInference
			onsetCovariateData = partition$onsetCovariateDataForInference
			durationCovariateData = partition$durationCovariateDataForInference

			#get the data for prior and set the prior hyperparameters
			prior = getHyperparametersViaQuantileRegression(responseDataForPrior=partition$responseDataForPrior, onsetCovariateDataForPrior=partition$onsetCovariateDataForPrior, durationCovariateDataForPrior=partition$durationCovariateDataForPrior, lowerQuantile=0.1, upperQuantile=0.9)

			#set the prior hyperparameters
			onsetHyperBeta = prior$onsetHyperBeta
			onsetHyperAnchor = prior$onsetHyperAnchor
			durationHyperBeta = prior$durationHyperBeta
			durationHyperAnchor = prior$durationHyperAnchor
			sigmaHyper = prior$sigmaHyper

			return( list(
				     responseData=responseData,
				     hyperparams_noCovariates=hyperparams_noCovariates, 
				     onsetCovariateData=onsetCovariateData,
				     durationCovariateData=durationCovariateData,
				     onsetHyperAnchor=onsetHyperAnchor,
				     onsetHyperBeta=onsetHyperBeta,
				     onsetHyperBetaMean=onsetHyperBetaMean,
				     onsetHyperBetaSD=onsetHyperBetaSD,
				     durationHyperAnchor=durationHyperAnchor,
				     durationHyperBeta=durationHyperBeta,
				     durationHyperBetaMean=durationHyperBetaMean,
				     durationHyperBetaSD=durationHyperBetaSD,
				     sigmaHyper=sigmaHyper
			)
			)
		}
	}

	if(type=="multistage-full") {
		nCovariates = ncol(onsetCovariateData)
		range = (maxResponse-minResponse)
		mr = mean(responseData)
		sdr = sd(responseData)
		oa = range/(nStages+1) 
		da = range/(nStages+1)
		sdx = rep(0,nCovariates)
		for(i in 1:nCovariates) {
			sdx[i] = sd(onsetCovariateData[,i])
		}
		#fill in default values if not set
		if(is.null(onsetHyperBeta)) { 
			cat("Automatically setting prior for onset model slope and SD.\n")
			onsetHyperBeta = data.frame(mean=rep(0,nCovariates), sd=1*sdr/sdx)
		}
		if(is.null(onsetHyperAnchor)) {
			cat("Automatically setting prior for onset model anchor.\n")
			onsetHyperAnchor=c(oa,1*sdr)	#sd arbitrarily set to 1/10 of expected duration
		}
		if(is.null(durationHyperBetaMean)) {
			cat("Automatically setting prior for duration model slope means.\n")
			durationHyperBetaMean=matrix(rep(0,(nStages-1)*nCovariates), nrow=nStages-1)
		}
		if(is.null(durationHyperBetaSD)) {
			cat("Automatically setting prior for duration model slope SDs.\n")
			durationHyperBetaSD=t(matrix(rep(1*sdr/sdx,(nStages-1)), nrow=nStages-1))
		}
		if(is.null(durationHyperAnchor)) {
			cat("Automatically setting prior for duration model anchors.\n")
			durationHyperAnchor = data.frame(mean = rep(da,nStages-1), sd=rep(1*sdr,nStages-1))	#sd arbitrarily set to 1/10 of expected duration
		}
		if(is.null(sigmaHyper)) {
			cat("Automatically setting prior for sigma.\n")
			sigmaHyper=c(0,0.1*sdr)
			#sigmaHyper=data.frame(mean=rep(0,nStages),sd=rep(0.1 * mr/sdr,nStages)) #half normal always positive (except for measure 0)
		}

		print("Hyperparameters: sd response data, oa, ob, da, db_m, db_sd, s")
		print(sdr)
		print(onsetHyperAnchor)
		print(onsetHyperBeta)
		print(durationHyperAnchor)
		print(durationHyperBetaMean)
		print(durationHyperBetaSD)
		print(sigmaHyper)

		if(nrow(durationHyperBetaSD) != nrow(durationHyperBetaMean)) {
			stop("The number of stages (rows) in the durationHyperBetaSD and durationHyperBetaMean matrices need to be the same.")
		}

		K = nCovariates
		S = nStages

		#check hyperparameter specifications
		if(nrow(onsetHyperBeta) != K || ncol(onsetHyperBeta) != 2 || nrow(durationHyperBetaMean) != (S-1) || ncol(durationHyperBetaMean) != K || nrow(durationHyperBetaSD) != (S-1) || ncol(durationHyperBetaSD) != K || length(onsetHyperAnchor) != 2 || nrow(durationHyperAnchor) != S-1 || ncol(durationHyperAnchor) != 2 || length(sigmaHyper) != 2) {
			cat("Hyperparameters should be as follows:\n onsetHyperAnchor: a 2-element vector with mean and sd of the onset time;\n onsetHyperBeta: a data frame with two columns (mean and sd of each covariate slope of the onset model), one covariate per row;\n durationHyperAnchor: a data frame with two columns (mean and sd of the duration), one pair for each stage(rows), excluding the last stage (S-1 rows);\n durationHyperBetaMean: a matrix with S-1 rows and K columns with each stage X covariate mean slope in cells;\n durationHyperBetaSD: a matrix with S-1 rows and K columns with each stage X covariate slope sd in cells;\n sigmaHyper: a data frame with two colums (mean and sd of each model sigma) and S rows, one for each stage.")
			stop("Please adjust hyperparameter inputs.")
		}

		#NULL's should all be set now. Checking datatypes
		if(!is.data.frame(onsetHyperBeta) || !is.matrix(durationHyperBetaMean) || !is.matrix(durationHyperBetaSD) || !is.vector(onsetHyperAnchor) || !is.data.frame(durationHyperAnchor) || !is.vector(sigmaHyper)) {
			cat("Expecting the following:\n\tonsetHyperBeta:\n\t\tA data frame with two columns. The first column is the mean hyperparameter for the onset slope coefficient for each covariate of the first stage. The second column is the standard deviation of the onset slope coefficient for each covariate of the first stage. The first row in the hyperparameters data frame is the first covariate corresponding to the first column in the covariate file, the second row with the second covariate, and so forth. \n\tdurationHyperBetaMean:\n\t\tA matrix with dimensions (number of stages -1) X (number of covariates). The expected slope value for each stage X covariate combination goes in the matrix cells. \n\tdurationHyperBetaSD:\n\t\tA matrix with dimensions (number of stages -1) X (number of covariates). The SD of each slope value for each stage X covariate combination goes in the matrix cells. \n\tonsetHyperAnchor:\n\t\tA two-element vector with the mean of the prior and the standard deviation of the prior for the onset anchor (the mean onset value when no covariate data are included) of the model for the first stage.\n\tdurationHyperAnchor:\n\t\tA data frame with two columns, with the mean of the prior and the standard deviation of the prior for the duration anchor (the mean duration value when no covariate data are included) in respective columns. Each row is for the duration model for each stage except for the last stage.\n\tsigmaHyper:\n\t\tA vector with two items, theman and the standard deviation of the onset sigma. \nSee documentation for additional information and examples")
			stop("Please provide appropriate inputs")
		}
		#all is good - return the information
		return( list(
			     responseData=responseData,
			     hyperparams_noCovariates=hyperparams_noCovariates, 
			     onsetCovariateData=onsetCovariateData,
			     durationCovariateData=durationCovariateData,
			     onsetHyperBeta=onsetHyperBeta,
			     onsetHyperBetaMean=onsetHyperBetaMean,
			     onsetHyperBetaSD=onsetHyperBetaSD,
			     onsetHyperAnchor=onsetHyperAnchor,
			     durationHyperBeta=durationHyperBeta,
			     durationHyperBetaMean=durationHyperBetaMean,
			     durationHyperBetaSD=durationHyperBetaSD,
			     durationHyperAnchor=durationHyperAnchor,
			     sigmaHyper=sigmaHyper
			)
		)
	}

	if(type=="full") {
		if(!is.data.frame(onsetHyperBeta) || !is.data.frame(durationHyperBeta) || !is.vector(onsetHyperAnchor) || !is.vector(durationHyperAnchor) || !is.vector(sigmaHyper)) {
			cat("Expecting the following:\n\tonsetHyperBeta:\n\t\tA data frame with two columns. The first column is the mean hyperparameter for the onset slope coefficient for each covariate. The second column is the standard deviation of the onset slope coefficient for each covariate. The first row in the hyperparameters data frame is the first covariate corresponding to the first column in the covariate file, the second row with the second covariate, and so forth. \n\tdurationHyperBeta:\n\t\tA data frame with two columns. The first column is the mean hyperparameter for the duration slope coefficient for each covariate. The second column is the standard deviation of the duration slope coefficient for each covariate. The first row in the hyperparameters data frame is the first covariate corresponding to the first column in the covariate file, the second row with the second covariate, and so forth. \nonsetHyperAnchor:\n\t\tA two-element vector with the mean of the prior and the standard deviation of the prior for the onset anchor (the mean onset value when no covariate data are included).\n\tdurationHyperAnchor:\n\t\tA two-element vector with the mean of the prior and the standard deviation of the prior for the duration anchor (the mean duration value when no covariate data are included).\n\tsigmaHyper:\n\t\tA two-element vector with the mean of the prior and the standard deviation of the prior for the sigma model parameter (variation in onset times and variation in cessation times).\nSee documentation for additional information and examples")
			stop("Please provide appropriate inputs")
		}
		K_O = ncol(onsetCovariateData)
		K_D = ncol(durationCovariateData)
		#check hyperparameter specifications
		if(nrow(onsetHyperBeta) != K_O || ncol(onsetHyperBeta) != 2 || nrow(durationHyperBeta) != K_D || ncol(durationHyperBeta) != 2 || length(onsetHyperAnchor) != 2 || length(durationHyperAnchor) != 2 || length(sigmaHyper) != 2) {
			return(list(error=TRUE, error_m="The input hyperparameter information should be a data frame with one row for each covariate, or provide a file name with corresponding file having the data frame as a text, tab-separated table. Each row should have two columns. The first column provides the mean value of the parameter's prior distribution, and the second column provides the SD of the parameter's prior distribution. The covariates in rows, top to bottom, should match the covariates in columns, left to right, in the input covariate data frames. Values should be in the original scale (e.g., units of days, or for slopes, total days changed over range of the covariate) of the observations. If files names are input, corresponding files should have headers. If data frames are input, these should have column labels of 'mean
				    and 'sd'."))

		}

		#all is good - return the information
		return( list(
			     responseData=responseData,
			     hyperparams_noCovariates=hyperparams_noCovariates, 
			     onsetCovariateData=onsetCovariateData,
			     durationCovariateData=durationCovariateData,
			     onsetHyperBeta=onsetHyperBeta,
			     onsetHyperBetaMean=onsetHyperBetaMean,
			     onsetHyperBetaSD=onsetHyperBetaSD,
			     onsetHyperAnchor=onsetHyperAnchor,
			     durationHyperBeta=durationHyperBeta,
			     durationHyperBetaMean=durationHyperBetaMean,
			     durationHyperBetaSD=durationHyperBetaSD,
			     durationHyperAnchor=durationHyperAnchor,
			     sigmaHyper=sigmaHyper
		    )
		)
	}

	if(type=="intercept-only") {
		if(sum(is.na(hyperparams_noCovariates)) || length(hyperparams_noCovariates) != 6) {
			stop("Expecting six hyperparameter values (mean and sd for mean onset, mean and sd for mean duration, mean and sd for sigma. Or, if you want hyperparameter values to be estimated for you, set 'partitionDataForPriors' to TRUE. Automated hyperparameter estimation is not recommended for these model types.")
		}
		if(hyperparams_noCovariates[2] < 0 || hyperparams_noCovariates[4] < 0 || hyperparams_noCovariates[6] < 0) {
			stop("SD hyperparameter values should be positive. Non-positive values detected in hyperparams_noCovariates.")
		}
		#all is good - return the information
		return( list(
			     responseData=responseData,
			     hyperparams_noCovariates=hyperparams_noCovariates, 
			     onsetCovariateData=onsetCovariateData,
			     durationCovariateData=durationCovariateData,
			     onsetHyperBeta=onsetHyperBeta,
			     onsetHyperBetaMean=onsetHyperBetaMean,
			     onsetHyperBetaSD=onsetHyperBetaSD,
			     onsetHyperAnchor=onsetHyperAnchor,
			     durationHyperBeta=durationHyperBeta,
			     durationHyperBetaMean=durationHyperBetaMean,
			     durationHyperBetaSD=durationHyperBetaSD,
			     durationHyperAnchor=durationHyperAnchor,
			     sigmaHyper=sigmaHyper
		)
		)
	}
	stop("Checking hyperparameters failed for unknown reason.")
}
