# phenoCollectR
R package for Bayesian analysis of phenological timing using biocollection data

# Installation 
```{r,R.options=list(max.print=20)}
install.packages("devtools")
devtools::install_github("david-j-hearn/phenoCollectR", dependencies = TRUE)
```

## Installation notes
`phenoCollectR` depends on STAN. The package will check if all dependencies to run STAN are available and will ask the user if STAN should be installed.

STAN installation (and checks) are performed using the utilities of the `cmdstanr` package. `phenoCollectR` will also install `cmdstanr` if necessary.

## Quick example:
```{r,R.options=list(max.print=20)}
## Get the file name with data for the blood root plant
file <- system.file("data", "Sanguinaria_canadensis.Full.txt", package = "phenoCollectR")
## Define the covariate names - remove up to all but 1
vars = c("Latitude", "Year", "Elevation", "AnnualMonthlyAverageTemp", "SpringMonthlyAverageTemp", "FirstQuarterMonthlyAverageTemp")
## Get the phenology data
data <- preparePhenologyData(dataFile=file, responseVariableName="DOY", onsetCovariateNames=vars, durationCovariateNames=vars, taxonName="Sanguinaria_canadensis", removeOutliers=TRUE)
## Run the Stan sampler
stanResult <- runStanPhenology(type="full", responseData = data$responseData, onsetCovariateData = data$onsetCovariateData, durationCovariateData = data$durationCovariateData, partitionDataForPriors = TRUE)
## Summarize the Stan run
stanSummary <- summarizePhenologyResults(stanRunResult = stanResult, taxonName = "Sanguinaria_canadensis",standardLinearModel = TRUE)
## Make posterior predictive graph
pp <- makePosteriorPredictivePlot(stanResult = stanResult, responseData = data$responseData, targetCovariateName = "SpringMonthlyAverageTemp", onsetCovariateData = data$onsetCovariateData, durationCovariateData = data$durationCovariateData)
## Display the posterior predictive graph		
print(pp)
```
