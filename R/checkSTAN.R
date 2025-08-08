## Functions that will run on package loading:
## Make some functions to check if Stan variables are configured:

## Ask "q" and get yes or no answer.
getAnswer <- function(q = "Do you want to continue? (yes/no): "){
  response <- readline(prompt = q)
  response <- match.arg(arg = response, choices = c("yes", "YES", "no", "NO"))
  if(response %in% c("yes", "YES")) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Check STAN and fix dependencies
#'
#' @description Function will conduct checks to verify if the system has the required dependencies to install and run STAN together with "phenoCollectR". The function will attempt to correct the dependencies and install STAN if necessary. The function is interactive and will prompt the user by printing information to the console.
#'
#' @return Return TRUE or FALSE (using the "invisible()" function) when the system is ready (or not) to do analyses using STAN.
#' 
#' @importFrom cmdstanr check_cmdstan_toolchain cmdstan_version install_cmdstan
#' @importFrom utils install.packages
#' @export
makeSTANchecks <- function(){
  ## Check if it is Windows
  if( isTRUE(.Platform$OS.type == "windows") ){
    message("Running phenoCollectR in Windows.")
    fix_tools <- TRUE
  } else{
    message("Running phenoCollectR in Linux or Mac.")
    fix_tools <- FALSE
  }
  message("Checking dependencies...")
  if( !requireNamespace("cmdstanr", quietly = TRUE) ){
    message("The cmdstanr package is not installed. Should try to install it?")
    if( getAnswer() ){
      ## Go ahead
      install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
      ## Check if successful:
      if( !requireNamespace("cmdstanr", quietly = TRUE) ){
        message("Failed to install cmdstanr package.")
        message("Please follow installation instructions from https://mc-stan.org/cmdstanr/")
        return( invisible(FALSE) )
      } else{
        message("Dependency package installed successfully!")
      }
    } else{
      message("Please install the cmdstanr package before running analyses.")
      message("Follow instructions on https://mc-stan.org/cmdstanr/")
      return( invisible(FALSE) )
    }
  }
  cmdstanr::check_cmdstan_toolchain(fix = fix_tools)
  if( is.null( cmdstanr::cmdstan_version(error_on_NA = FALSE) ) ){
    message("CmdStan program not found. Should try to install it?")
    if( getAnswer() ){
      ## Go ahead and install STAN
      cmdstanr::install_cmdstan()
    } else{
      message("Please install CmdStan before running analyses.")
      return( invisible(FALSE) )
    }
  }
  message("All dependencies found. phenoCollectR is ready!")
  return( invisible(TRUE) )
}

## This is a TRUE or FALSE return function that silently checks for STAN dependencies.
#' @importFrom cmdstanr check_cmdstan_toolchain cmdstan_version install_cmdstan
#' @noRd
makeSTANpassChecks <- function(){
  ## has_* is FALSE if missing
  has_cmdstanr <- requireNamespace("cmdstanr", quietly = TRUE)
  has_toolchain <- cmdstanr::check_cmdstan_toolchain(fix = FALSE, quiet = TRUE)
  has_cmdstan <- !is.null( cmdstanr::cmdstan_version(error_on_NA = FALSE) )
  passed <- all(has_cmdstanr, has_toolchain, has_cmdstan)
  ## TRUE if all checks passed.
  return( passed )
}
