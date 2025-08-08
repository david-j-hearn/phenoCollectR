#' @importFrom stats pnorm
obs_gp_nc <- function(mu_O, duration, sigma) {
  mu_C <- mu_O + duration
  mu_C_m1_over_sigma = (mu_C - 1) / sigma;
  mu_C_over_sigma = mu_C / sigma;
  mu_O_m1_over_sigma = (mu_O - 1) / sigma;
  mu_O_over_sigma = mu_O / sigma;
  
  exp_sum = -exp(-0.5 * mu_C_m1_over_sigma^2) + exp(-0.5 * mu_C_over_sigma^2) +
    exp(-0.5 * mu_O_m1_over_sigma^2) - exp(-0.5 * mu_O_over_sigma^2);
  scaled_exp_sum = sqrt(1.0/(2 * pi)) * sigma * exp_sum;
  mu_C_terms = -(mu_C - 1) * (pnorm(mu_C_m1_over_sigma) - 0.5) +
    mu_C * (pnorm(mu_C_over_sigma) - 0.5);
  mu_O_terms = -(mu_O - 1) * pnorm(-mu_O_m1_over_sigma) +
    mu_O * pnorm(-mu_O_over_sigma);
  
  -0.5 + scaled_exp_sum + mu_C_terms + mu_O_terms;    
}

robs_gp <- function(n, mu_O, duration, sigma, max_rejection = 10000) {
  res <- rep(NA_real_, n)
  mu_C <- mu_O + duration
  if(length(mu_O) == 1) {
    mu_O <- rep(mu_O, n)
  }
  if(length(mu_C) == 1) {
    mu_C <- rep(mu_C, n)
  }
  if(length(sigma) == 1) {
    sigma <- rep(sigma, n)
  }
  
  for(i in 1:max_rejection) {
    na_indices <- is.na(res)
    n_na <- sum(na_indices)
    if(n_na == 0) {
      break
    }
    t <- runif(n_na)
    gp_t_all <- rnorm(n_na, t, sigma[na_indices])
    seen <- mu_O[na_indices] < gp_t_all & gp_t_all < mu_C[na_indices]
    t[!seen] <- NA_real_
    res[na_indices] <- t
  }
  if(any(is.na(res))) {
    warning("Could not produce viable observations for some indices by rejection sampling.")
  }
  res
}

#' @importFrom brms get_dpar
#' @noRd
posterior_predict_obs_gp <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  duration <- brms::get_dpar(prep, "duration", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(robs_gp(prep$ndraws, mu, duration, sigma))
}

#' @importFrom brms get_dpar
#' @noRd
posterior_epred_obs_gp <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  duration <- brms::get_dpar(prep, "duration")
  return(mu + 0.5 * duration)
}

#' @importFrom brms get_dpar
#' @importFrom dplyr if_else
#' @importFrom stats pnorm
#' @noRd
log_lik_obs_gp <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  duration <- brms::get_dpar(prep, "duration", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  t <- prep$data$Y[i]
  
  mu_O <- mu
  mu_C <- mu + duration
  base_ll <- 
    dplyr::if_else(t > mu_O, 
      brms:::log_diff_exp(pnorm(mu_C, t, sigma, log.p = TRUE), pnorm(mu_O, t, sigma, log.p = TRUE)),
      brms:::log_diff_exp(pnorm(2 * t - mu_O, t, sigma, log.p = TRUE), pnorm(2 * t - mu_C, t, sigma, log.p = TRUE))
    )
  nc <- obs_gp_nc(mu_O, duration, sigma)                            
  ll <- base_ll - log(nc)

  ll
}

#' @importFrom brms custom_family stanvar
#' @importFrom dplyr if_else
#' @noRd
obs_gp_brmsfamily <- function(link = "identity", link_duration = "log", link_sigma = "log") {
  family <- brms::custom_family(
    "obs_gp",
    dpars = c("mu", "duration", "sigma"),
    links = c(link, link_duration, link_sigma),
    lb = c(NA, 0, 0),
    ub = c(NA, NA, NA),
    type = "real",
    log_lik = log_lik_obs_gp,
    posterior_predict = posterior_predict_obs_gp,
    posterior_epred = posterior_epred_obs_gp
  )
  family$stanvars <- brms::stanvar(
    scode = r"(
  // Implemented in C++ in gp.hpp
  real log_nc(real mu_O, real mu_C, real sigma);

  // The likelihood exluding the normalisation constant (log_nc)
  real obs_gp_lpdf(real t, real mu, real duration, real sigma) {
    real mu_O = mu;
    real mu_C = mu_O + duration;
    real ll;
    if(duration <= 0) {
       reject("Duration must be positive but is ", duration, ".");
    }
    if(t > mu_O) {
      // observed tiem after onset, CDF is going to be small (and thus precise).
      // straightforward computation
      ll = log_diff_exp(normal_lcdf(mu_C | t, sigma), normal_lcdf(mu_O | t, sigma));
    } else {
      // observed tiem before onset, CDF might be imprecise.
      // flipping around t to get small values again
      ll = log_diff_exp(normal_lcdf(2 * t - mu_O | t, sigma), normal_lcdf(2 * t - mu_C | t, sigma));
    }
    return(ll - log_nc(mu_O, mu_C, sigma));
  }
)",
    block = "functions"
  )
  
  gp_file <- system.file("stan", "gp.hpp", package = "phenoCollectR")
  family$stan_model_args <- list(
    stanc_options = list("allow-undefined"),
    user_header = gp_file)
  return(family)
}