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

dobs_gp <- function(t, mu_O, duration, sigma, log = FALSE) {
  mu_C <- mu_O + duration
  base_ll <- 
    dplyr::if_else(t > mu_O, 
                   brms:::log_diff_exp(pnorm(mu_C, t, sigma, log = TRUE), pnorm(mu_O, t, sigma, log = TRUE)),
                   brms:::log_diff_exp(pnorm(2 * t - mu_O, t, sigma, log = TRUE), pnorm(2 * t - mu_C, t, sigma, log = TRUE))
    )
  nc <- obs_gp_nc(mu_O, duration, sigma)                            
  ll <- base_ll - log(nc)
  if(log) {
    ll
  } else {
    exp(ll)
  }
}

robs_gp <- function(n, mu_O, duration, sigma, max_rejection = 10000) {
  if(length(mu_O) == 1) {
    mu_O <- rep(mu_O, n)
  }
  if(length(duration) == 1) {
    duration <- rep(duration, n)
  }
  if(length(sigma) == 1) {
    sigma <- rep(sigma, n)
  }
  
  res <- rep(NA_real_, n)
  
  prob_reject_estimate_time <- pmin(1, 
                                    pmax(0, -mu_O + sigma) / (duration + 2 * sigma) + 
                                      pmax(0, -1 + mu_O + duration + sigma) / (duration + 2 * sigma)) 
  prob_reject_estimate_obs <- pmin(1, 
                                    pmax(0, mu_O - sigma) + 
                                      pmax(0, 1 - mu_O - duration - sigma)) 
  time_indices <- prob_reject_estimate_time <= prob_reject_estimate_obs
    res[time_indices] <- robs_gp_time(sum(time_indices), mu_O[time_indices], duration[time_indices], sigma[time_indices])
  res[!time_indices] <- robs_gp_observation(sum(!time_indices), mu_O[!time_indices], duration[!time_indices], sigma[!time_indices])
  res
}

# Sample observation times and check whether it is within onset/cessation.
# Will be faster choice if onset - cessation has a large window outside of [0, 1]
robs_gp_observation <- function(n, mu_O, duration, sigma, max_rejection = 10000) {
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

# Sample times within onset - cessation and the check whether they are within
# [0, 1].
# Will be faster choice if onset-cessation lies wholly (or at least motly) within [0, 1] 
robs_gp_time <- function(n, mu_O, duration, sigma, max_rejection = 10000) {
  res <- rep(NA_real_, n)
  if(length(mu_O) == 1) {
    mu_O <- rep(mu_O, n)
  }
  if(length(duration) == 1) {
    duration <- rep(duration, n)
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
    onset <- rnorm(n_na, mean = mu_O[na_indices], sd = sigma[na_indices])
    cessation <- onset + duration[na_indices] 
    t <- runif(n_na, onset, cessation)
    seen <- t > 0 & t < 1
    t[!seen] <- NA_real_
    res[na_indices] <- t
  }
  if(any(is.na(res))) {
    warning("Could not produce viable observations for some indices by rejection sampling.")
  }
  res
}

posterior_predict_obs_gp <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  duration <- brms::get_dpar(prep, "duration", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(robs_gp(prep$ndraws, mu, duration, sigma))
}

posterior_epred_obs_gp <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  duration <- brms::get_dpar(prep, "duration")
  return(mu + 0.5 * duration)
}

log_lik_obs_gp <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  duration <- brms::get_dpar(prep, "duration", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  t <- prep$data$Y[i]
  
  dobs_gp(t, mu, duration, sigma, log = TRUE)
}

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