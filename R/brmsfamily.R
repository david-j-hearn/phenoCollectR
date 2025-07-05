robs_gp <- function(n, mu_O, duration, sigma) {
  res <- rep(NA_real_, n)
  mu_C <- mu_O + duration
  repeat {
    na_indices <- is.na(res)
    n_na <- sum(na_indices)
    if(n_na == 0) {
      break
    }
    t <- runif(n_na)
    gp_t_all <- rnorm(n_na, t, sigma)
    seen <- mu_O < gp_t_all & gp_t_all < mu_C
    t[!seen] <- NA_real_
    res[na_indices] <- t
  }
  while(any(is.na(res))) {
    
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

obs_gp_brmsfamily <- function(link = "identity", link_duration = "log", link_sigma = "log") {
  family <- brms::custom_family(
    "obs_gp",
    dpars = c("mu", "duration", "sigma"),
    links = c(link, link_duration, link_sigma),
    lb = c(NA, 0, 0),
    ub = c(NA, NA, NA),
    type = "real",
    # log_lik = log_lik_obs_gp, # Not implemented so far, but possible
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