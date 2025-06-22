#include <stan/model/model_header.hpp>

template <typename T0__, typename T1__, typename T2__,
          stan::require_all_t<stan::is_stan_scalar<T0__>,
                              stan::is_stan_scalar<T1__>,
                              stan::is_stan_scalar<T2__>>* = nullptr>
  stan::promote_args_t<T0__, T1__, T2__>
  log_nc(const T0__& mu_O, const T1__& mu_C, const T2__& sigma, const int&
  debug, std::ostream* pstream__);

template <>
  double
  log_nc(const double& mu_O, const double& mu_C, const double& sigma, const int&
  debug, std::ostream* pstream__) {
    double mu_C_m1_over_sigma = (mu_C - 1) / sigma;
    double mu_C_over_sigma = mu_C / sigma;
    double mu_O_m1_over_sigma = (mu_O - 1) / sigma;
    double mu_O_over_sigma = mu_O / sigma;

    double exp_sum = -exp(-0.5 * stan::math::square(mu_C_m1_over_sigma)) + exp(-0.5 * stan::math::square(mu_C_over_sigma)) +
      exp(-0.5 * stan::math::square(mu_O_m1_over_sigma)) - exp(-0.5 * stan::math::square(mu_O_over_sigma));
    double scaled_exp_sum = sqrt(1.0/(2 * stan::math::pi())) * sigma * exp_sum;
    double mu_C_terms = -(mu_C - 1) * (stan::math::std_normal_cdf(mu_C_m1_over_sigma) - 0.5) +
      mu_C * (stan::math::std_normal_cdf(mu_C_over_sigma) - 0.5);
    double mu_O_terms = -(mu_O - 1) * stan::math::std_normal_cdf(-mu_O_m1_over_sigma) +
      mu_O * stan::math::std_normal_cdf(-mu_O_over_sigma);

    double res = log(-0.5 + scaled_exp_sum + mu_C_terms + mu_O_terms);

    return res;
}


template <>
stan::math::var
  log_nc(const stan::math::var& mu_O, const stan::math::var& mu_C, const double& sigma, const int&
    debug, std::ostream* pstream__) {

    using namespace std;

    double mu_C_m1_over_sigma = (mu_C.val() - 1) / sigma;
    double mu_C_over_sigma = mu_C.val() / sigma;
    double mu_O_m1_over_sigma = (mu_O.val() - 1) / sigma;
    double mu_O_over_sigma = mu_O.val() / sigma;


    double exp_sum = -exp(-0.5 * stan::math::square(mu_C_m1_over_sigma)) + exp(-0.5 * stan::math::square(mu_C_over_sigma)) +
      exp(-0.5 * stan::math::square(mu_O_m1_over_sigma)) - exp(-0.5 * stan::math::square(mu_O_over_sigma));
    double scaled_exp_sum = sqrt(1.0/(2 * stan::math::pi())) * sigma * exp_sum;
    double mu_C_terms = -(mu_C.val() - 1) * (stan::math::std_normal_cdf(mu_C_m1_over_sigma) - 0.5) +
      mu_C.val() * (stan::math::std_normal_cdf(mu_C_over_sigma) - 0.5);
    double mu_O_terms = -(mu_O.val() - 1) * stan::math::std_normal_cdf(-mu_O_m1_over_sigma) +
      mu_O.val() * stan::math::std_normal_cdf(-mu_O_over_sigma);

    double nc = -0.5 + scaled_exp_sum + mu_C_terms + mu_O_terms;
    double value = log(nc);

    double raw_grad_mu_O =  stan::math::std_normal_cdf(mu_O_m1_over_sigma) - stan::math::std_normal_cdf(mu_O_over_sigma);
    double raw_grad_mu_C = -stan::math::std_normal_cdf(mu_C_m1_over_sigma) + stan::math::std_normal_cdf(mu_C_over_sigma);

    double grad_mu_O = raw_grad_mu_O / nc;
    double grad_mu_C = raw_grad_mu_C / nc;

    if(debug == 1) {
      *pstream__ << "muC overs: " << mu_C_over_sigma << ", " << mu_C_m1_over_sigma <<
        " - muO: " << mu_O_over_sigma << ", " << mu_O_m1_over_sigma << std::endl;
      *pstream__ << "Raw grad: " << raw_grad_mu_O << ", " << raw_grad_mu_C << std::endl;
      *pstream__ << "Grad: " << grad_mu_O << ", " << grad_mu_C << std::endl;
      *pstream__ << "Val: " << value << std::endl;
    }

    std::vector<stan::math::var> operands{ mu_O, mu_C };
    std::vector<double> gradients{ grad_mu_O, grad_mu_C };


    return stan::math::precomputed_gradients(value, operands, gradients);
  }



template <>
stan::math::var
log_nc(const stan::math::var& mu_O, const stan::math::var& mu_C, const stan::math::var& sigma, const int&
  debug, std::ostream* pstream__) {

  using namespace std;
  using stan::math::square;
  using stan::math::std_normal_cdf;

  double mu_C_m1_over_sigma = (mu_C.val() - 1) / sigma.val();
  double mu_C_over_sigma = mu_C.val() / sigma.val();
  double mu_O_m1_over_sigma = (mu_O.val() - 1) / sigma.val();
  double mu_O_over_sigma = mu_O.val() / sigma.val();


  double exp_sum = -exp(-0.5 * square(mu_C_m1_over_sigma)) + exp(-0.5 * square(mu_C_over_sigma)) +
    exp(-0.5 * square(mu_O_m1_over_sigma)) - exp(-0.5 * square(mu_O_over_sigma));
  double scaled_exp_sum = sqrt(1.0/(2 * stan::math::pi())) * sigma.val() * exp_sum;
  double mu_C_terms = -(mu_C.val() - 1) * (std_normal_cdf(mu_C_m1_over_sigma) - 0.5) +
    mu_C.val() * (std_normal_cdf(mu_C_over_sigma) - 0.5);
  double mu_O_terms = -(mu_O.val() - 1) * std_normal_cdf(-mu_O_m1_over_sigma) +
    mu_O.val() * std_normal_cdf(-mu_O_over_sigma);

  double nc = -0.5 + scaled_exp_sum + mu_C_terms + mu_O_terms;
  double value = log(nc);

  // Gradient w.r.t mus
  double raw_grad_mu_O =  std_normal_cdf(mu_O_m1_over_sigma) - std_normal_cdf(mu_O_over_sigma);
  double raw_grad_mu_C = -std_normal_cdf(mu_C_m1_over_sigma) + std_normal_cdf(mu_C_over_sigma);

  //Gradient w.r.t sigma
  double mu_C_sq = square(mu_C.val());
  double mu_O_sq = square(mu_O.val());
  double inv_2_sigma_sq = 1.0 / (2 * square(sigma.val()));
  double disgma_term1 = -exp(-mu_O_sq * inv_2_sigma_sq);
  double disgma_term2 =  exp((-mu_O_sq - 1 + 2 * mu_O.val()) * inv_2_sigma_sq);
  double disgma_term3 =  exp(-mu_C_sq * inv_2_sigma_sq);
  double disgma_term4 = -exp((-mu_C_sq - 1 + 2 * mu_C.val()) * inv_2_sigma_sq);
  double disgma_numerator = disgma_term1 + disgma_term2 + disgma_term3 + disgma_term4;
  double disgma_denominator = sqrt(2 * stan::math::pi());

  double raw_grad_sigma =  disgma_numerator / disgma_denominator;

  double grad_mu_O = raw_grad_mu_O / nc;
  double grad_mu_C = raw_grad_mu_C / nc;
  double grad_sigma = raw_grad_sigma / nc;

  if(debug == 1) {
    *pstream__ << "muC overs: " << mu_C_over_sigma << ", " << mu_C_m1_over_sigma <<
      " - muO: " << mu_O_over_sigma << ", " << mu_O_m1_over_sigma << std::endl;
    *pstream__ << "grad sigma terms: " << disgma_term1 << ", " << disgma_term2 << ", " <<
      disgma_term3 << ", " << disgma_term4 <<
      " denom: " << disgma_denominator << std::endl;

    *pstream__ << "Raw grad: " << raw_grad_mu_O << ", " << raw_grad_mu_C <<  ", " << raw_grad_sigma << std::endl;
    *pstream__ << "Grad: " << grad_mu_O << ", " << grad_mu_C << ", " << grad_sigma << std::endl;
    *pstream__ << "Val: " << value << std::endl;
  }

  std::vector<stan::math::var> operands{ mu_O, mu_C, sigma };
  std::vector<double> gradients{ grad_mu_O, grad_mu_C, grad_sigma };


  return stan::math::precomputed_gradients(value, operands, gradients);
}


template <typename T0__, typename T1__, typename T2__,
          stan::require_all_t<stan::is_stan_scalar<T0__>,
                              stan::is_stan_scalar<T1__>,
                              stan::is_stan_scalar<T2__>>* = nullptr>
  stan::promote_args_t<T0__, T1__, T2__>
  log_nc(const T0__& mu_O, const T1__& mu_C, const T2__& sigma , std::ostream* pstream__)
{
  return(log_nc(mu_O, mu_C, sigma, 0, pstream__));
}

