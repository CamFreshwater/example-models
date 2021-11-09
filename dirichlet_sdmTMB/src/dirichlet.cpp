#include <TMB.hpp>

template<class Type>
Type ddirmultinom(vector<Type> obs, vector<Type> p, Type phi, int do_log) {
  int dim = obs.size();
  Type N = obs.sum();
  // e.g. Eq. 4 Thorson et al. 2017 Fish. Res.; phi here = Beta in paper
  // or Eq. B.2 in https://www.sciencedirect.com/science/article/pii/S0165783621000953#sec0190
  // this is the 'saturating' version
  Type ll = lgamma(N + Type(1.0)) + lgamma(phi) - lgamma(N + phi); // eq. outside of summations
  for (int a = 0; a < dim; a++) {
    ll += -lgamma(obs(a) + Type(1.0)) +
      lgamma(obs(a) + phi * (p(a) + Type(1.0e-15))) - // 1e-15 for robustness to 0s
      lgamma(phi * (p(a) + Type(1.0e-15)));
  }
  if (do_log) return ll;
  else return exp(ll);
}

template<class Type>
  Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(paa_obs);
  DATA_VECTOR(Neff);

  PARAMETER(log_phi);
  PARAMETER_VECTOR(p);

  int n_t = paa_obs.rows();
  Type nll = Type(0.0);

  // paa_pred usually produced by dynamics of assessment model, but here:
  vector<Type> paa_pred = exp(p) / exp(p).sum(); // simplex; sum to 1

  for (int i = 0; i < n_t; i++) {
    vector<Type> temp_n = Neff(i) * paa_obs.row(i);
    nll -= ddirmultinom(temp_n, paa_pred, exp(log_phi), true);
  }
  return nll;
  }
