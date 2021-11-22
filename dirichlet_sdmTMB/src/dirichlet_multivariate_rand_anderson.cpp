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
  // load namespace with multivariate distributions
  using namespace density;

  DATA_MATRIX(paa_obs);
  DATA_VECTOR(Neff);
  // DATA_IVECTOR(rfac);    // vector of random factor levels
  // DATA_INTEGER(n_rfac);  // number of random factor levels
  
  PARAMETER(log_phi);
  PARAMETER_VECTOR(p);
  // PARAMETER_MATRIX(z_rfac);  // matrix of random intercepts (n_rfac x n_cat)
  // PARAMETER(log_sigma_rfac); // among random intercept SD

  // intermediates
  int n_t = paa_obs.rows();
  int n_cat = paa_obs.cols();

  Type nll = Type(0.0);

  vector<Type> paa_pred = exp(p) / exp(p).sum(); // simplex; sum to 1

  for (int i = 0; i < n_t; i++) {
    // adj p_pred to include random effects
    // vector<Type> adj_p(n_cat);
    // for(int k = 0; k < n_cat; k++) {
    //   adj_p(k) = p(k) + z_rfac(rfac(i), k);
    // }
    
    // vector<Type> paa_pred = exp(adj_p) / exp(adj_p).sum(); // simplex; sum to 1
    vector<Type> temp_n = Neff(i) * paa_obs.row(i);

    nll -= ddirmultinom(temp_n, paa_pred, exp(log_phi), true);
  }

  // make covariance matrix
  matrix<Type> sigma_mat(n_cat, n_cat);
  for (int j = 0; j < n_cat; j++) {
    for (int jj = 0; jj < n_cat; jj++) {
      if (j == jj) {
        sigma_mat(j, jj) = 1;
      } else {
        sigma_mat(j, jj) = 0;
      }
    }
  }

  // sigma_mat.fill(0);
  // sigma_mat.diagonal(1); //+= exp(log_sigma_rfac);

  REPORT(sigma_mat);
  ADREPORT(sigma_mat);

  // Probability of multivariate random intercepts
  // for (int h = 0; h < n_rfac; h++) {
  //   vector<Type> z_rfac_vec = z_rfac.row(h);
  //   nll -= density::MVNORM(sigma_mat)(z_rfac_vec);
  // }

  // Probability of random intercepts
  // for (int h = 0; h < n_rfac; h++) {
  //   nll -= dnorm(z_rfac(h), Type(0.0), exp(log_sigma_rfac), true);
  // }


  return nll;
  }
