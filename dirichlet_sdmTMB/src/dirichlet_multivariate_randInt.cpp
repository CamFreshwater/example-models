#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() () 
{
  // load namespace with multivariate distributions
  using namespace density;

  // data inputs
  DATA_MATRIX(y_obs);    // response matrix (a n-by-k matrix)
  DATA_MATRIX(fx_cov);   // covariate matrix (a n-by-j matrix, first column is 1 to account for intercept)
  DATA_IVECTOR(rfac);    // vector of random factor levels
  DATA_INTEGER(n_rfac);  // number of random factor levels
  DATA_MATRIX(pred_cov);    // model matrix for predictions
  // DATA_IVECTOR(pred_factor1k_i); // vector of predicted random intercepts

  // parameter inputs
  PARAMETER_MATRIX(z_ints); // parameter matrix
  PARAMETER_MATRIX(z_rfac);  // matrix of random intercepts (n_rfac x n_cat)
  PARAMETER(log_sigma_rfac); // among random intercept SD

  // The dimensions of data
  int n_obs = y_obs.rows();         // number of observations
  int n_cat = y_obs.cols();         // number of categories
  // int n_levels = pred_cov.rows();   // number of covariates to make predictions on

  // Matrix for intermediate objects
  matrix<Type> total_eff(n_obs, n_cat); // matrix of combined fixed/random eff

  // calculate effects
  matrix<Type> fx_eff = fx_cov * z_ints;

  for (int i = 0; i < n_obs; ++i) {
    for(int k = 0; k < n_cat; k++) {
      total_eff(i, k) = fx_eff(i, k) + z_rfac(rfac(i), k);
    }
  }

  matrix<Type> gamma = exp(total_eff.array()); // add random effect
  vector<Type> n_plus = y_obs.rowwise().sum(); // row sum of response
  vector<Type> gamma_plus = gamma.rowwise().sum(); // row sum of gamma
  
  Type jll = 0; // initialize joint log-likelihood
  for(int i = 0; i <= (n_obs - 1); i++){
    jll = jll + lgamma((n_plus(i) + 1));
    jll = jll + lgamma(gamma_plus(i));
    jll = jll - lgamma((n_plus(i) + gamma_plus(i)));
    for(int k = 0; k <= (n_cat - 1); k++){
      jll += lgamma((y_obs(i, k) + gamma(i, k)));
      jll -= lgamma(gamma(i, k));
      jll -= lgamma((y_obs(i, k) + 1));
    }
  }

  Type jnll;
  jnll = -jll;


  // make covariance matrix
  matrix<Type> cov_mat(n_cat, n_cat);
  for (int j = 0; j < n_cat; j++) {
    for (int jj = 0; jj < n_cat; jj++) {
      if (j == jj) {
        cov_mat(j, jj) = exp(log_sigma_rfac) * exp(log_sigma_rfac);
      } else {
        cov_mat(j, jj) = 0;
      }
    }
  }

  // Probability of multivariate random intercepts
  for (int h = 0; h < n_rfac; h++) {
    vector<Type> z_rfac_vec = z_rfac.row(h);
    // jnll -= density::MVNORM(sigma_mat)(z_rfac_vec);
    MVNORM_t<Type> neg_log_dmvnorm(cov_mat);
    jnll += neg_log_dmvnorm(z_rfac_vec);
  }

  // N_0_Sigma.cov();                   // Returns covariance matrix (Sigma in this case)
  // REPORT(N_0_Sigma.cov());           // Report back to R

  Type sigma_rfac = exp(log_sigma_rfac);
  ADREPORT(sigma_rfac);
  
  // calculate predictions
  // matrix<Type> pred_eff_fe(n_levels, n_cat);    //pred fixed effects on log scale
  // matrix<Type> pred_eff(n_levels, n_cat);    //pred FE + RE on log scale
  // matrix<Type> pred_gamma(n_levels, n_cat);  //transformed pred effects 
  // vector<Type> pred_gamma_plus(n_levels);        
  // vector<Type> pred_theta(n_levels); 
  // matrix<Type> pred_pi(n_levels, n_cat);      // predicted counts in real 
  // vector<Type> pred_n_plus(n_levels); 
  // matrix<Type> pred_pi_prop(n_levels, n_cat); // predicted counts as ppn.

  // pred_eff_fe = pred_cov * z_ints; 

  // for (int i = 0; i < n_obs; ++i) {
  //   for(int k = 0; k < n_cat; k++) {
  //     pred_eff(i, k) = pred_eff_fe(i, k) + z_rfac(pred_factor1k_i(i), k);
  //   }
  // }

  
  // pred_gamma = exp(pred_eff.array());
  // pred_gamma_plus = pred_gamma.rowwise().sum();
  // pred_theta = 1 / (pred_gamma_plus + 1);
  // for(int m = 0; m < n_levels; m++) {
  //   for(int k = 0; k < n_cat; k++) {
  //     pred_pi(m, k) = pred_gamma(m, k) / pred_theta(m);
  //   }
  // }
  // pred_n_plus = pred_pi.rowwise().sum();
  // for(int m = 0; m < n_levels; m++) {
  //   for(int k = 0; k < n_cat; k++) {
  //     pred_pi_prop(m, k) = pred_pi(m, k) / pred_n_plus(m);
  //   }
  // }

  // REPORT(pred_gamma);
  // ADREPORT(pred_gamma);
  // REPORT(pred_pi);
  // ADREPORT(pred_pi);
  // REPORT(pred_pi_prop);
  // ADREPORT(pred_pi_prop);
  
  // Return negative loglikelihood
  return jnll;
  
}
