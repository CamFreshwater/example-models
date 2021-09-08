READ ME

Directory includes scripts and data for developing a negative binomial / Dirichlet-multinomial model that includes a Gaussian random fields component

- NB/Dirichlet from stockseasonr
- GRF from gfranges repo (https://github.com/pbs-assess/gfranges/blob/421e4af3f649865175ef9aa0953a12c7b50255a1/analysis/VOCC/vocc_regression.cpp#L85-L101)

- Goal is to estimate spatiotemporal variability in abundance as well as composition




  // Matern:
  Type range = sqrt(Type(8.0)) / exp(ln_kappa);
  vector<Type> sigma_O(n_k);
  for(int k = 0; k < n_k; k++) {
    sigma_O(k) = 1 / sqrt(Type(4.0) * M_PI * exp(Type(2.0) * ln_tau_O(k)) *
      exp(Type(2.0) * ln_kappa));
  }
  Eigen::SparseMatrix<Type> Q; // Precision matrix
  Q = R_inla::Q_spde(spde, exp(ln_kappa));

  // ------------------ INLA projections ---------------------------------------

  // Here we are projecting the spatial random effects to the
  // locations of the data using the INLA 'A' matrices.
  array<Type> omega_sk_A(A_sk.rows(), A_sk.cols());
  for (int k = 0; k < n_k; k++)
    omega_sk_A.col(k) = A_sk * vector<Type>(omega_sk.col(k));
  vector<Type> omega_sk_A_vec(n_i);
