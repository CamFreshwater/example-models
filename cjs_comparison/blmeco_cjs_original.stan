// Stan Code for CJS Model w/ Random Effects from Bayesian Data
// Analysis in Ecology Using Linear Models with R, BUGS, and Stan

data {
	int<lower = 2> K; // capture events
	int<lower = 0> I; // number of individuals
	int<lower = 0, upper = 1> CH[I, K]; // CH[i,k]: individual i captured at K
	int<lower = 0> nfam; // number of families
	int<lower = 0, upper = nfam> family[I];// index of group variable
	vector[I] carez; // duration of parental care, z-trans.
	int<lower = 1,upper = 4> year[I]; // index of year
	vector[K] agec; // age of fledgling, centered
}

transformed data {
	int<lower=0,upper=K+1> last[I]; // last[i]: ind i last capture
	last = rep_array(0,I);
	for (i in 1:I) {
		for (k in 1:K) {
			if (CH[i,k] == 1) {
				if (k > last[i]) last[i] = k;
			}
		}
	}
}

parameters {
	real b0[4]; // intercepts per year for p
	real b1[4]; // slope for age per year for p
	real a[K - 1]; // intercept of phi
	real a1; // coef of phi
	real<lower = 0> sigmaphi; // between family sd in logit(phi)
	real<lower = 0> sigmayearphi; // between-year sd in logit(phi)
	real<lower = 0> sigmap; // between family sd in logit(p)
	real fameffphi[nfam]; // family effects for phi
	real fameffp[nfam]; // family effects for p
	real yeareffphi[4]; // year effect on phi
}

transformed parameters {
	real<lower = 0, upper = 1> p[I, K]; // capture probability
	real<lower = 0, upper = 1> phi[I, K - 1]; // survival probability
	real<lower = 0, upper = 1> chi[I, K + 1]; // probability that an ndividual  
									   // is never recap. after its last cap.
	{
		int k;
		for(ii in 1:I) {
			for(tt in 1:(K - 1)) {
				// linear predictor with random effect for phi:
				// add fixed and random effects here
				phi[ii,tt] = inv_logit(a[tt] + (a1 * carez[ii]) + 
					(sigmayearphi * yeareffphi[year[ii]]) + (sigmaphi * 
					fameffphi[family[ii]]));
			}
		}
		for(i in 1:I) {
			// linear predictor with random effect
			// for p: add fixed and random effects here
			p[i,1] = 1; // first occasion is marking occasion
			for(kk in 2:K) {
				p[i, kk] = inv_logit(b0[year[i]] + (b1[year[i]] * agec[kk]) + 
					(sigmap * fameffp[family[i]]));
			}

			// probability that an individual is never recaptured after its
			// last capture
			chi[i, K + 1] = 1.0;
			k = K;
			while (k > 1) {
				chi[i, k] = (1 - phi[i, k-1]) + phi[i, k-1] * (1 - p[i, k]) * 
					chi[i, k+1];
				k = k - 1;
			}
			chi[i,1] = (1 - p[i, 1]) * chi[i, 2];
		}
	}
}

model {
	// priors
	for(j in 1:4) {
		b0[j] ~ normal(0, 5);
		b1[j] ~ normal(0, 5);
	}
	for(v in 1:(K-1)) {
		a[v] ~ normal(0, 5);
	}
	a1 ~ normal(0, 5);
	sigmaphi ~ student_t(2, 0, 1);
	sigmayearphi ~ student_t(2, 0, 1);
	sigmap ~ student_t(2, 0, 1);
	
	// random effects
	for(g in 1:nfam) {
		fameffphi[g] ~ normal(0, 1);
		fameffp[g] ~ normal(0,1);
	}
	for(ye in 1:4) {
		yeareffphi[ye] ~ normal(0, 1);
	}

	// likelihood
	for (i in 1:I) {
		if (last[i] > 0) {
			for (k in 1:last[i]) {
				if(k > 1) 1 ~ bernoulli(phi[i, k-1]);
				CH[i,k] ~ bernoulli(p[i, k]);
			}
		}
		// deprecated function; replace w/ below 
		// increment_log_prob(log(chi[i, last[i] + 1]));
		target += log(chi[i, last[i] + 1]);
	}
}
