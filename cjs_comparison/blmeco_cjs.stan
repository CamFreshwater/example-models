// Stan Code for CJS Model w/ Random Effects from Bayesian Data
// Analysis in Ecology Using Linear Models with R, BUGS, and Stan
// Modified to match stan reference example

data {
	int<lower = 2> K; // capture events
	int<lower = 0> I; // number of individuals
	int<lower = 0, upper = 1> CH[I, K]; // CH[i,k]: individual i captured at K
	int<lower = 0> ngroup; // number of groups
	int<lower = 0, upper = ngroup> group[I];// index of group variable
}

transformed data {
	// real beta1 = 0;      // Corner constraint
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
	//real a[K - 1];			// recapture-specific intercept of phi
    vector<lower=0,upper=1>[ngroup] p_g;    // Group-spec. recapture
    vector<lower=0,upper=1>[ngroup] phi_g;    // Group-spec. surv
	//real beta2; 			// difference in male/female survival
}

transformed parameters {
	real<lower = 0, upper = 1> p[I, K]; // capture probability
	real<lower = 0, upper = 1> phi[I, K - 1]; // survival probability
	real<lower = 0, upper = 1> chi[I, K + 1]; // probability that an ndividual  
									   // is never recap. after its last cap.

	//vector[ngroup] beta;
	//beta[1] = beta1;
  	//beta[2] = beta2;
	
	{
		int k;
		for(ii in 1:I) {
			for(tt in 1:(K - 1)) {
				phi[ii,tt] = phi_g[group[ii]]; // inv_logit(a[tt] + beta[group[ii]]);
			}
		}
		for(i in 1:I) {
			p[i,1] = 1; // first occasion is marking occasion
			for(kk in 2:K) {
				p[i, kk] = p_g[group[i]];
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
	// a ~ normal(0, 10);
	// beta2 ~ normal(0, 10);
	phi_g ~ uniform(0, 1);
	p_g ~ uniform(0, 1);
	
	// likelihood
	for (i in 1:I) {
		if (last[i] > 0) {
			for (k in 1:last[i]) {
				if(k > 1) 1 ~ bernoulli(phi[i, k-1]);
				CH[i,k] ~ bernoulli(p[i, k]);
			}
		}
		target += log(chi[i, last[i] + 1]);
	}
}
