library(tidyverse)

set.seed(42)


# random intercepts for groups (e.g. sites)
n_sites <- 5
n_obs_per_site <- 100
sd_site <- 0.5
n <- n_sites * n_obs_per_site
site_dat <- data.frame(site = seq(1, n_sites, 1),
                       site_mean = rnorm(mean = 0, sd = sd_site,
                                         n = n_sites))

J <- 3 #n categories (e.g. stocks)
P <- 2 #n categorical fixed effects (or covariates - 1 if some continuous)
N <- sample(c(30:70), n, replace = TRUE) #sample size per observation

# input data frame and design matrix
dat <- data.frame(strata = sample(1:P, n, replace = T)) %>%
  mutate(strata_f = as.factor(strata),
         site = sample(1:n_sites, n, replace = T),

         site_f = as.factor(site),
         N = N)

# fixed effects
X <- model.matrix(~strata_f + 0, dat)
# generate matrix of coefficients has fixed covariates (rows) and
# categories (cols)
beta0 <- matrix(rnorm((ncol(X)) * J), ncol(X), J)

fix_eff <- X %*% beta0
colnames(fix_eff) <- paste("beta_j", seq(1, J, by = 1), sep = "")

# combine fixed and random effects
dat2 <- dat %>%
  left_join(., site_dat, by = "site") %>%
  cbind(., fix_eff)


# generate data
Gamma = exp(fix_eff + dat2$site_mean) #fixed effects
Gamma_plus = apply(Gamma, 1, sum) #sum of fixed_effects
theta = 1 / (Gamma_plus + 1)
pi = apply(Gamma, 2, function(x) {x / theta})

pi_plus <- apply(pi, 1, sum)
pred_pi <- pi / pi_plus


# Gamma and pi vary among random intercepts
unique(Gamma)
unique(pi)

# predicted proportions (i.e. pred_pi) converge on ~same values despite different
# random intercepts
unique(pred_pi)
