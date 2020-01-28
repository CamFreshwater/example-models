## 7. Estimation of survival probabilities using capture-recapture data
## 7.5. Models with individual variation
## 7.5.2. Random group effects

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data2 <- read_rdump(here::here("BPA", "Ch.07", "cjs_group.data.R"))

## Parameters monitored
params <- c("mean_phi", "mean_p", "phi_g", "sigma", "beta")

## MCMC settings
ni <- 4000
nt <- 1
nb <- 3000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i) {
    list(phi_g = runif(length(unique(stan_data2$group)), 0, 1),
         p_g = runif(length(unique(stan_data2$group)), 0, 1))})

## Call Stan from R
cjs_group_raneff <- stan(here::here("BPA", "Ch.07",
                                    "cjs_group_raneff.stan"),
                         data = stan_data2, init = inits, pars = params,
                         chains = nc, iter = ni, warmup = nb, thin = nt,
                         control = list(adapt_delta = 0.9),
                         seed = 1,
                         open_progress = FALSE)
saveRDS(cjs_group_raneff, here::here("BPA", "Ch.07", "gen_data",
                                     "group_raneff.rds"))

## Summarize posteriors
print(cjs_group_raneff, digits = 3)


# Adjusted model with time-varying detection probabilities
params2 <- c("mean_phi", "mean_p", "phi_g", "sigma", "beta")

inits2 <- lapply(1:nc, function(i) {
  list(gamma_phi = rnorm(stan_data2$n_occasions - 1),
       #gamma_pp = rnorm(stan_data2$n_occasions - 1),
       phi_g = runif(length(unique(stan_data2$group)), 0, 1))
  })

cjs_groupRE_tempFE_mod <- stan_model(here::here("BPA", "Ch.07",
                                  "cjs_group_raneff_temp_fixeff.stan"))
cjs_groupRE_tempFE  <- sampling(cjs_groupRE_tempFE_mod, data = stan_data,
                                init = inits2,
                      pars = params2,
                      chains = nc, iter = ni, warmup = nb, thin = nt,
                      control = list(adapt_delta = 0.9),
                      seed = 1, open_progress = FALSE)
saveRDS(cjs_groupRE_tempFE, here::here("BPA", "Ch.07", "gen_data",
                                     "group_raneff_temp_fixeff.rds"))

print(cjs_groupRE_tempFE, digits = 3)

