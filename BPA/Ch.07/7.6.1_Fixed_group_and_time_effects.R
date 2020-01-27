## 7. Estimation of survival probabilities using capture-recapture data
## 7.6. Models with time and group effects
## 7.6.1. Fixed group and time effects

## Additions by CF include a second time-varying estimate for detection
# probability p

library(rstan)
library(tidyverse)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump(here::here("BPA", "Ch.07", "cjs_add.data.R"))

## Parameters monitored
params <- c("phi_g1", "phi_g2", "p_g", "beta")

## MCMC settings
ni <- 3000
nt <- 2
nb <- 1000
nc <- 4

## Initial values
inits1 <- lapply(1:nc, function(i) {
    list(gamma = rnorm(stan_data$n_occasions - 1),
         beta = c(0, rnorm(1)),
         p_g = runif(length(unique(stan_data$group)), 0, 1))})

## Call Stan from R
cjs_add  <- stan(here::here("BPA", "Ch.07", "cjs_add.stan"),
                 data = stan_data, init = inits1, pars = params,
                 chains = nc, iter = ni, warmup = nb, thin = nt,
                 seed = 1,
                 open_progress = FALSE)
cjs_add_trim  <- stan(here::here("BPA", "Ch.07", "cjs_add_trim.stan"),
                 data = stan_data, init = inits1, pars = params,
                 chains = nc, iter = ni, warmup = nb, thin = nt,
                 seed = 1,
                 open_progress = FALSE)
print(cjs_add, digits = 3)
print(cjs_add_trim, digits = 3)


# Adjusted model with time-varying detection probabilities
params2 <- c("phi_g1", "phi_g2", "p_g1", "p_g2", "beta_phi", "beta_p")

inits <- lapply(1:nc, function(i) {
  list(gamma_phi = rnorm(stan_data$n_occasions - 1),
       gamma_pp = rnorm(stan_data$n_occasions - 1),
       beta_phi = c(0, rnorm(1)),
       beta_p = c(0, rnorm(1)))})

addX_mod <- stan_model(here::here("BPA", "Ch.07", "cjs_add2.stan"))
cjs_add2  <- sampling(addX_mod, data = stan_data, init = inits,
                      pars = params2,
                      chains = nc, iter = ni, warmup = nb, thin = nt,
                      seed = 1, open_progress = FALSE)


## Summarize posteriors
print(cjs_add2, digits = 3)

# Extract phi estimates
gen_int <- function(mcmc_est, parm, group) {
  mu_phi <- apply(mcmc_est, 2, mean)
  low_int <- apply(mcmc_est, 2, function(x) quantile(x, 0.05))
  up_int <- apply(mcmc_est, 2, function(x) quantile(x, 0.95))
  data.frame(stage = seq(1, length(mu_phi), by = 1),
             mu = mu_phi,
             low = low_int,
             high = up_int,
             par = parm,
             group = group)
}

phi1 <- rstan::extract(cjs_add)$phi_g1
dat1 <- gen_int(phi1, parm = "phi", group = "1")
phi2 <- rstan::extract(cjs_add)$phi_g2
dat2 <- gen_int(phi2, parm = "phi", group = "2")
plot_dat <- rbind(dat1, dat2)

ggplot(plot_dat, aes(x = as.factor(group), y = mu)) +
  geom_pointrange(aes(ymin = low, ymax = high)) +
  facet_wrap(~as.factor(stage))

phi1b <- rstan::extract(cjs_add2)$phi_g1
dat1b <- gen_int(phi1b, parm = "phi", group = "1")
phi2b <- rstan::extract(cjs_add2)$phi_g2
dat2b <- gen_int(phi2b, parm = "phi", group = "2")
plot_datb <- rbind(dat1b, dat2b)

ggplot(plot_datb, aes(x = as.factor(group), y = mu)) +
  geom_pointrange(aes(ymin = low, ymax = high)) +
  facet_wrap(~as.factor(stage))
