## CJS Model Structure Comparison
# Fit two approximately equivalent CJS models to shared data to evaluate
# which structure is more appropriate.
# Model 1 (blmeco_cjs) is based on swallows example model from Korner-Nievergelt
# et al. 2015. Model is less versatile, but more intuitive and clear how p
# could be fixed for time steps of choice.
# Model 2 (stan_ref_cjs) is based on the individal CJS models from Ch. 7 of Stan
# reference guide.

library(blmeco)
library(rstan)
library(tidyverse)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)


stan_data <- read_rdump(here::here("BPA", "Ch.07", "cjs_add.data.R"))

## MCMC settings
ni <- 1000
nt <- 1
nb <- 250
nc <- 4


# Model 1 - blmeco
blmeco_data <- stan_data
names(blmeco_data) <- c("CH", "K", "ngroup", "group", "I", "df", "R")

inits_blmeco <- lapply(1:nc, function(i) {
  list(#a = rnorm(stan_data$n_occasions - 1),
       #beta = c(0, rnorm(1)),
       phi_g = runif(length(unique(stan_data$group)), 0, 1),
       p_g = runif(length(unique(stan_data$group)), 0, 1))})

blmeco_params <- c("phi_g", "p_g") #c("a", "beta", "p_g")
blmeco_cjs_mod <- here::here("cjs_comparison", "blmeco_cjs.stan")
blmeco_cjs  <- stan(blmeco_cjs_mod,
                    data = blmeco_data, init = inits_blmeco[1],
                    pars = blmeco_params,
                    chains = 1, iter = 100, warmup = 25, thin = nt,
                    seed = 1,
                    open_progress = FALSE)


# Model 2 - stan reference
inits_stan <- lapply(1:nc, function(i) {
  list(gamma = rnorm(stan_data$n_occasions - 1),
       beta = c(0, rnorm(1)),
       p_g = runif(length(unique(stan_data$group)), 0, 1))})

stan_params <- c("phi_g1", "phi_g2", "p_g", "beta")

stan_cjs_mod <- here::here("cjs_comparison", "stan_ref_cjs.stan")
stan_cjs  <- stan(stan_cjs_mod,
                  data = stan_data, init = inits_stan[1], pars = stan_params,
                  # chains = nc, iter = ni, warmup = nb, thin = nt,
                  chains = 1, iter = 100, warmup = 25, thin = nt,
                  seed = 1,
                  open_progress = FALSE)



#### Swallows example

data("survival_swallows")
dat <- survival_swallows

blmeco_cjs_mod1 <- here::here("cjs_comparison", "blmeco_cjs_original.stan")
mod <- stan(blmeco_cjs_mod1, data = dat, chains = 1, iter = 1000)

