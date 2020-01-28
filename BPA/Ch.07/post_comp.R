## Compare posteriors generated with other models

require(rstan)
require(tidyverse)

# additive fixed effects (cjs_add2.stan); only model to include time-varying p
fe_fit <- readRDS(here::here("BPA", "Ch.07", "gen_data", "add2.rds"))
# random group effects (cjs_group_raneff.stan)
re_fit <- readRDS(here::here("BPA", "Ch.07", "gen_data", "group_raneff.rds"))
# random group effects + temporal fixed (cjs_group_raneff_temp_fixeff.stan)
mix_fit <- readRDS(here::here("BPA", "Ch.07", "gen_data",
                              "group_raneff_temp_fixeff.rds"))

# Helper function generate posterior intervals
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


## Compare average group specific survival
# Extract phis and calculate mean where necessary
make_dat <- function(mcmc_est, model, group, take_mean = TRUE) {
  if (take_mean == TRUE) {
    dat <- mcmc_est %>%
      apply(., 1, mean)
  } else {
    dat <- mcmc_est
  }
  data.frame(trial = seq(from = 1, to = length(dat), by = 1),
             est = dat,
             model = model,
             group = group)
}

phi1_fe <- make_dat(rstan::extract(fe_fit)$phi_g1, model = "fix", group = "1")
phi2_fe <- make_dat(rstan::extract(fe_fit)$phi_g2, model = "fix", group = "2")
phi1_re <- make_dat(rstan::extract(re_fit)$phi_g[ , 1], model = "rand",
                    group = "1", take_mean = FALSE)
phi2_re <- make_dat(rstan::extract(re_fit)$phi_g[ , 2], model = "rand",
                    group = "2", take_mean = FALSE)
phi1_mix <- make_dat(rstan::extract(mix_fit)$phi_g[ , 1, ], model = "mix",
                     group = "1")
phi2_mix <- make_dat(rstan::extract(mix_fit)$phi_g[ , 2, ], model = "mix",
                     group = "2")
post_list <- list(phi1_fe, phi2_fe, phi1_re, phi2_re, phi1_mix, phi2_mix)

plot_dat <- do.call(rbind, post_list) %>%
  group_by(model, group) %>%
  summarize(mean_phi = mean(est),
            low = quantile(est, 0.05),
            high = quantile(est, 0.95))

ggplot(plot_dat, aes(x = model, y = mean_phi)) +
  geom_pointrange(aes(ymin = low, ymax = high)) +
  facet_wrap(~as.factor(group))


tt <- rstan::extract(mix_fit)$phi_g[ , 2, ]
t1 <- rstan::extract(mix_fit)$phi_g[ , 1, ]
