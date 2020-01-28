## 7. Estimation of survival probabilities using capture-recapture data
## 7.6. Models with time and group effects
## 7.6.2. Fixed group and random time effects

## Note these models are a jumping off point to account for covariance among
# phi and p but as currently developed are not applicable to our tagging work

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump(here::here("BPA", "Ch.07", "cjs_add.data.R"))

## Parameters monitored
params <- c("eta_phi", "p_g", "Sigma", "mean_phi")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i) {
    list(p_g = runif(length(unique(stan_data$group)), 0, 1),
         Omega = matrix(c(1, 0, 0, 1), ncol = 2))})

## Call Stan from R
cjs_temp_corr <- stan(here::here("BPA", "Ch.07", "cjs_temp_corr.stan"),
                      data = stan_data, init = inits, pars = params,
                      chains = nc, iter = ni, warmup = nb, thin = nt,
                      seed = 1,
                      open_progress = FALSE)

## Summarize posteriors
print(cjs_temp_corr, digits = 3)

## Secondary version with random probability effects
inits2 <- lapply(1:nc, function(i) {
  list(p_g = runif(length(unique(stan_data$group)), 0, 1),
       Omega = matrix(c(1, 0, 0, 1), ncol = 2))})
params2 <- c("eta_phi", "eta_p", "Sigma_phi", "Sigma_p", "mean_phi", "mean_p")

cjs_temp_corr2 <- stan_model(here::here("BPA", "Ch.07", "cjs_temp_corr2.stan"))
fit_temp_corr2  <- sampling(cjs_temp_corr2, data = stan_data,
                      #init = inits,
                      pars = params2,
                      chains = nc, iter = ni, warmup = nb, thin = nt,
                      seed = 1, open_progress = FALSE)
print(fit_temp_corr2, digits = 3)
saveRDS(fit_temp_corr2, here::here("BPA", "Ch.07", "gen_data", "hier2.rds"))


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

phi_re <- rstan::extract(fit_temp_corr2)$eta_phi %>%
  boot::inv.logit()
# calculate cumulative survival
cum_surv1 <- t(apply(phi_re[, , 1], 1, function(x) cumprod(x)))
cum_surv2 <- t(apply(phi_re[, , 2], 1, function(x) cumprod(x)))

dat1 <- gen_int(cum_surv1, parm = "phi", group = "1")
dat2 <- gen_int(cum_surv2, parm = "phi", group = "2")
plot_dat <- rbind(dat1, dat2)

ggplot(plot_dat, aes(x = as.factor(stage), y = mu)) +
  geom_pointrange(aes(ymin = low, ymax = high)) +
  facet_wrap(~as.factor(group))
