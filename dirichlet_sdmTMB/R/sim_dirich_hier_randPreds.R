# Simulate and fit unordered dirichlet response models with addition of
# random intercepts; based on models from following repo:
# https://github.com/carriegu/STAT-840/tree/master/DMRegression
# Identical to sim_dirich_hier.R but model includes predictions for various RE
# levels

library(dirmult)
library(TMB)
library(tidyverse)
library(DirichletReg)

set.seed(42)

# random intercepts for groups (e.g. sites)
n_sites <- 5
n_obs_per_site <- 100
sd_site <- 0.5
n <- n_sites * n_obs_per_site
site_seq <- seq(1, n_sites, 1)
site_dat <- data.frame(site = site_seq,
                       site_mean = rnorm(mean = 0, sd = sd_site,
                                         n = n_sites))

# fixed covariate data
J <- 4 #n categories (e.g. stocks)
P <- 3 #n categorical fixed effects (or covariates - 1 if some continuous)
N <- sample(c(5:20), n, replace = TRUE) #sample size per observation

# input data frame and design matrix
dat <- data.frame(strata = sample(1:P, n, replace = T)) %>%
  mutate(strata_f = as.factor(strata),
         site = sample(1:n_sites, n, replace = T),
         site_f = as.factor(site),
         N = N)
# fixed effects
X <- model.matrix(~strata_f + 0, dat)
# matrix of coefficients has fixed covariates (rows) and categories (cols)
beta0 <- matrix(rnorm((ncol(X)) * J), ncol(X), J)
# model matrix for predictions
pred_dat <- dat %>%
  arrange(strata_f) %>%
  select(strata_f) %>%
  distinct()
pred_cov <- model.matrix(~strata_f + 0, pred_dat)
pred_dat_long <- expand.grid(strata_f = pred_dat$strata_f,
                             site = site_seq)



# function to simulate dirichlet data based on parameters and matrix
f_sim <- function(trial = 1) {

  # simulate initial beta matrix to generate fixed effects
  fix_eff <- X %*% beta0
  colnames(fix_eff) <- paste("beta_j", seq(1, J, by = 1), sep = "")

  # combine fixed and random effects
  dat2 <- dat %>%
    left_join(., site_dat, by = "site") %>%
    cbind(., fix_eff)

  # reparametrization for data generation
  Gamma = exp(fix_eff + dat2$site_mean) #fixed effects
  Gamma_plus = apply(Gamma, 1, sum) #sum of fixed_effects
  theta = 1/(Gamma_plus + 1)
  pi = apply(Gamma, 2, function(x) {x / theta})

  # Simulate the responses Y from package 'dirmult'
  # set.seed(123)
  Y = simPop(J = 1, n = N[1], pi = pi[1,], theta = theta[1])$data
  for(jj in c(2:n)){
    Y = rbind(Y,
              simPop(J = 1, n = N[jj], pi = pi[jj,], theta = theta[jj])$data)
  }

  # Switch Y to non-integer
  add_noise <- function(nn) {
    out <- nn
    keep <- which(nn > 1)
    xx <- rnorm(length(nn[keep]), 0, 0.2)
    noise = xx - mean(xx)
    out[keep] = out[keep] + noise
    return(out)
  }
  Y2 = matrix(NA, nrow = nrow(Y), ncol = ncol(Y))
  for (i in seq_len(nrow(Y))) {
    Y2[i, ] <- add_noise(Y[i, ])
  }

  # create the observations:
  dat3 <- cbind(dat2, Y)

  return(list("trial" = trial, "obs" = Y, "noisey_obs" = Y2,
              "fixed_fac" = dat3$strata_f,
              "rand_mean" = site_dat$site_mean, "rand_fac" = dat3$site_f,
              "full_data" = dat3, "fixed_ppns" = pi
  ))
}

# simulate
n_trials <- 2
sims <- vector(mode = "list", length = n_trials)
for (i in 1:n_trials) {
  sims[[i]] <- f_sim(trial = i)
}
sims_in <- sims[[1]]

## look at raw data (NEEDS TO BE COMPLETED)
# glimpse(sims[[1]]$full_data)
# dum <- sims[[1]]$full_data %>%
#   pivot_longer(cols = `1`:`4`, names_to = "group", values_to = "count") %>%
#   group_by(strata_f, site_f, N) %>%
#   summarize(obs_ppn = sum(count) / N)
# mean_dum <- dum %>%
#   group_by(strata_f, site_f, group) %>%
#   summarize(mean_obs_ppn = mean(obs_ppn))
#
# ggplot(dum) +
#   geom_bar(aes(x = strata_f, y = count, fill = group),
#            stat = "identity") +
#   ggsidekick::theme_sleek() +
#   facet_wrap(~site_f)


# pred_eff <- pred_cov %*% beta0
# z_rfac <- site_dat$site_mean
# pred_factor2k_i <- re_preds
# pred_dat_re <- do.call(rbind, replicate(5, pred_cov, simplify = FALSE))
# pred_dat_re <- pred_dat_re %*% beta0
# pred_eff_re <- matrix(NA, nrow = nrow(pred_dat_re), ncol = J)

# for (p in 1:(n_sites * P)) { # RE x FE levels
#   for (k in 1:J) { #categories
#     pred_eff_re[p, k] = pred_dat_re[p, k] + z_rfac[re_preds[p]]
#   }
# }


# compile models
compile(here::here("dirichlet_sdmTMB", "src",
                   "dirichlet_randInt_randPreds.cpp"))
dyn.load(dynlib(here::here("dirichlet_sdmTMB", "src",
                           "dirichlet_randInt_randPreds")))


#fit models with all data
fit_list_hier <- map(sims, function(sims_in) {
  Y_in <- sims_in$obs
  rfac <- as.numeric(sims_in$rand_fac) - 1 #subtract 1 for indexing in c++

  # variables for estimating random effect predictions
  n_rfac <- length(unique(rfac))
  pred_dat_re <- do.call(rbind, replicate(5, pred_cov, simplify = FALSE))
  re_preds <- rep(seq(0, max(rfac), by = 1), each = nrow(pred_cov))

  #initial parameter values
  beta_in <- matrix(rnorm((ncol(X)) * J), ncol(X), J)
  rand_int_in <- rep(0, times = length(unique(rfac)))

  if(length(re_preds) != nrow(pred_dat_re)) {
    warning("Random effects misspecified")
  }

  obj <- MakeADFun(
    data = list(fx_cov = X,
                y_obs = Y_in,
                pred_cov = pred_dat_re,
                rfac = rfac,
                n_rfac = n_rfac,
                pred_factor2k_i = re_preds,
                n_levels = length(re_preds),
                re_predictions = 1
    ),
    parameters = list(z_ints = beta_in,
                      z_rfac = rand_int_in,
                      log_sigma_rfac = 0
    ),
    random = c("z_rfac"),
    DLL = "dirichlet_randInt_randPreds"
  )

  opt <- nlminb(obj$par, obj$fn, obj$gr)
  sdr <- sdreport(obj)
  ssdr <- summary(sdr)

  re_preds_out <- ssdr[rownames(ssdr) %in% "pred_pi_prop", ]
  re_preds_dat <- data.frame(pred_ppn = re_preds_out[, "Estimate"]) %>%
    cbind(., do.call(rbind, replicate(J, pred_dat_re, simplify = FALSE))) %>%
    mutate(site = rep(re_preds, times = J),
           group = rep(seq(1, J, by = 1), each = length(re_preds)))


  fix_eff <- data.frame(
    par = "beta",
    par_n = seq(1, length(beta0), by = 1),
    est = ssdr[rownames(ssdr) %in% "z_ints" , 1],
    true = as.vector(beta0)
  )
  rand_eff <- data.frame(
    par = c(rep("alpha", n_rfac), "sigma_a"),
    par_n = c(seq(1, n_rfac, by = 1), 1),
    est = c(ssdr[rownames(ssdr) %in% "z_rfac" , 1],
            ssdr[rownames(ssdr) %in% "sigma_rfac" , 1]),
    true = c(unique(sims_in$rand_mean), sd_site)
  )

  rbind(fix_eff, rand_eff) %>%
    mutate(tran = sims_in$trans,
           trial = sims_in$trial)

  # list(pred_eff = dum, ssdr = ssdr)
})

