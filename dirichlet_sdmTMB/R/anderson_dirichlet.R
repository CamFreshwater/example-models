library(tidyverse)
library(TMB)

set.seed(10293)

# functions
rdirichlet <- function(p, phi) {
  alpha <- p * phi + 1e-7
  obs <- rgamma(length(p), alpha, 1)
  obs <- obs / sum(obs)
  obs
}

rdirichlet_linear <- function(p, theta, N) {
  alpha <- p * theta * N + 1e-7
  obs <- rgamma(length(p), alpha, 1)
  obs <- obs / sum(obs)
  obs
}

# rdirichlet_linear2 <- function(p_n, theta) {
#   alpha <- p_n * theta + 1e-7
#   obs <- rgamma(length(p), alpha, 1)
#   obs <- obs / sum(obs)
#   obs
# }

rdirmultinom <- function(N, p, phi, type = "saturating") {
  if (type[[1]] != "saturating") {
    dp <- rdirichlet_linear(p, phi, N)
  } else {
    dp <- rdirichlet(p, phi) # note bug in WHAM simulate(); should be done once per multinomial
  }
  obs <- rmultinom(n = 1, prob = dp, size = N)
  obs
}


# fixed effects
n_cat <- 3
log_pp <- rnorm(n_cat, 0, 1) #runif(50)
phi <- 100      # overdispersion parameter
N <- 200        # 200 samples per observation
n_obs <- 40

# original pars
pp_old <- runif(n_cat)
p_old <- pp_old / sum(pp_old)


## univariate random effects
n_r_ints <- 5
r_sigma <- 0.5
r_ints <- rnorm(n_r_ints, 0, r_sigma)


paa_uni <- purrr::map(r_ints, function (xx) {
  # pp <- exp(pp_old + xx)
  # p <- pp / sum(pp) # 3 categories
  #
  # # adj_phi <- exp(log(phi) + xx) # make phi a function of rand int. (IGNORE)
  # set.seed(123)
  # naa <- matrix(ncol = n_cat, nrow = n_obs, data = 0) # 40 observations
  # for (i in 1:nrow(naa)) naa[i, ] <- rdirmultinom(N, p, phi = phi)
  #
  # paa_uni <- naa / N

  set.seed(123)
  naa <- matrix(ncol = n_cat, nrow = n_obs, data = 0) # 40 observations
  for (i in 1:nrow(naa)) naa[i, ] <- rdirmultinom(N, p_old, phi = phi)

  paa_old <- naa / N
}) %>%
  do.call(rbind, .)

# compare pp and p
# purrr::map(r_ints, function (xx) {
#    pp <- exp(log_pp + xx)
#    p <- pp / sum(pp)
#    # adj_phi <- exp(log(phi) + xx)
#    list(pp, p)
# })


# multivariate random effects (i.e. n_r_ints x n_cat matrix of random intercepts)
# v1 - common variance, zero covariance (controlled by phi)
r_sig_mat <- diag(x = r_sigma, nrow = n_cat, ncol = n_cat)
r_int_mat <- MASS::mvrnorm(n_r_ints, rep(0, n_cat), r_sig_mat)
# split into list for purrr
r_int_list <- as.list(data.frame(t(r_int_mat)))

paa_multi <- purrr::map(r_int_list, function (xx) {
  pp <- exp(log_pp + xx)
  p <- pp / sum(pp) # 3 categories

  set.seed(123)
  naa <- matrix(ncol = n_cat, nrow = n_obs, data = 0) # 40 observations
  for (i in 1:nrow(naa)) naa[i, ] <- rdirmultinom(N, p, phi = phi)

  naa / N
}) %>%
  do.call(rbind, .)

purrr::map(r_int_list, function (xx) {
   pp <- exp(log_pp + xx)
   p <- pp / sum(pp)
   # adj_phi <- exp(log(phi) + xx)
   list(pp, p)
})

r_fac <- rep(seq(0, n_r_ints - 1, by = 1), each = n_obs)


# check data
widen_foo <- function(dat) {
  cbind(dat, r_fac) %>%
    as.data.frame() %>%
    pivot_longer(cols = g_1:g_3, names_to = "group")
}
colnames(paa_uni) <- paste("g", seq(1, n_cat, by = 1), sep = "_")
paa_uni2 <- widen_foo(paa_uni) %>% mutate(dat = "uni")
colnames(paa_multi) <- paste("g", seq(1, n_cat, by = 1), sep = "_")
paa_multi2 <- widen_foo(paa_multi) %>% mutate(dat = "multi")

rbind(paa_uni2, paa_multi2) %>%
  ggplot(.) +
  geom_boxplot(aes(x = group, y = value, fill = dat)) +
  facet_wrap(~r_fac)
# more variability with multivariate REs


## compile
compile(here::here("dirichlet_sdmTMB", "src", "dirichlet.cpp"))
dyn.load(dynlib(here::here("dirichlet_sdmTMB", "src", "dirichlet")))
compile(
  here::here("dirichlet_sdmTMB", "src", "dirichlet_multivariate_rand_anderson.cpp")
)
dyn.load(dynlib(
  here::here("dirichlet_sdmTMB", "src", "dirichlet_multivariate_rand_anderson"))
)


data <- list(paa_obs = paa_old, Neff = rep(N, nrow(paa_old)))#, rfac = r_fac,
             # n_rfac = n_r_ints)
parameters <- list(log_phi = rnorm(1, log(phi), 0.5),
                   p = rep(1, ncol(paa_old))#,
                   # z_rfac = rnorm(n_r_ints, 0, 1),
                   # log_sigma_rfac = rnorm(1, 0.5, 0.5)
                   )

obj <- MakeADFun(data, parameters, DLL="dirichlet_multivariate_rand_anderson")
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
summary(sdr)

obj <- MakeADFun(data, parameters, DLL="dirichlet")
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)






p_hat <- as.numeric(opt$par[-1])
paa_pred <- exp(p_hat) / sum(exp(p_hat))
phi
exp(as.numeric(opt$par[1]))
plot(paa_pred, p);abline(0, 1)
