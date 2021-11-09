library(tidyverse)
library(TMB)

# compile(here::here("dirichlet_sdmTMB", "src", "dirichlet.cpp"))
# dyn.load(dynlib(here::here("dirichlet_sdmTMB", "src", "dirichlet")))
compile(
  here::here("dirichlet_sdmTMB", "src", "dirichlet_rand_anderson.cpp")
  )
dyn.load(dynlib(
  here::here("dirichlet_sdmTMB", "src", "dirichlet_rand_anderson"))
  )


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

set.seed(1234)
rdirichlet_linear(p, phi, N)
set.seed(1234)
rdirichlet_linear2(pp * N, phi)


# add REs
n_r_ints <- 5
r_ints <- rnorm(n_r_ints, 0, 2)

n_cat <- 3
log_pp <- rnorm(n_cat, 0, 1) #runif(50)
phi <- 200      # overdispersion parameter
N <- 200        # 200 samples per observation
n_obs <- 40

# original pars
# pp_old <- runif(n_cat)
# p_old <- pp_old / sum(pp_old)


paa <- purrr::map(r_ints, function (xx) {
  pp <- exp(log_pp + xx)
  p <- pp / sum(pp) # 3 categories

  # adj_phi <- exp(log(phi) + xx) # make phi a function of rand int. (IGNORE)

  set.seed(1234)
  naa <- matrix(ncol = n_cat, nrow = n_obs, data = 0) # 40 observations
  for (i in 1:nrow(naa)) naa[i, ] <- rdirmultinom(N, p, phi = phi)

  naa / N

}) %>%
  do.call(rbind, .)


# compare pp and p
purrr::map(r_ints, function (xx) {
   pp <- exp(log_pp + xx)
   p <- pp / sum(pp)
   # adj_phi <- exp(log(phi) + xx)
   list(pp, p)
})



r_fac <- rep(seq(0, n_r_ints - 1, by = 1), each = n_obs)

# check data
colnames(paa2) <- paste("g", seq(1, n_cat, by = 1), sep = "_")
cbind(paa2, r_fac) %>%
  as.data.frame() %>%
  pivot_longer(cols = g_1:g_3, names_to = "group") %>%
  ggplot(.) +
  geom_boxplot(aes(x = group, y = value)) +
  facet_wrap(~r_fac)

data <- list(paa_obs = paa, Neff = rep(N, nrow(paa)), rfac = r_fac,
             n_rfac = n_r_ints)
parameters <- list(log_phi = rnorm(1, log(phi), 0.5),
                   p = rep(1, ncol(paa)),
                   z_rfac = rnorm(n_r_ints, 0, 1),
                   log_sigma_rfac = rnorm(1, 0.5, 0.5)
                   )

obj <- MakeADFun(data, parameters, DLL="dirichlet_rand_anderson")
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)

p_hat <- as.numeric(opt$par[-1])
paa_pred <- exp(p_hat) / sum(exp(p_hat))
phi
exp(as.numeric(opt$par[1]))
plot(paa_pred, p);abline(0, 1)
