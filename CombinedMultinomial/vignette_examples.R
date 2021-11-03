## Based on combined multinomial repo and vignette:
# file:///C:/github/CombinedMultinomial/vignettes/Vignette_CombinedMultinomial.html

library(CombinedMultinomial)
library(HMP)

## can't load package so source functions
files.sources = list.files(here::here("CombinedMultinomial", "package_functions"))
for (i in seq_along(files.sources)) {
  source(
    here::here("CombinedMultinomial", "package_functions", files.sources[[i]])
    )
}


# make design matrix
set.seed(12345)
Q <- 3; lvl.cov <- list(2);
FunDes <- DesMat(Q,lvl.cov)
FunDes

# simulate data
N <- 100     # number of observations
Q <- 3       # number of categories
S <- 2000    # counts per observation
s.u <- 0.88  # sd of RE
theta <- 0.2 # overdispersion parameter
n.cov <- 1   # number of fixed covariates
Des <- FunDes; lvl.cov <- list(2); method = "DMM";
beta <- c(-1,0.1,0.5,0.8,-2) # presumably coefficients for fixed effects?

SimDat <- DataGen(N = N, Q = Q, S = S, beta = beta, theta = theta, s.u = s.u,
                  n.cov = n.cov, lvl.cov = lvl.cov, Des = Des, method = method)

head(SimDat)
