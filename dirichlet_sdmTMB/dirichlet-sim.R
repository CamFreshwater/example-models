## Initial development using simulated data

# rDirichlet with MGLM package
library(MGLM)

n = 10 # number of sampling events
size = 20 # number of observations
#alpha = # vector of parameters = to n categories




obs <- rdirmn(10, 55, alpha = c(0.3, 0.29, 0.35))
obs_mat <- obs / rowSums(obs)
# replace 0 with small values
obs_mat[obs_mat == 0] <- 0.00001
obs_mat <- obs_mat / rowSums(obs_mat)
colnames(obs_mat) <- paste("group", seq(1, 3, by = 1))

obs_dat <- as.data.frame(obs_mat)
obs_dat$Y <- DR_data(obs_dat[,1:3])


mod <- DirichletReg::DirichReg(Y ~ 1, data = obs_dat)

# if (is.vector(alpha) && missing(n)) {
#   stop("When alpha is a vector, must give n.")
# }
# if (is.vector(alpha) && length(alpha) > 1) {
#   alpha <- matrix(alpha, n, length(alpha), byrow = TRUE)
#   if (length(size) == 1) {
#     size <- rep(size, n)
#   }
#   else if (length(size) != n) {
#     stop("The length of size variable doesn't match with sample size")
#   }
# }
# else if (missing(n)) {
#   n <- nrow(alpha)
# }
# else if (length(n) != 1) {
#   stop("n should be a scalar.")
# }
# else if (n != nrow(alpha)) {
#   stop("The sizes of the input don't match. \n\t\talpha must be a vector or a matrix with n rows")
# }
# k <- ncol(alpha)
# if (k <= 1) {
#   stop("The multivariate response data need to have more than one category.")
# }
# if (!is.vector(size)) {
#   stop("Size must be a scalar, or a column vector matches with \n\t\tthe number of rows of alpha")
# }
# else if (length(size) != n)
#   stop("The length of size should match with n")
# G <- sapply(c(alpha), "rgamma", n = 1)
# G <- matrix(G, nrow = n)
# prob <- G/rowSums(G)
# ridx <- rowSums(G) == 0
# if (any(ridx)) {
#   if (sum(ridx) > 1) {
#     prob[ridx, ] <- t(apply((alpha[ridx, ]/rowSums(alpha[ridx,
#     ])), 1, function(x) return(rmultinom(n = 1, size = 1,
#                                          prob = x))))
#   }
#   else {
#     prob[ridx, ] <- rmultinom(n = 1, size = 1, alpha[ridx,
#     ]/sum(alpha[ridx, ]))
#   }
# }
# rdm <- t(sapply(1:n, function(i, size, prob) return(
#   rmultinom(n = 1, size = size[i], prob = prob[i, ])), size, prob)
#   )
# return(rdm)


# rDirichlet with dirmult package
library(dirmult)

K = 5 # number of categories
n_event = 10 # number of sampling events
n_obs = 20 # number of observations per sampling event
# pi = vector of allele probs
theta = 0.03 # overdispersion parameter

if (length(n) == 1) {
  n <- rep(n, J)
}

# if (is.null(pi)) {
  pi <- rnorm(K, mean = 14, sd = 4)
# } else {
#   K <- length(pi)
  # }

pi <- pi/sum(pi)

# P <- rdirichlet(J, pi * (1 - theta)/theta)
# dirmult::rdirichlet <- function (n = 1, alpha)
# {
#   Gam <- matrix(0, n, length(alpha))
#   for (i in 1:length(alpha)) Gam[, i] <- rgamma(n, shape = alpha[i])
#   Gam/rowSums(Gam)
# }

alpha <- pi * (1 - theta) / theta
P <- matrix(0, J, length(alpha))
for (i in 1:length(alpha)) P[, i] <- rgamma(J, shape = alpha[i])
P <- P/rowSums(P)

X <- matrix(0, J, K)

for (i in 1:J) X[i, ] <- rmultinom(1, n[i], P[i, ])

list(theta = theta, pi = pi, data = X)





set.seed(123)
n <- 200
d <- 4
alpha <- rep(1, d)
m <- 50
Y <- rdirmn(n, m, alpha)
