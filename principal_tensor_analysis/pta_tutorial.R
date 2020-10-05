## Principal Tensor Analysis Tutorial
# Oct. 5, 2020

library(ade4)
library(PTAk)
library(tidyverse)

load(here::here("principal_tensor_analysis", "IBTS_Tensor.Rdata"))


# PCA --------------------------------------------------------------------------

# average over time to reduce dimensionality
IBTS_space <- apply(IBTS_tensor,c(3,1),mean)
# log to remove outliers
IBTS_logspace <- log(IBTS_space+1)

# fit pca
pca_space <- dudi.pca(IBTS_logspace, scale = TRUE, center = TRUE)
inertia.dudi(pca_space)

# visualize pca
s.label(pca_space$li, xax=1, yax=2) #areas
s.label(pca_space$co, xax=1, yax=2, clabel = 0.4) #species

# cluster species based on spatial distribution
#1. Compute the distance between species
dist_species <- dist(pca_space$co, method = "euclidean")

#2. Build a tree with Ward method
den <- hclust(dist_species, method = "ward.D2")

#3. Plot the dendogram
plot(den, hang=-1, ax = T, ann=T, xlab="", sub="", cex=0.6)
nclust <- 5
rect.hclust(den, k = nclust, border="dimgrey")

#4. Create the clusters
clust_space <- as.factor(cutree(den, k = nclust))
s.class(pca_space$co,fac=clust_space, col=rainbow(nclust),xax=1,yax=2)


# PTA --------------------------------------------------------------------------

# log transform the uncollapsed data
IBTS_logtensor <- log(IBTS_tensor+1)

# scale by species
IBTS_logscale <- array(0, dim = dim(IBTS_tensor))

#Loop scanning each species
for (i in 1:nrow(IBTS_tensor)){
  ma <- mean(IBTS_logtensor[i,,])
  sa <- sd(IBTS_logtensor[i,,])
  IBTS_logscale[i,,] <- (IBTS_logtensor[i,,] - ma) / sa
}
dimnames(IBTS_logscale) <- dimnames(IBTS_tensor)

pta <- PTA3(IBTS_logscale, nbPT = 3, nbPT2 = 3, minpct = 0.1)
summary.PTAk(pta, testvar = 0)

# create scree plot
out <- !substr(pta[[3]]$vsnam, 1, 1) == "*"
gct <- (pta[[3]]$pct * pta[[3]]$ssX / pta[[3]]$ssX[1])[out]
barplot(sort(gct, decreasing = TRUE), xlab="PT",
        ylab="Percentage of variance")
# apparent bend after 4th PT

# identify tensor product names based on number selected
tp_keep <- which((pta[[3]]$pct * pta[[3]]$ssX / pta[[3]]$ssX[1]) %in% sort(gct, decreasing = TRUE)[1:4])
pta[[3]]$vsnam[tp_keep]

# plot (mod refers dimensions, nb to the components)
plot(pta, mod=c(2, 3), nb1 = 1, nb2 = 11, xpd=NA, lengthlabels = 4)
plot(pta, mod=1, nb1 = 1, nb2 = 11, lengthlabels = 3)

# these components have the same temporal mode as vs111 (i.e. same temporal
# component) which can be observed by identical distribution of years
plot(pta, mod=c(2,3), nb1 = 1, nb2 = 6, xpd=NA, lengthlabels = 4)
plot(pta, mod=1, nb1 = 1, nb2 = 6, xpd=NA, lengthlabels = 4)
plot(pta, mod=c(2,3), nb1 = 1, nb2 = 7, xpd=NA, lengthlabels = 4)
plot(pta, mod=1, nb1 = 1, nb2 = 7, xpd=NA, lengthlabels = 4)

# Conduct clustering as before
coo <- t(pta[[1]]$v[tp_keep, ])
labkeep <- paste0(pta[[3]]$vsnam[tp_keep], " - ", round((100 * (pta[[3]]$d[tp_keep])^2)/pta[[3]]$ssX[1],1), "%")

#1. Compute the distance between species
dist1 <- dist(coo, method = "euclidean")

#2. Build a tree with Ward linkage
den <- hclust(dist1,method = "ward.D2")

#3. Plot the dendogram
plot(den, hang=-1, ax = T, ann=F, xlab="", sub="",labels = FALSE)

#Choose the number of clusters
nclust <- 6
rect.hclust(den, k=nclust, border=rainbow(nclust)[c(6,5,2,4,3,1)])
