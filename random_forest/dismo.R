## Boosted Regression Trees
# Based on dismo vignette: https://geocompr.robinlovelace.net/eco.html
# Example of cloud forest plant community structure combining ordination with
# spatial random forest
# Jan. 24, 2022


library(tidyverse)
library(dismo)
data(Anguilla_train)
head(Anguilla_train)


angaus.tc5.lr01 <- gbm.step(data=Anguilla_train, gbm.x = 3:13, gbm.y = 2,
                            family = "bernoulli", tree.complexity = 5,
                            learning.rate = 0.01, bag.fraction = 0.5)
angaus.tc5.lr005 <- gbm.step(data=Anguilla_train, gbm.x = 3:13, gbm.y = 2,
                            family = "bernoulli", tree.complexity = 5,
                            learning.rate = 0.005, bag.fraction = 0.5)


# evaluate simplified model
angaus.simp <- gbm.simplify(angaus.tc5.lr005, n.drops = 5)

# fit a model with optimal predictors dropped as defined above (note not used
# subsequently)
angaus.tc5.lr005.simp <- gbm.step(Anguilla_train,
                                  gbm.x=angaus.simp$pred.list[[1]], gbm.y=2,
                                  tree.complexity=5, learning.rate=0.005)

gbm.plot(angaus.tc5.lr005, n.plots=11, plot.layout=c(4, 3), write.title = FALSE)


# evaluate interactions
find.int <- gbm.interactions(angaus.tc5.lr005)
find.int$rank.list

gbm.perspec(angaus.tc5.lr005, 7, 1, y.range=c(15,20), z.range=c(0,0.6))


## PREDICTIONS -----------------------------------------------------------------

data(Anguilla_test)
library(gbm)
preds <- predict.gbm(angaus.tc5.lr005, Anguilla_test,
                     n.trees=angaus.tc5.lr005$gbm.call$best.trees,
                     type="response")

# check performance
calc.deviance(obs=Anguilla_test$Angaus_obs, pred=preds, calc.mean=TRUE)

d <- cbind(Anguilla_test$Angaus_obs, preds)
pres <- d[d[,1]==1, 2]
abs <- d[d[,1]==0, 2]
e <- evaluate(p=pres, a=abs)
e

# predict to vector of trees
angaus.5000 <- gbm.fixed(data=Anguilla_train, gbm.x=3:13, gbm.y=2,
                         learning.rate=0.005, tree.complexity=5, n.trees=5000)
## [1] fitting gbm model with a fixed number of 5000 trees for Angaus
## [1] total deviance = 1006.33
## [1] residual deviance = 203.21
tree.list <- seq(100, 5000, by=100)

# generate matrix of predictions (ncol = tree.list, nrow = pred DF)
pred <- predict.gbm(angaus.5000, Anguilla_test, n.trees=tree.list, "response")

# check deviance
angaus.pred.deviance <- rep(0,50)
for (i in 1:50) {
  angaus.pred.deviance[i] <- calc.deviance(Anguilla_test$Angaus_obs,
                                           pred[,i], calc.mean=TRUE)
}

plot(tree.list, angaus.pred.deviance, ylim=c(0.7,1), xlim=c(-100,5000),
     type='l', xlab="number of trees", ylab="predictive deviance",
     cex.lab=1.5)


## SPATIAL PREDICTIONS ---------------------------------------------------------

data(Anguilla_grids)
plot(Anguilla_grids)

Method <- factor('electric', levels = levels(Anguilla_train$Method))
add <- data.frame(Method)
p <- predict(Anguilla_grids, angaus.tc5.lr005, const=add,
             n.trees=angaus.tc5.lr005$gbm.call$best.trees, type="response")
p <- mask(p, raster(Anguilla_grids, 1))
plot(p, main='Angaus - BRT prediction')
