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


gbm.plot(angaus.tc5.lr005, n.plots=11, plot.layout=c(4, 3), write.title = FALSE)
