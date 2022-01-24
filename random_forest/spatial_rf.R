## Spatial Random Forest
# Based on blog post: https://geocompr.robinlovelace.net/eco.html
# Example of cloud forest plant community structure combining ordination with
# spatial random forest
# Jan. 24, 2022


library(tidyverse)
library(sf)
library(raster)
library(RQGIS)
library(mlr)
library(dplyr)
library(vegan)


data("study_area", "random_points", "comm", "dem", "ndvi",
     package = "spDataLarge")
data("ep", package = "spDataLarge")

random_points[, names(ep)] = raster::extract(ep, random_points)


# process community data
# presence-absence matrix
pa = decostand(comm, "pa")  # 100 rows (sites), 69 columns (species)
# keep only sites in which at least one species was found
pa = pa[rowSums(pa) != 0, ]  # 84 rows, 69 columns

set.seed(25072018)
nmds = metaMDS(comm = pa, k = 4, try = 500)
nmds$stress

elev = dplyr::filter(random_points, id %in% rownames(pa)) %>%
  dplyr::pull(dem)
# rotating NMDS in accordance with altitude (proxy for humidity)
rotnmds = MDSrotate(nmds, elev)
# extracting the first two axes
sc = scores(rotnmds, choices = 1:2)
# plotting the first axis against altitude
plot(y = sc[, 1], x = elev, xlab = "elevation in m",
     ylab = "First NMDS axis", cex.lab = 0.8, cex.axis = 0.8)


# construct response-predictor matrix
# id- and response variable
rp = data.frame(id = as.numeric(rownames(sc)), sc = sc[, 1])
# join the predictors (dem, ndvi and terrain attributes)
rp = inner_join(random_points, rp, by = "id")


## random forest model
# extract the coordinates into a separate data frame
coords = sf::st_coordinates(rp) %>%
  as.data.frame() %>%
  rename(x = X, y = Y)
# only keep response and predictors which should be used for the modeling
rp = dplyr::select(rp, -id, -spri) %>%
  st_drop_geometry()


## generate mlr building tasks
# create task
task = makeRegrTask(data = rp, target = "sc", coordinates = coords)
# learner
lrn_rf = makeLearner(cl = "regr.ranger", predict.type = "response")


## tune model hyperparameters with spatial cross-validation
# spatial partitioning
perf_level = makeResampleDesc("SpCV", iters = 5)

# specifying random search
ctrl = makeTuneControlRandom(maxit = 50L)

# specifying the search space
ps = makeParamSet(
  makeIntegerParam("mtry", lower = 1, upper = ncol(rp) - 1),
  makeNumericParam("sample.fraction", lower = 0.2, upper = 0.9),
  makeIntegerParam("min.node.size", lower = 1, upper = 10)
)

set.seed(02082018)
tune = tuneParams(learner = lrn_rf,
                  task = task,
                  resampling = perf_level,
                  par.set = ps,
                  control = ctrl,
                  measures = mlr::rmse)


## generate predictions
# learning using the best hyperparameter combination
lrn_rf = makeLearner(cl = "regr.ranger",
                     predict.type = "response",
                     mtry = tune$x$mtry,
                     sample.fraction = tune$x$sample.fraction,
                     min.node.size = tune$x$min.node.size)
# doing the same more elegantly using setHyperPars()
# lrn_rf = setHyperPars(makeLearner("regr.ranger", predict.type = "response"),
#                       par.vals = tune$x)
# train model
model_rf = train(lrn_rf, task)
# to retrieve the ranger output, run:
# mlr::getLearnerModel(model_rf)
# which corresponds to:
# ranger(sc ~ ., data = rp,
#        mtry = tune$x$mtry,
#        sample.fraction = tune$x$sample.fraction,
#        min.node.sie = tune$x$min.node.size)


# prediction dataframe
# convert raster stack into a data frame
new_data = as.data.frame(as.matrix(ep))
# apply the model to the data frame
pred_rf = predict(model_rf, newdata = new_data)
# put the predicted values into a raster
pred = dem
# replace altitudinal values by rf-prediction values
pred[] = pred_rf$data$response
