## caret vignette
# Based on caret vignette:
# https://topepo.github.io/caret/pre-processing.html
# Jan. 24, 2022

library(caret)
library(earth)
library(mlbench)


## DATA CLEANING ---------------------------------------------------------------

# caret requires dummy variables and includes helper function to
# convert
data(etitanic)
head(model.matrix(survived ~ ., data = etitanic))

dummies <- dummyVars(survived ~ ., data = etitanic)
head(predict(dummies, newdata = etitanic))


# near zero variance values
data(mdrr)
data.frame(table(mdrrDescr$nR11))
nzv <- nearZeroVar(mdrrDescr, saveMetrics= TRUE)
nzv[nzv$nzv,][1:10,]
nzv <- nearZeroVar(mdrrDescr)
filteredDescr <- mdrrDescr[, -nzv]
dim(filteredDescr)


# correlated predictors (remove above a cutoff)
descrCor <-  cor(filteredDescr)
sum(abs(descrCor[upper.tri(descrCor)]) > .999)
summary(descrCor[upper.tri(descrCor)])

highlyCorDescr <- findCorrelation(descrCor, cutoff = .75)
filteredDescr <- filteredDescr[,-highlyCorDescr]
descrCor2 <- cor(filteredDescr)
summary(descrCor2[upper.tri(descrCor2)])


# centering and scaling
set.seed(96)
inTrain <- sample(seq(along = mdrrClass), length(mdrrClass)/2)

training <- filteredDescr[inTrain,]
test <- filteredDescr[-inTrain,]
trainMDRR <- mdrrClass[inTrain]
testMDRR <- mdrrClass[-inTrain]

preProcValues <- preProcess(training, method = c("center", "scale"))

trainTransformed <- predict(preProcValues, training)
testTransformed <- predict(preProcValues, test)


# integrated example
library(AppliedPredictiveModeling)
data(schedulingData)
str(schedulingData)

# apply Y-J transform, keep factors as factors
pp_hpc <- preProcess(schedulingData[, -8],
                     method = c("center", "scale", "YeoJohnson"))
pp_hpc
transformed <- predict(pp_hpc, newdata = schedulingData[, -8])
head(transformed)

# check distribution of predictor
mean(schedulingData$NumPending == 0)

# add check for zero-variance predictors
pp_no_nzv <- preProcess(schedulingData[, -8],
                        method = c("center", "scale", "YeoJohnson", "nzv"))
pp_no_nzv
predict(pp_no_nzv, newdata = schedulingData[1:6, -8])


## DATA SPLITTING --------------------------------------------------------------

# create balanced datasets
set.seed(3456)
trainIndex <- createDataPartition(iris$Species, p = .8,
                                  list = FALSE,
                                  times = 1)
irisTrain <- iris[ trainIndex,] %>% glimpse()
irisTest  <- iris[-trainIndex,] %>% glimpse()

# splitting time series using createTimeSlices
createTimeSlices(1:9, 5, 1, fixedWindow = FALSE)
createTimeSlices(1:9, 5, 1, fixedWindow = TRUE)
createTimeSlices(1:9, 5, 3, fixedWindow = TRUE)
createTimeSlices(1:9, 5, 3, fixedWindow = FALSE)

# splitting based on groups
set.seed(3527)
subjects <- sample(1:20, size = 80, replace = TRUE)
table(subjects)
folds <- groupKFold(subjects, k = 15)


## FITTING ---------------------------------------------------------------------

data(Sonar)
str(Sonar[, 1:10])

# stratified splitting of data
set.seed(998)
inTraining <- createDataPartition(Sonar$Class, p = .75, list = FALSE)
training <- Sonar[ inTraining,]
testing  <- Sonar[-inTraining,]

# specify resampling
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)

# fit model
set.seed(825)
gbmFit1 <- train(Class ~ ., data = training,
                 method = "gbm",
                 trControl = fitControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = FALSE)
gbmFit1

# adjust grid space for hyperparameter tuning
gbmGrid <-  expand.grid(interaction.depth = c(1, 5, 9),
                        n.trees = (1:30)*50,
                        shrinkage = 0.1,
                        n.minobsinnode = 20)
nrow(gbmGrid)
set.seed(825)
gbmFit2 <- train(Class ~ ., data = training,
                 method = "gbm",
                 trControl = fitControl,
                 verbose = FALSE,
                 ## Now specify the exact models
                 ## to evaluate:
                 tuneGrid = gbmGrid)
gbmFit2

# visualize resampling profile
trellis.par.set(caretTheme())
plot(gbmFit2)

# rebuild model with slightly different optimization criterion
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using
                           ## the following function
                           summaryFunction = twoClassSummary)

set.seed(825)
gbmFit3 <- train(Class ~ ., data = training,
                 method = "gbm",
                 trControl = fitControl,
                 verbose = FALSE,
                 tuneGrid = gbmGrid,
                 ## Specify which metric to optimize
                 metric = "ROC")
gbmFit3


## PREDICTIONS -----------------------------------------------------------------

# type options for predict standardized to be "class" or "prob"
predict(gbmFit3, newdata = head(testing))
predict(gbmFit3, newdata = head(testing), type = "prob")


## RESAMPLING DISTRIBUTIONS ----------------------------------------------------

trellis.par.set(caretTheme())
densityplot(gbmFit3, pch = "|")

# between model comparison
set.seed(825)
svmFit <- train(Class ~ ., data = training,
                method = "svmRadial",
                trControl = fitControl,
                preProc = c("center", "scale"),
                tuneLength = 8,
                metric = "ROC")
svmFit

set.seed(825)
rdaFit <- train(Class ~ ., data = training,
                method = "rda",
                trControl = fitControl,
                tuneLength = 4,
                metric = "ROC")
rdaFit
