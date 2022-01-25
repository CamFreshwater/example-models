## caret vs. tidymodels
# Based on dismo vignette:
# https://towardsdatascience.com/caret-vs-tidymodels-how-to-use-both-packages-together-ee3f85b381c
# Based on bike sharing
# Jan. 24, 2022

library(tidymodels)
library(caret)
library(lubridate)
library(tidyverse)
library(moments)
library(corrr)
library(randomForest)


bike <- read.csv(here::here("random_forest", "bike_share_data", "hour.csv")) %>%
  mutate(instant = NULL, yr = yr + 2011) %>%
  rename(
    date = dteday,
    year = yr,
    month = mnth,
    hour = hr,
    weather = weathersit,
    humidity = hum,
    total = cnt
  )


## Observe raw data
bike %>%
  pivot_longer(
    cols = c(casual, registered, total),
    names_to = "usertype",
    values_to = "count"
  ) %>%
  ggplot(aes(count, colour = usertype)) +
  geom_density() +
  labs(
    title = "Distribution of the number of rental bikes",
    x = "Number per hour", y = "Density"
  ) +
  scale_colour_discrete(
    name = "User type",
    breaks = c("casual", "registered", "total"),
    labels = c("Non-registered", "Registered", "Total")
  )

bike_cor <- bike %>%
  select(where(is.numeric)) %>%
  correlate() %>%
  rearrange(absolute = FALSE) %>%
  shave()

rplot(bike_cor, print_cor = TRUE)


## Transform variables
bike_all <- bike %>%
  select(-casual, -registered)

# Transform with cubic root
bike_all$total <- bike_all$total^(1 / 3)

bike_all$season <- factor(
  bike_all$season,
  levels = c(1, 2, 3, 4),
  labels = c("spring", "summer", "autumn", "winter")
)
bike_all$holiday <- factor(
  bike_all$holiday,
  levels = c(0, 1), labels = c(FALSE, TRUE)
)
bike_all$workingday <- factor(
  bike_all$workingday,
  levels = c(0, 1), labels = c(FALSE, TRUE)
)
bike_all$weather <- factor(
  bike_all$weather,
  levels = c(1, 2, 3, 4),
  labels = c("clear", "cloudy", "rainy", "heavy rain"),
  ordered = TRUE
)


## Split training and testing
set.seed(25)
split <- initial_split(bike_all, prop = 0.8)
train_data <- training(split)
train_data %>% dim()

test_data <- testing(split)
test_data %>% dim()

train_cv <- vfold_cv(train_data, v = 10)


## Additional data cleaning
prep_recipe <- recipe(total ~ ., data = train_data) %>%
  step_rm(year, month, weekday) %>%
  step_date(date) %>%
  step_corr(all_numeric(), threshold = 0.8) %>%
  step_dummy(all_nominal())


## Custom fnctions
# Generate prediction tables
predict_table <- function(model, data, tidy_flag) {
  if (tidy_flag == TRUE) {
    result <- model %>%
      predict(data) %>%
      rename(pred = .pred) %>%
      mutate(
        actual = data$total,
        pred_real = pred^3,
        actual_real = actual^3
      )
  } else {
    result <- model %>%
      predict(data) %>%
      as_tibble_col(column_name = "pred") %>%
      mutate(
        actual = data$total,
        pred_real = pred^3,
        actual_real = actual^3
      )
  }
  result
}

# Extract RMSE for models
pull_rmse <- function(result_table) {
  rmse_result <- rmse(result_table, pred, actual) %>%
    pull(.estimate)
  rmse_result_real <- rmse(result_table, pred_real, actual_real) %>%
    pull(.estimate)
  result <- c(rmse = rmse_result, real_rmse = rmse_result_real)
}


base_train_pred <-
  tibble(
    actual = train_data$total,
    actual_real = train_data$total^3
  ) %>%
  mutate(pred = mean(actual), pred_real = mean(actual_real))
base_test_pred <-
  tibble(
    actual = test_data$total,
    actual_real = test_data$total^3
  ) %>%
  mutate(pred = mean(actual), pred_real = mean(actual_real))
base_train_rmse <- pull_rmse(base_train_pred)
print(base_train_rmse)
base_test_rmse <- pull_rmse(base_test_pred)
print(base_test_rmse)


# Cost complexity for decision tree parameter
tree_cp <- seq(0.01, 0.1, 0.01)

set.seed(26)
tree_tidy_time1 <- Sys.time()

# Specify model
tree_engine <-
  decision_tree(mode = "regression", cost_complexity = tune()) %>%
  set_engine("rpart")

# Set workflow (Preprocess & model)
tree_workflow <-
  workflow() %>%
  add_recipe(prep_recipe) %>%
  add_model(tree_engine)

# Tune parameters with cross-validation
tree_tune <- tune_grid(
  tree_workflow,
  resamples = train_cv,
  grid = data.frame(cost_complexity = tree_cp),
  metrics = metric_set(rmse)
)

# Fit again with the best parameter
tree_best <-
  finalize_workflow(tree_workflow, select_best(tree_tune)) %>%
  fit(train_data)
