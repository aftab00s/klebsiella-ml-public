# Run from the project root: Rscript -e 'testthat::test_dir("tests/testthat")'

source(here::here("R", "config.R"))
source(here::here("R", "preprocessing.R"))

# A small made-up training set: numbers with one gap, and categories with one gap.
train <- data.frame(
  outcome = factor(rep(c("No", "Yes"), 5), levels = c("No", "Yes")),
  number = c(1, 2, 3, 4, NA, 6, 7, 8, 9, 10),
  category = factor(c("a", "b", NA, "b", "a", "b", "a", "b", "a", "b")),
  constant = rep(5, 10)
)

test_set <- data.frame(
  outcome = factor(c("No", "Yes"), levels = c("No", "Yes")),
  number = c(NA, 1000),
  category = factor(c(NA, "a")),
  constant = c(5, 5)
)

result <- preprocess_train_test(train, "outcome", c("No", "Yes"), test_set)

test_that("missing test values are filled using the training median only", {
  training_median <- median(train$number, na.rm = TRUE)
  expect_equal(result$test$number[1], training_median)

  # A huge value in the test set must not change anything learned from training.
  other_test_set <- test_set
  other_test_set$number[2] <- 999999
  other_result <- preprocess_train_test(train, "outcome", c("No", "Yes"), other_test_set)
  expect_equal(other_result$train, result$train)
  expect_equal(other_result$test$number[1], training_median)
})

test_that("no missing values remain after preprocessing", {
  expect_false(anyNA(result$train))
  expect_false(anyNA(result$test))
})

test_that("output has the right dimensions and matching columns", {
  expect_equal(nrow(result$train), nrow(train))
  expect_equal(nrow(result$test), nrow(test_set))
  expect_equal(names(result$train), names(result$test))
  expect_false("constant" %in% names(result$train))
})
