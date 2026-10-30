# Model fitting: logistic regression, multinomial regression and random forest.

# Fit logistic regression inside a safe wrapper so one failed fold does not stop
# the whole run. Warnings are suppressed because small folds often trigger
# separation warnings; diagnose_glm_model() records those problems in the notes.
fit_glm_safely <- function(formula, data) {
  suppressWarnings(
    tryCatch(
      {
        model_data <- as.data.frame(data)
        fit <- stats::glm(
          formula = formula,
          data = model_data,
          family = stats::binomial(),
          control = stats::glm.control(maxit = 100)
        )
        # stepAIC refits the model later by re-reading fit$call. The data only
        # exists inside this function, so it is stored in the call to keep it
        # available after this function returns.
        fit$call$data <- model_data
        fit
      },
      error = function(e) NULL
    )
  )
}

# Stepwise selection in both directions, falling back to the full model if it
# fails.
fit_stepwise_model <- function(null_model, full_model) {
  if (is.null(full_model)) {
    return(null_model)
  }

  suppressWarnings(
    tryCatch(
      MASS::stepAIC(
        object = null_model,
        scope = list(lower = stats::formula(null_model), upper = stats::formula(full_model)),
        direction = "both",
        trace = FALSE
      ),
      error = function(e) full_model
    )
  )
}

# Generate probabilities from the fitted GLM and fall back to the observed
# event rate if a fold-specific prediction call fails.
predict_with_glm <- function(model, newdata, fallback_rate) {
  if (is.null(model)) {
    return(rep(fallback_rate, nrow(newdata)))
  }

  predictions <- suppressWarnings(
    tryCatch(
      stats::predict(model, newdata = newdata, type = "response"),
      error = function(e) rep(fallback_rate, nrow(newdata))
    )
  )

  clip_probabilities(predictions)
}

# Fit multinomial regression with a stored model frame so later coefficient
# extraction remains stable on this R installation.
fit_multinom_safely <- function(model_formula, model_data) {
  model_formula <- stats::as.formula(model_formula)
  model_data <- as.data.frame(model_data)

  suppressWarnings(
    tryCatch(
      {
        # do.call builds the multinom() call from this list. As in fit_glm_safely,
        # the formula and data are stored in the call for later use.
        fit <- do.call(
          nnet::multinom,
          list(
            formula = model_formula,
            data = model_data,
            trace = FALSE,
            MaxNWts = 50000,
            model = TRUE
          )
        )
        fit$call$formula <- model_formula
        fit$call$data <- model_data
        fit
      },
      error = function(e) NULL
    )
  )
}

# If a class is missing from a fold's training data, the model returns no column
# for it. This puts every expected class back so all folds have the same columns.
align_probability_columns <- function(prob_matrix, outcome_levels) {
  if (is.null(prob_matrix)) {
    return(matrix(0, nrow = 0, ncol = length(outcome_levels), dimnames = list(NULL, outcome_levels)))
  }

  prob_matrix <- as.matrix(prob_matrix)

  if (is.null(colnames(prob_matrix))) {
    if (length(outcome_levels) == 2 && ncol(prob_matrix) == 1) {
      colnames(prob_matrix) <- outcome_levels[2]
    } else if (ncol(prob_matrix) == length(outcome_levels)) {
      colnames(prob_matrix) <- outcome_levels
    }
  }

  aligned <- matrix(
    0,
    nrow = nrow(prob_matrix),
    ncol = length(outcome_levels),
    dimnames = list(rownames(prob_matrix), outcome_levels)
  )

  shared_columns <- intersect(colnames(prob_matrix), outcome_levels)
  if (length(shared_columns) > 0) {
    aligned[, shared_columns] <- prob_matrix[, shared_columns, drop = FALSE]
  }

  row_sums <- rowSums(aligned)
  empty_rows <- row_sums == 0
  if (any(empty_rows)) {
    aligned[empty_rows, ] <- 1 / length(outcome_levels)
    row_sums <- rowSums(aligned)
  }

  aligned / row_sums
}

# Convert multinomial probabilities into a clean prediction object with
# explicit class labels and fallback probabilities if fitting fails.
predict_with_multinom <- function(model, newdata, outcome_levels, fallback_probs) {
  if (is.null(model)) {
    probs <- matrix(
      rep(fallback_probs, each = nrow(newdata)),
      nrow = nrow(newdata),
      byrow = FALSE,
      dimnames = list(NULL, outcome_levels)
    )
  } else {
    raw_probs <- suppressWarnings(
      tryCatch(
        stats::predict(model, newdata = newdata, type = "probs"),
        error = function(e) NULL
      )
    )

    if (is.null(raw_probs)) {
      probs <- matrix(
        rep(fallback_probs, each = nrow(newdata)),
        nrow = nrow(newdata),
        byrow = FALSE,
        dimnames = list(NULL, outcome_levels)
      )
    } else {
      if (is.vector(raw_probs)) {
        raw_probs <- matrix(raw_probs, nrow = nrow(newdata))
      }
      probs <- align_probability_columns(raw_probs, outcome_levels)
    }
  }

  predicted_class <- factor(
    outcome_levels[max.col(probs, ties.method = "first")],
    levels = outcome_levels
  )

  list(probabilities = probs, predicted_class = predicted_class)
}

# Extract interpretable binary-model coefficients including odds ratios and
# confidence intervals for dissertation reporting.
extract_glm_coefficients <- function(model, positive_class) {
  if (is.null(model) || length(stats::coef(model)) <= 1) {
    return(tibble::tibble())
  }

  summary(model)$coefficients |>
    as.data.frame() |>
    tibble::rownames_to_column("term") |>
    dplyr::rename(
      estimate = Estimate,
      std_error = `Std. Error`,
      z_value = `z value`,
      p_value = `Pr(>|z|)`
    ) |>
    dplyr::filter(term != "(Intercept)") |>
    dplyr::mutate(
      odds_ratio = exp(estimate),
      conf_low = exp(estimate - 1.96 * std_error),
      conf_high = exp(estimate + 1.96 * std_error),
      interpretation = ifelse(
        odds_ratio > 1,
        paste("Higher odds of", positive_class),
        paste("Lower odds of", positive_class)
      )
    ) |>
    dplyr::arrange(p_value, dplyr::desc(abs(z_value)))
}

# Extract class-specific multinomial coefficients so phenotype associations can
# still be described even though the endpoint is multiclass.
extract_multinom_coefficients <- function(model) {
  if (is.null(model)) {
    return(tibble::tibble())
  }

  model_summary <- summary(model)
  coefficients <- model_summary$coefficients
  std_errors <- model_summary$standard.errors

  if (is.null(dim(coefficients))) {
    coefficients <- matrix(coefficients, nrow = 1, dimnames = list(model$lev[2], names(coefficients)))
    std_errors <- matrix(std_errors, nrow = 1, dimnames = list(model$lev[2], names(std_errors)))
  }

  coef_df <- as.data.frame(coefficients) |>
    tibble::rownames_to_column("class")
  se_df <- as.data.frame(std_errors) |>
    tibble::rownames_to_column("class")

  coef_long <- tidyr::pivot_longer(coef_df, cols = -class, names_to = "term", values_to = "estimate")
  se_long <- tidyr::pivot_longer(se_df, cols = -class, names_to = "term", values_to = "std_error")

  dplyr::left_join(coef_long, se_long, by = c("class", "term")) |>
    dplyr::filter(term != "(Intercept)") |>
    dplyr::mutate(
      z_value = estimate / std_error,
      p_value = 2 * stats::pnorm(abs(z_value), lower.tail = FALSE),
      relative_risk_ratio = exp(estimate),
      conf_low = exp(estimate - 1.96 * std_error),
      conf_high = exp(estimate + 1.96 * std_error)
    ) |>
    dplyr::arrange(class, p_value, dplyr::desc(abs(z_value)))
}

# Build the grid of random forest settings to try: mtry, ntree and nodesize.
make_rf_grid <- function(n_predictors, outcome_type = "binary") {
  candidate_values <- c(
    1,
    2,
    round(sqrt(n_predictors)),
    round(n_predictors / 10),
    round(n_predictors / 5)
  )

  mtry_values <- sort(unique(pmax(1, pmin(n_predictors, candidate_values))))
  ntree_values <- c(500L, 1000L)
  nodesize_values <- if (outcome_type == "binary") c(1L, 3L, 5L) else c(1L, 3L)

  expand.grid(
    mtry = mtry_values,
    ntree = ntree_values,
    nodesize = nodesize_values,
    KEEP.OUT.ATTRS = FALSE
  ) |>
    tibble::as_tibble() |>
    dplyr::arrange(mtry, ntree, nodesize)
}

# Fit random forest with inverse-frequency class weights so minority classes are
# not ignored during training.
fit_rf_safely <- function(train_df, outcome_name, mtry, ntree, nodesize = 1L, proximity = FALSE) {
  predictors <- setdiff(names(train_df), outcome_name)
  y <- droplevels(train_df[[outcome_name]])

  if (length(predictors) == 0 || nlevels(y) < 2) {
    return(NULL)
  }

  class_counts <- table(y)
  class_weights <- sum(class_counts) / class_counts
  class_weights <- class_weights / sum(class_weights)

  suppressWarnings(
    tryCatch(
      randomForest::randomForest(
        x = train_df[, predictors, drop = FALSE],
        y = y,
        mtry = min(mtry, length(predictors)),
        ntree = ntree,
        nodesize = nodesize,
        importance = TRUE,
        proximity = proximity,
        keep.forest = TRUE,
        classwt = class_weights
      ),
      error = function(e) NULL
    )
  )
}

# Predict from random forest and keep a probability matrix aligned to the full
# expected outcome space.
predict_with_rf <- function(model, newdata, outcome_levels, fallback_probs) {
  if (is.null(model)) {
    probs <- matrix(
      rep(fallback_probs, each = nrow(newdata)),
      nrow = nrow(newdata),
      dimnames = list(NULL, outcome_levels)
    )
  } else {
    raw_probs <- suppressWarnings(
      tryCatch(
        predict(model, newdata = newdata, type = "prob"),
        error = function(e) NULL
      )
    )

    if (is.null(raw_probs)) {
      probs <- matrix(
        rep(fallback_probs, each = nrow(newdata)),
        nrow = nrow(newdata),
        dimnames = list(NULL, outcome_levels)
      )
    } else {
      if (is.vector(raw_probs)) {
        raw_probs <- matrix(raw_probs, nrow = nrow(newdata))
      }
      probs <- align_probability_columns(raw_probs, outcome_levels)
    }
  }

  predicted_class <- factor(
    outcome_levels[max.col(probs, ties.method = "first")],
    levels = outcome_levels
  )

  list(probabilities = probs, predicted_class = predicted_class)
}

# Turn fold-level random-forest predictions into a tuning summary table so the
# best mtry, ntree, and nodesize setting can be selected transparently.
summarise_rf_grid <- function(prediction_df, spec) {
  # Make one group per combination of settings, so each is scored separately.
  split_key <- interaction(
    prediction_df$mtry,
    prediction_df$ntree,
    prediction_df$nodesize,
    drop = TRUE
  )
  split_predictions <- split(prediction_df, split_key)

  grid_summary <- dplyr::bind_rows(lapply(split_predictions, function(df) {
    metric_object <- if (spec$type == "binary") {
      compute_classification_metrics(
        actual = df$actual,
        predicted = df$predicted_class,
        probability = df$prob_positive,
        model_name = "Random forest",
        outcome_type = spec$type,
        positive_class = spec$positive_class
      )
    } else {
      compute_classification_metrics(
        actual = df$actual,
        predicted = df$predicted_class,
        probability = NULL,
        model_name = "Random forest",
        outcome_type = spec$type
      )
    }

    metric_object$summary |>
      dplyr::mutate(
        mtry = unique(df$mtry),
        ntree = unique(df$ntree),
        nodesize = unique(df$nodesize),
        outcome = spec$display_name,
        model_type = "Random forest",
        n_used = length(unique(df$sample_id)),
        n_predictions = nrow(df)
      )
  }))

  if (spec$type == "binary") {
    grid_summary |>
      dplyr::mutate(notes = "") |>
      dplyr::arrange(
        dplyr::desc(balanced_accuracy),
        dplyr::desc(auc),
        dplyr::desc(accuracy),
        dplyr::desc(kappa),
        mtry,
        ntree,
        nodesize
      )
  } else {
    grid_summary |>
      dplyr::mutate(notes = "") |>
      dplyr::arrange(
        dplyr::desc(balanced_accuracy),
        dplyr::desc(accuracy),
        dplyr::desc(kappa),
        mtry,
        ntree,
        nodesize
      )
  }
}

# Extract and rank random-forest importance in one consistent table regardless
# of whether the endpoint is binary or multiclass.
extract_rf_importance <- function(model) {
  if (is.null(model)) {
    return(tibble::tibble())
  }

  importance_df <- randomForest::importance(model) |>
    as.data.frame() |>
    tibble::rownames_to_column("predictor")

  if ("MeanDecreaseGini" %in% names(importance_df)) {
    importance_df$importance_score <- importance_df$MeanDecreaseGini
  } else {
    numeric_cols <- vapply(importance_df, is.numeric, logical(1))
    importance_df$importance_score <- rowMeans(importance_df[, numeric_cols, drop = FALSE], na.rm = TRUE)
  }

  importance_df |>
    dplyr::arrange(dplyr::desc(importance_score))
}

# Collapse diagnostic notes into one readable cell for dissertation appendices.
collapse_notes <- function(notes) {
  notes <- unique(notes[!is.na(notes) & nzchar(notes)])
  if (length(notes) == 0) "" else paste(notes, collapse = " | ")
}

# Flag common GLM problems without stopping the pipeline. This is intentionally
# conservative: warnings are exported for interpretation rather than used to
# discard models automatically.
diagnose_glm_model <- function(model) {
  if (is.null(model)) {
    return("Model fit failed; metrics use event-rate fallback for affected folds.")
  }

  notes <- character()
  if (!isTRUE(model$converged)) {
    notes <- c(notes, "Model did not report convergence.")
  }

  coefficients <- suppressWarnings(stats::coef(model))
  if (any(is.na(coefficients))) {
    notes <- c(notes, "One or more coefficients are aliased or singular.")
  }
  if (any(abs(coefficients[is.finite(coefficients)]) > 20)) {
    notes <- c(notes, "Very large coefficient detected; possible separation or sparse levels.")
  }

  collapse_notes(notes)
}

# Run stepwise GLM selection safely in a named direction. A failed selection
# falls back to the starting model and says so in the notes.
fit_stepwise_glm_safely <- function(start_model, lower_formula, upper_formula, direction, fallback_model = NULL) {
  if (is.null(start_model)) {
    return(list(model = fallback_model, notes = "Stepwise selection skipped because the starting model failed."))
  }

  model <- suppressWarnings(
    tryCatch(
      MASS::stepAIC(
        object = start_model,
        scope = list(lower = lower_formula, upper = upper_formula),
        direction = direction,
        trace = FALSE
      ),
      error = function(e) NULL
    )
  )

  if (is.null(model)) {
    return(list(
      model = fallback_model,
      notes = collapse_notes(c("Stepwise selection failed.", diagnose_glm_model(fallback_model)))
    ))
  }

  list(model = model, notes = diagnose_glm_model(model))
}

# Fit the five binary GLM variants on a supplied training
# set. The full model is screened first so fold-specific dimensionality remains
# plausible on a small biomedical dataset.
fit_binary_glm_variant_set <- function(train_df, spec) {
  outcome_name <- spec$cleaned_outcome
  selection <- select_screened_predictors(train_df, outcome_name, spec$type, config$max_binary_df)

  null_formula <- stats::as.formula(paste(outcome_name, "~ 1"))
  full_formula <- if (length(selection$selected_predictors) > 0) {
    stats::reformulate(selection$selected_predictors, response = outcome_name)
  } else {
    null_formula
  }

  null_model <- fit_glm_safely(null_formula, train_df)
  full_model <- fit_glm_safely(full_formula, train_df)

  forward_fit <- fit_stepwise_glm_safely(
    start_model = null_model,
    lower_formula = null_formula,
    upper_formula = full_formula,
    direction = "forward",
    fallback_model = null_model
  )
  backward_fit <- fit_stepwise_glm_safely(
    start_model = full_model,
    lower_formula = null_formula,
    upper_formula = full_formula,
    direction = "backward",
    fallback_model = full_model %||% null_model
  )
  both_fit <- fit_stepwise_glm_safely(
    start_model = null_model,
    lower_formula = null_formula,
    upper_formula = full_formula,
    direction = "both",
    fallback_model = full_model %||% null_model
  )

  variants <- list(
    "Null model" = list(model = null_model, notes = diagnose_glm_model(null_model)),
    "Forward selection model" = forward_fit,
    "Backward selection model" = backward_fit,
    "Stepwise both-direction model" = both_fit,
    "Full screened candidate model" = list(model = full_model, notes = diagnose_glm_model(full_model))
  )

  variants <- lapply(variants, function(variant) {
    variant$formula_text <- if (is.null(variant$model)) {
      NA_character_
    } else {
      paste(deparse(stats::formula(variant$model)), collapse = " ")
    }
    variant
  })

  list(variants = variants, selection = selection)
}

# Compare the five GLM variants through the same fold plan used by the primary
# binary workflow. Failed fold fits are documented and use the training event
# rate as a transparent fallback so the run completes.
run_binary_glm_model_comparison <- function(dataset_df, spec, fold_plan) {
  outcome_name <- spec$cleaned_outcome
  outcome_levels <- levels(dataset_df[[outcome_name]])
  variant_names <- c(
    "Null model",
    "Forward selection model",
    "Backward selection model",
    "Stepwise both-direction model",
    "Full screened candidate model"
  )

  prediction_records <- list()
  note_records <- list()

  for (fold_id in seq_along(fold_plan$folds)) {
    test_index <- fold_plan$folds[[fold_id]]
    train_index <- setdiff(seq_len(nrow(dataset_df)), test_index)

    processed <- preprocess_train_test(
      dataset_df[train_index, , drop = FALSE],
      outcome_name,
      outcome_levels,
      dataset_df[test_index, , drop = FALSE]
    )

    train_df <- processed$train
    test_df <- processed$test
    variant_set <- fit_binary_glm_variant_set(train_df, spec)
    fallback_rate <- mean(train_df[[outcome_name]] == spec$positive_class)

    for (model_name in variant_names) {
      variant <- variant_set$variants[[model_name]]
      predicted_prob <- predict_with_glm(variant$model, test_df, fallback_rate)
      predicted_class <- factor(
        ifelse(predicted_prob >= 0.5, spec$positive_class, spec$negative_class),
        levels = outcome_levels
      )

      prediction_records[[paste(fold_id, model_name, sep = "_")]] <- tibble::tibble(
        sample_id = rownames(test_df),
        fold = fold_id,
        outcome = spec$display_name,
        model_name = model_name,
        actual = test_df[[outcome_name]],
        predicted_class = predicted_class,
        predicted_probability = predicted_prob
      )

      note_records[[paste(fold_id, model_name, sep = "_")]] <- tibble::tibble(
        model_name = model_name,
        notes = variant$notes
      )
    }
  }

  all_predictions <- dplyr::bind_rows(prediction_records)
  all_notes <- dplyr::bind_rows(note_records)

  full_processed <- preprocess_train_test(dataset_df, outcome_name, outcome_levels)$train
  final_variant_set <- fit_binary_glm_variant_set(full_processed, spec)

  dplyr::bind_rows(lapply(variant_names, function(model_name) {
    model_predictions <- all_predictions[all_predictions$model_name == model_name, ]

    metric_object <- tryCatch(
      compute_classification_metrics(
        actual = model_predictions$actual,
        predicted = model_predictions$predicted_class,
        probability = model_predictions$predicted_probability,
        model_name = model_name,
        outcome_type = spec$type,
        positive_class = spec$positive_class
      ),
      error = function(e) NULL
    )

    if (is.null(metric_object)) {
      metrics <- tibble::tibble(
        accuracy = NA_real_,
        sensitivity = NA_real_,
        specificity = NA_real_,
        balanced_accuracy = NA_real_,
        kappa = NA_real_,
        auc = NA_real_,
        p_value_vs_nir = NA_real_
      )
      metric_note <- "Metric calculation failed."
    } else {
      metrics <- metric_object$summary |>
        dplyr::select(accuracy, sensitivity, specificity, balanced_accuracy, kappa, auc, p_value_vs_nir)
      metric_note <- ""
    }

    final_variant <- final_variant_set$variants[[model_name]]
    note_text <- collapse_notes(c(
      all_notes$notes[all_notes$model_name == model_name],
      final_variant$notes,
      metric_note
    ))

    tibble::tibble(
      outcome = spec$display_name,
      model_name = model_name,
      # The metrics describe the whole procedure across folds, where selected
      # variables can differ. This formula comes from a separate full-data fit.
      formula = final_variant$formula_text,
      n_used = nrow(dataset_df)
    ) |>
      dplyr::bind_cols(metrics) |>
      dplyr::mutate(notes = note_text)
  }))
}

# Run the full binary GLM workflow: fold-wise preprocessing, predictor
# screening, null/full/stepwise modelling, and final interpretability outputs.
run_binary_glm_pipeline <- function(dataset_df, spec, fold_plan) {
  outcome_name <- spec$cleaned_outcome
  outcome_levels <- levels(dataset_df[[outcome_name]])

  cv_predictions <- list()
  screening_tables <- list()

  for (fold_id in seq_along(fold_plan$folds)) {
    test_index <- fold_plan$folds[[fold_id]]
    train_index <- setdiff(seq_len(nrow(dataset_df)), test_index)

    train_raw <- dataset_df[train_index, , drop = FALSE]
    test_raw <- dataset_df[test_index, , drop = FALSE]

    # Rebuild the preprocessing inside every fold so no information from the
    # test fold leaks into imputation or level handling.
    processed <- preprocess_train_test(train_raw, outcome_name, outcome_levels, test_raw)
    train_df <- processed$train
    test_df <- processed$test

    selection <- select_screened_predictors(train_df, outcome_name, spec$type, config$max_binary_df)
    selected_predictors <- selection$selected_predictors

    # Fit null, full and stepwise models, keeping the full model to a screened
    # predictor set that small folds can support.
    null_formula <- stats::as.formula(paste(outcome_name, "~ 1"))
    full_formula <- if (length(selected_predictors) > 0) {
      stats::reformulate(selected_predictors, response = outcome_name)
    } else {
      null_formula
    }

    null_model <- fit_glm_safely(null_formula, train_df)
    full_model <- fit_glm_safely(full_formula, train_df)
    step_model <- fit_stepwise_model(null_model, full_model)

    fallback_rate <- mean(train_df[[outcome_name]] == spec$positive_class)
    predicted_prob <- predict_with_glm(step_model, test_df, fallback_rate)
    predicted_class <- factor(
      ifelse(predicted_prob >= 0.5, spec$positive_class, spec$negative_class),
      levels = outcome_levels
    )

    cv_predictions[[fold_id]] <- tibble::tibble(
      sample_id = rownames(test_df),
      fold = fold_id,
      actual = test_df[[outcome_name]],
      glm_probability = predicted_prob,
      glm_prediction = predicted_class
    )

    screening_tables[[fold_id]] <- selection$score_table |>
      dplyr::mutate(fold = fold_id)
  }

  cv_predictions <- dplyr::bind_rows(cv_predictions) |>
    dplyr::arrange(sample_id)

  cv_metrics <- compute_classification_metrics(
    actual = cv_predictions$actual,
    predicted = cv_predictions$glm_prediction,
    probability = cv_predictions$glm_probability,
    model_name = "GLM (stepwise)",
    outcome_type = spec$type,
    positive_class = spec$positive_class
  )

  # Refit the chosen binary model on the full cleaned dataset so final
  # coefficients and odds ratios use all available information.
  full_processed <- preprocess_train_test(dataset_df, outcome_name, outcome_levels)$train
  full_selection <- select_screened_predictors(full_processed, outcome_name, spec$type, config$max_binary_df)

  null_formula <- stats::as.formula(paste(outcome_name, "~ 1"))
  full_formula <- if (length(full_selection$selected_predictors) > 0) {
    stats::reformulate(full_selection$selected_predictors, response = outcome_name)
  } else {
    null_formula
  }

  null_model <- fit_glm_safely(null_formula, full_processed)
  full_model <- fit_glm_safely(full_formula, full_processed)
  step_model <- fit_stepwise_model(null_model, full_model)

  coefficients <- extract_glm_coefficients(step_model, spec$positive_class)

  model_fit <- tibble::tibble(
    model = c("Null model", "Full model", "Stepwise model"),
    formula = c(
      paste(deparse(stats::formula(null_model)), collapse = " "),
      paste(deparse(stats::formula(full_model)), collapse = " "),
      paste(deparse(stats::formula(step_model)), collapse = " ")
    ),
    aic = c(stats::AIC(null_model), stats::AIC(full_model), stats::AIC(step_model))
  )

  list(
    cv_predictions = cv_predictions,
    cv_metrics = cv_metrics,
    screening = dplyr::bind_rows(screening_tables),
    coefficients = coefficients,
    model_fit = model_fit,
    selected_predictors = full_selection$selected_predictors,
    final_model = step_model
  )
}

# Run the binary random-forest workflow over the same folds and a small tuning
# grid so the binary outcomes are evaluated in parallel with the GLM models.
run_binary_rf_pipeline <- function(dataset_df, spec, fold_plan) {
  outcome_name <- spec$cleaned_outcome
  outcome_levels <- levels(dataset_df[[outcome_name]])
  full_processed <- preprocess_train_test(dataset_df, outcome_name, outcome_levels)$train
  predictor_count <- ncol(full_processed) - 1
  rf_grid <- make_rf_grid(predictor_count, spec$type)

  all_predictions <- list()

  # Evaluate each RF grid setting by rebuilding the same fold-wise preprocessing
  # used in the GLM workflow; this keeps tuning and imputation leakage-safe.
  for (grid_id in seq_len(nrow(rf_grid))) {
    mtry_value <- rf_grid$mtry[[grid_id]]
    ntree_value <- rf_grid$ntree[[grid_id]]
    nodesize_value <- rf_grid$nodesize[[grid_id]]
    fold_predictions <- list()

    for (fold_id in seq_along(fold_plan$folds)) {
      test_index <- fold_plan$folds[[fold_id]]
      train_index <- setdiff(seq_len(nrow(dataset_df)), test_index)

      train_raw <- dataset_df[train_index, , drop = FALSE]
      test_raw <- dataset_df[test_index, , drop = FALSE]

      processed <- preprocess_train_test(train_raw, outcome_name, outcome_levels, test_raw)
      train_df <- processed$train
      test_df <- processed$test

      class_prior <- prop.table(table(train_df[[outcome_name]]))
      fallback_probs <- setNames(rep(0, length(outcome_levels)), outcome_levels)
      fallback_probs[names(class_prior)] <- as.numeric(class_prior)

      rf_model <- fit_rf_safely(
        train_df,
        outcome_name,
        mtry = mtry_value,
        ntree = ntree_value,
        nodesize = nodesize_value,
        proximity = FALSE
      )
      prediction <- predict_with_rf(rf_model, test_df[, setdiff(names(test_df), outcome_name), drop = FALSE], outcome_levels, fallback_probs)

      fold_predictions[[fold_id]] <- tibble::tibble(
        sample_id = rownames(test_df),
        fold = fold_id,
        mtry = mtry_value,
        ntree = ntree_value,
        nodesize = nodesize_value,
        actual = test_df[[outcome_name]],
        predicted_class = prediction$predicted_class,
        prob_positive = prediction$probabilities[, spec$positive_class]
      )
    }

    all_predictions[[paste(mtry_value, ntree_value, nodesize_value, sep = "_")]] <- dplyr::bind_rows(fold_predictions)
  }

  all_predictions_df <- dplyr::bind_rows(all_predictions)
  rf_grid_summary <- summarise_rf_grid(all_predictions_df, spec)
  # Limitation: the best setting is chosen on the same folds that give the
  # reported score (no nested cross-validation), so that score is optimistic.
  best_mtry <- rf_grid_summary$mtry[[1]]
  best_ntree <- rf_grid_summary$ntree[[1]]
  best_nodesize <- rf_grid_summary$nodesize[[1]]

  best_predictions <- all_predictions_df |>
    dplyr::filter(mtry == best_mtry, ntree == best_ntree, nodesize == best_nodesize) |>
    dplyr::arrange(sample_id)

  cv_metrics <- compute_classification_metrics(
    actual = best_predictions$actual,
    predicted = best_predictions$predicted_class,
    probability = best_predictions$prob_positive,
    model_name = "Random forest",
    outcome_type = spec$type,
    positive_class = spec$positive_class
  )

  # Refit the best random-forest setting on the full cleaned dataset to obtain
  # final variable importance and the proximity-based MDS plot.
  final_model <- fit_rf_safely(
    full_processed,
    outcome_name,
    mtry = best_mtry,
    ntree = best_ntree,
    nodesize = best_nodesize,
    proximity = TRUE
  )
  importance <- extract_rf_importance(final_model)
  mds_object <- create_rf_mds_plot(
    final_model,
    full_processed,
    outcome_name,
    paste(spec$display_name, "random forest MDS plot")
  )

  list(
    grid = rf_grid_summary,
    best_mtry = best_mtry,
    best_ntree = best_ntree,
    best_nodesize = best_nodesize,
    cv_predictions = best_predictions |>
      dplyr::rename(
        rf_prediction = predicted_class,
        rf_probability = prob_positive
      ),
    cv_metrics = cv_metrics,
    final_model = final_model,
    importance = importance,
    mds = mds_object
  )
}

# Run the multiclass regression workflow for phenotype with guarded scoring and
# class-specific coefficient extraction.
run_multiclass_model_pipeline <- function(dataset_df, spec, fold_plan) {
  outcome_name <- spec$cleaned_outcome
  outcome_levels <- levels(dataset_df[[outcome_name]])

  cv_predictions <- list()
  screening_tables <- list()

  for (fold_id in seq_along(fold_plan$folds)) {
    test_index <- fold_plan$folds[[fold_id]]
    train_index <- setdiff(seq_len(nrow(dataset_df)), test_index)

    train_raw <- dataset_df[train_index, , drop = FALSE]
    test_raw <- dataset_df[test_index, , drop = FALSE]

    processed <- preprocess_train_test(train_raw, outcome_name, outcome_levels, test_raw)
    train_df <- processed$train
    test_df <- processed$test

    selection <- select_screened_predictors(train_df, outcome_name, spec$type, config$max_multiclass_df)
    selected_predictors <- selection$selected_predictors

    model_formula <- if (length(selected_predictors) > 0) {
      stats::reformulate(selected_predictors, response = outcome_name)
    } else {
      stats::as.formula(paste(outcome_name, "~ 1"))
    }

    # Use multinomial regression instead of binary GLM because the phenotype
    # endpoint has more than two classes.
    multinom_model <- fit_multinom_safely(model_formula, train_df)
    class_prior <- prop.table(table(train_df[[outcome_name]]))
    fallback_probs <- setNames(rep(0, length(outcome_levels)), outcome_levels)
    fallback_probs[names(class_prior)] <- as.numeric(class_prior)

    prediction <- predict_with_multinom(
      multinom_model,
      test_df,
      outcome_levels = outcome_levels,
      fallback_probs = fallback_probs
    )

    probability_df <- as.data.frame(prediction$probabilities)
    names(probability_df) <- paste0("multinom_prob_", tolower(names(probability_df)))

    cv_predictions[[fold_id]] <- dplyr::bind_cols(
      tibble::tibble(
        sample_id = rownames(test_df),
        fold = fold_id,
        actual = test_df[[outcome_name]],
        multinom_prediction = prediction$predicted_class
      ),
      probability_df
    )

    screening_tables[[fold_id]] <- selection$score_table |>
      dplyr::mutate(fold = fold_id)
  }

  cv_predictions <- dplyr::bind_rows(cv_predictions) |>
    dplyr::arrange(sample_id)

  # Score only on the indices allowed by the fold plan so singleton classes do
  # not create invalid validation summaries.
  scoring_predictions <- cv_predictions |>
    dplyr::filter(sample_id %in% rownames(dataset_df)[fold_plan$scoring_indices])

  cv_metrics <- compute_classification_metrics(
    actual = scoring_predictions$actual,
    predicted = scoring_predictions$multinom_prediction,
    probability = NULL,
    model_name = "Multinomial model",
    outcome_type = spec$type
  )

  full_processed <- preprocess_train_test(dataset_df, outcome_name, outcome_levels)$train
  full_selection <- select_screened_predictors(full_processed, outcome_name, spec$type, config$max_multiclass_df)

  final_formula <- if (length(full_selection$selected_predictors) > 0) {
    stats::reformulate(full_selection$selected_predictors, response = outcome_name)
  } else {
    stats::as.formula(paste(outcome_name, "~ 1"))
  }

  final_model <- fit_multinom_safely(final_formula, full_processed)
  coefficients <- extract_multinom_coefficients(final_model)
  model_fit <- tibble::tibble(
    model = "Multinomial model",
    formula = paste(deparse(final_formula), collapse = " "),
    aic = if (!is.null(final_model)) stats::AIC(final_model) else NA_real_
  )

  list(
    cv_predictions = cv_predictions,
    scoring_predictions = scoring_predictions,
    cv_metrics = cv_metrics,
    screening = dplyr::bind_rows(screening_tables),
    coefficients = coefficients,
    model_fit = model_fit,
    selected_predictors = full_selection$selected_predictors,
    final_model = final_model
  )
}

# Run the multiclass random-forest workflow, including fold-wise tuning and
# final importance/MDS outputs.
run_multiclass_rf_pipeline <- function(dataset_df, spec, fold_plan) {
  outcome_name <- spec$cleaned_outcome
  outcome_levels <- levels(dataset_df[[outcome_name]])
  full_processed <- preprocess_train_test(dataset_df, outcome_name, outcome_levels)$train
  predictor_count <- ncol(full_processed) - 1
  rf_grid <- make_rf_grid(predictor_count, spec$type)

  all_predictions <- list()

  for (grid_id in seq_len(nrow(rf_grid))) {
    mtry_value <- rf_grid$mtry[[grid_id]]
    ntree_value <- rf_grid$ntree[[grid_id]]
    nodesize_value <- rf_grid$nodesize[[grid_id]]
    fold_predictions <- list()

    for (fold_id in seq_along(fold_plan$folds)) {
      test_index <- fold_plan$folds[[fold_id]]
      train_index <- setdiff(seq_len(nrow(dataset_df)), test_index)

      train_raw <- dataset_df[train_index, , drop = FALSE]
      test_raw <- dataset_df[test_index, , drop = FALSE]

      processed <- preprocess_train_test(train_raw, outcome_name, outcome_levels, test_raw)
      train_df <- processed$train
      test_df <- processed$test

      class_prior <- prop.table(table(train_df[[outcome_name]]))
      fallback_probs <- setNames(rep(0, length(outcome_levels)), outcome_levels)
      fallback_probs[names(class_prior)] <- as.numeric(class_prior)

      rf_model <- fit_rf_safely(
        train_df,
        outcome_name,
        mtry = mtry_value,
        ntree = ntree_value,
        nodesize = nodesize_value,
        proximity = FALSE
      )
      prediction <- predict_with_rf(
        rf_model,
        test_df[, setdiff(names(test_df), outcome_name), drop = FALSE],
        outcome_levels,
        fallback_probs
      )

      # Carry class probabilities forward even in the multiclass setting so the
      # prediction export remains interpretable.
      probability_df <- as.data.frame(prediction$probabilities)
      names(probability_df) <- paste0("rf_prob_", tolower(names(probability_df)))

      fold_predictions[[fold_id]] <- dplyr::bind_cols(
        tibble::tibble(
          sample_id = rownames(test_df),
          fold = fold_id,
          mtry = mtry_value,
          ntree = ntree_value,
          nodesize = nodesize_value,
          actual = test_df[[outcome_name]],
          predicted_class = prediction$predicted_class
        ),
        probability_df
      )
    }

    all_predictions[[paste(mtry_value, ntree_value, nodesize_value, sep = "_")]] <- dplyr::bind_rows(fold_predictions)
  }

  all_predictions_df <- dplyr::bind_rows(all_predictions)
  scoring_predictions <- all_predictions_df |>
    dplyr::filter(sample_id %in% rownames(dataset_df)[fold_plan$scoring_indices])

  grid_summary <- summarise_rf_grid(scoring_predictions, spec)

  # Same limitation as the binary workflow: tuning is not nested.
  best_mtry <- grid_summary$mtry[[1]]
  best_ntree <- grid_summary$ntree[[1]]
  best_nodesize <- grid_summary$nodesize[[1]]
  best_predictions <- all_predictions_df |>
    dplyr::filter(mtry == best_mtry, ntree == best_ntree, nodesize == best_nodesize) |>
    dplyr::arrange(sample_id)

  best_scoring_predictions <- best_predictions |>
    dplyr::filter(sample_id %in% rownames(dataset_df)[fold_plan$scoring_indices])

  cv_metrics <- compute_classification_metrics(
    actual = best_scoring_predictions$actual,
    predicted = best_scoring_predictions$predicted_class,
    probability = NULL,
    model_name = "Random forest",
    outcome_type = spec$type
  )

  final_model <- fit_rf_safely(
    full_processed,
    outcome_name,
    mtry = best_mtry,
    ntree = best_ntree,
    nodesize = best_nodesize,
    proximity = TRUE
  )
  importance <- extract_rf_importance(final_model)
  mds_object <- create_rf_mds_plot(
    final_model,
    full_processed,
    outcome_name,
    paste(spec$display_name, "random forest MDS plot")
  )

  list(
    grid = grid_summary,
    best_mtry = best_mtry,
    best_ntree = best_ntree,
    best_nodesize = best_nodesize,
    cv_predictions = best_predictions |>
      dplyr::rename(rf_prediction = predicted_class),
    scoring_predictions = best_scoring_predictions |>
      dplyr::rename(rf_prediction = predicted_class),
    cv_metrics = cv_metrics,
    final_model = final_model,
    importance = importance,
    mds = mds_object
  )
}
