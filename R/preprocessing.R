# Reading the data, cleaning it, and building the cross-validation folds.

# `x %||% y` means: use y if x is NULL, otherwise use x.
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# Make column names safe to use in formulas: lowercase, with underscores.
clean_column_names <- function(x) {
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  make.unique(x, sep = "_")
}

# Prefer the known analysis sheet names first, but still fall back safely
# if the workbook structure changes in the future.
resolve_sheet_name <- function(path, preferred_sheet) {
  sheets <- readxl::excel_sheets(path)
  if (length(sheets) == 0) {
    stop("No worksheets were found in the workbook.")
  }

  matching_sheet <- preferred_sheet[preferred_sheet %in% sheets]
  if (length(matching_sheet) > 0) {
    return(matching_sheet[[1]])
  }

  if (length(sheets) == 1) {
    return(sheets[[1]])
  }

  warning(
    sprintf(
      "Preferred sheet '%s' was not found. Falling back to '%s'.",
      preferred_sheet,
      sheets[[1]]
    ),
    call. = FALSE
  )

  sheets[[1]]
}

# Recode raw outcome values into the requested dissertation labels while
# keeping factor levels explicit and stable across later modelling steps.
normalize_outcome_by_spec <- function(x, spec) {
  x_chr <- trimws(tolower(as.character(x)))
  x_chr[x_chr %in% c("", "na", "n/a", "nan", "null")] <- NA_character_

  map_names <- tolower(names(spec$label_map))
  map_values <- unname(spec$label_map)
  normalized <- map_values[match(x_chr, map_names)]
  normalized <- as.character(normalized)
  normalized[is.na(normalized)] <- NA_character_

  if (spec$type == "binary") {
    levels_out <- c(spec$negative_class, spec$positive_class)
  } else {
    levels_out <- unname(spec$label_map)
  }

  factor(normalized, levels = levels_out)
}

# Read the workbook once, clean the structure, and align the real column
# names to the standard names used in the combined pipeline.
read_raw_data <- function(path, preferred_sheet) {
  if (!file.exists(path)) {
    stop("Input workbook not found: ", path)
  }

  selected_sheet <- resolve_sheet_name(path, preferred_sheet)
  raw_df <- readxl::read_excel(path, sheet = selected_sheet) |>
    as.data.frame(check.names = FALSE)

  original_names <- names(raw_df)
  cleaned_names <- clean_column_names(original_names)
  names(raw_df) <- cleaned_names

  if (!config$sample_id %in% names(raw_df)) {
    stop("The workbook must contain 'Sample ID'.")
  }

  # Use sample identifiers as row names so fold predictions can be traced back
  # to isolates without keeping the ID as a predictor.
  rownames(raw_df) <- make.unique(as.character(raw_df[[config$sample_id]]))
  raw_df[[config$sample_id]] <- NULL

  # Trim character input and convert blank strings to missing values before any
  # outcome-specific filtering is performed.
  for (nm in names(raw_df)) {
    if (is.character(raw_df[[nm]])) {
      raw_df[[nm]] <- trimws(raw_df[[nm]])
      raw_df[[nm]][raw_df[[nm]] == ""] <- NA_character_
    } else if (is.logical(raw_df[[nm]])) {
      raw_df[[nm]] <- factor(ifelse(raw_df[[nm]], "TRUE", "FALSE"))
    }
  }

  # Map the real workbook outcome columns onto the standard outcome names used
  # throughout the rest of the script.
  for (spec in outcome_specs) {
    candidates <- clean_column_names(spec$raw_column_candidates)
    match_name <- intersect(candidates, names(raw_df))
    if (length(match_name) == 0) {
      stop("Outcome column for ", spec$display_name, " was not found.")
    }
    if (spec$cleaned_outcome != match_name[[1]]) {
      names(raw_df)[names(raw_df) == match_name[[1]]] <- spec$cleaned_outcome
    }
  }

  list(
    data = raw_df,
    sheet = selected_sheet,
    original_names = original_names
  )
}

# Build one outcome-specific dataset at a time so missing ventilation values do
# not unnecessarily remove rows from gentamicin or phenotype analyses.
prepare_outcome_dataset <- function(raw_df, spec, all_outcomes) {
  analysis_df <- raw_df
  analysis_df[[spec$cleaned_outcome]] <- normalize_outcome_by_spec(analysis_df[[spec$cleaned_outcome]], spec)

  # Count invalid but non-missing codes separately so exclusions can be reported
  # transparently in the workbook notes.
  invalid_mask <- is.na(analysis_df[[spec$cleaned_outcome]]) & !is.na(raw_df[[spec$cleaned_outcome]])
  invalid_count <- sum(invalid_mask)

  analysis_df <- analysis_df[!is.na(analysis_df[[spec$cleaned_outcome]]), , drop = FALSE]

  if (nrow(analysis_df) == 0) {
    stop("No rows remain after filtering valid outcomes for ", spec$display_name, ".")
  }

  # Remove the other outcome columns from the current analysis set to reduce
  # leakage and keep each endpoint analytically separate.
  keep_columns <- setdiff(names(analysis_df), setdiff(all_outcomes, spec$cleaned_outcome))
  analysis_df <- analysis_df[, keep_columns, drop = FALSE]

  for (nm in names(analysis_df)) {
    if (nm == spec$cleaned_outcome) {
      next
    }
    if (is.character(analysis_df[[nm]])) {
      analysis_df[[nm]] <- as.factor(analysis_df[[nm]])
    }
  }

  list(
    data = analysis_df,
    invalid_count = invalid_count,
    n_rows = nrow(analysis_df),
    class_counts = table(analysis_df[[spec$cleaned_outcome]])
  )
}

# Apply the same training-only preprocessing rules to every model so missing
# data, unseen categories, and unstable predictors are handled consistently.
preprocess_train_test <- function(train_df, outcome_name, outcome_levels, test_df = NULL) {
  train_processed <- train_df
  test_processed <- test_df
  predictor_names <- setdiff(names(train_processed), outcome_name)

  for (nm in predictor_names) {
    train_x <- train_processed[[nm]]
    test_x <- if (!is.null(test_processed)) test_processed[[nm]] else NULL

    # Treat categorical predictors by keeping missing values explicit and
    # reserving a "novel" level for unseen categories at prediction time.
    if (is.factor(train_x) || is.character(train_x)) {
      train_chr <- trimws(as.character(train_x))
      train_chr[is.na(train_chr) | train_chr == ""] <- "missing"

      training_levels <- unique(train_chr)
      final_levels <- unique(c(training_levels, "missing", "novel"))
      train_processed[[nm]] <- factor(train_chr, levels = final_levels)

      if (!is.null(test_processed)) {
        test_chr <- trimws(as.character(test_x))
        test_chr[is.na(test_chr) | test_chr == ""] <- "missing"
        test_chr[!(test_chr %in% training_levels) & test_chr != "missing"] <- "novel"
        test_processed[[nm]] <- factor(test_chr, levels = final_levels)
      }
    # Impute numeric predictors with the training-set median so no test-set
    # information leaks back into the fit.
    } else if (is.numeric(train_x) || is.integer(train_x)) {
      train_median <- stats::median(train_x, na.rm = TRUE)
      if (!is.finite(train_median)) {
        train_median <- 0
      }

      train_x[is.na(train_x)] <- train_median
      train_processed[[nm]] <- train_x

      if (!is.null(test_processed)) {
        test_x[is.na(test_x)] <- train_median
        test_processed[[nm]] <- test_x
      }
    }
  }

  # Remove constant predictors before modelling because they add no signal and
  # can trigger factor/contrast errors in small datasets.
  keep_non_constant <- predictor_names[
    vapply(
      train_processed[predictor_names],
      function(x) length(unique(stats::na.omit(x))) >= 2,
      logical(1)
    )
  ]

  train_processed <- train_processed[, c(outcome_name, keep_non_constant), drop = FALSE]
  if (!is.null(test_processed)) {
    test_processed <- test_processed[, c(outcome_name, keep_non_constant), drop = FALSE]
  }

  # Remove near-zero variance predictors to reduce instability and noise,
  # especially important in p >> n biomedical feature sets.
  if (length(keep_non_constant) > 0) {
    nzv_metrics <- caret::nearZeroVar(
      train_processed[, keep_non_constant, drop = FALSE],
      saveMetrics = TRUE
    )
    keep_nzv <- rownames(nzv_metrics)[!nzv_metrics$nzv]
  } else {
    keep_nzv <- character(0)
  }

  train_processed <- train_processed[, c(outcome_name, keep_nzv), drop = FALSE]
  if (!is.null(test_processed)) {
    test_processed <- test_processed[, c(outcome_name, keep_nzv), drop = FALSE]
  }

  train_processed[[outcome_name]] <- factor(train_processed[[outcome_name]], levels = outcome_levels)
  if (!is.null(test_processed)) {
    test_processed[[outcome_name]] <- factor(test_processed[[outcome_name]], levels = outcome_levels)
  }

  list(
    train = train_processed,
    test = test_processed,
    predictor_names = keep_nzv
  )
}

# Screen predictors univariately before multivariable modelling so the final
# regression models stay small enough to be estimable on a very small sample.
score_predictors_univariately <- function(train_df, outcome_name, outcome_type) {
  outcome <- train_df[[outcome_name]]
  predictors <- setdiff(names(train_df), outcome_name)

  if (length(predictors) == 0) {
    return(tibble::tibble())
  }

  # Keep only predictors with at least minimal observed variation so sparse
  # one-off categories do not dominate screening in tiny folds.
  # vapply checks each column in turn and must return one TRUE or FALSE for it.
  stable_predictors <- predictors[
    vapply(train_df[predictors], function(x) {
      observed <- stats::na.omit(x)

      if (length(unique(observed)) < 2) {
        return(FALSE)
      }

      if (is.factor(x) || length(unique(observed)) <= 5) {
        return(all(table(observed) >= 2))
      }

      TRUE
    }, logical(1))
  ]

  if (length(stable_predictors) == 0) {
    return(tibble::tibble())
  }

  # lapply runs this function once per predictor and gives a list of one-row
  # tables; bind_rows stacks them into a single table.
  dplyr::bind_rows(lapply(stable_predictors, function(pred_name) {
    x <- train_df[[pred_name]]

    # Choose a simple test that suits the predictor and outcome type. If a test
    # fails, return a p-value of 1, which puts the predictor at the bottom of the
    # ranking so it is effectively ignored.
    p_value <- tryCatch(
      {
        if (outcome_type == "binary") {
          if (is.factor(x) || length(unique(x)) <= 5) {
            contingency <- table(x, outcome)
            if (nrow(contingency) < 2 || ncol(contingency) < 2) {
              1
            } else if (all(dim(contingency) == c(2, 2))) {
              stats::fisher.test(contingency)$p.value
            } else {
              suppressWarnings(stats::chisq.test(contingency, simulate.p.value = TRUE, B = 2000)$p.value)
            }
          } else {
            suppressWarnings(stats::wilcox.test(x ~ outcome)$p.value)
          }
        } else {
          if (is.factor(x) || length(unique(x)) <= 5) {
            contingency <- table(x, outcome)
            if (nrow(contingency) < 2 || ncol(contingency) < 2) {
              1
            } else {
              suppressWarnings(stats::chisq.test(contingency, simulate.p.value = TRUE, B = 2000)$p.value)
            }
          } else {
            suppressWarnings(stats::kruskal.test(x ~ outcome)$p.value)
          }
        }
      },
      error = function(e) 1
    )

    tibble::tibble(predictor = pred_name, screening_p_value = p_value)
  })) |>
    dplyr::arrange(screening_p_value, predictor)
}

# Convert screening results into a constrained predictor set by checking the
# effective degrees of freedom of the candidate design matrix.
select_screened_predictors <- function(train_df, outcome_name, outcome_type, max_df_limit) {
  score_table <- score_predictors_univariately(train_df, outcome_name, outcome_type)

  if (nrow(score_table) == 0) {
    return(list(selected_predictors = character(0), score_table = tibble::tibble(), max_df = 0L))
  }

  if (outcome_type == "binary") {
    outcome_counts <- table(train_df[[outcome_name]])
    max_df <- max(1L, min(max_df_limit, floor(min(outcome_counts) / 2)))
  } else {
    max_df <- max_df_limit
  }

  selected_predictors <- character(0)

  for (pred_name in score_table$predictor) {
    candidate_predictors <- c(selected_predictors, pred_name)
    candidate_formula <- stats::reformulate(candidate_predictors, response = outcome_name)
    design_matrix <- tryCatch(stats::model.matrix(candidate_formula, data = train_df), error = function(e) NULL)

    if (is.null(design_matrix)) {
      next
    }

    # Degrees of freedom are counted from the model matrix because a categorical
    # predictor uses one column per level after the first, not just one column.
    effective_df <- ncol(design_matrix) - 1
    if (effective_df <= max_df) {
      selected_predictors <- candidate_predictors
    }
  }

  if (length(selected_predictors) == 0) {
    for (pred_name in score_table$predictor) {
      single_formula <- stats::reformulate(pred_name, response = outcome_name)
      design_matrix <- tryCatch(stats::model.matrix(single_formula, data = train_df), error = function(e) NULL)
      if (!is.null(design_matrix) && (ncol(design_matrix) - 1) <= max_df) {
        selected_predictors <- pred_name
        break
      }
    }
  }

  list(
    selected_predictors = selected_predictors,
    score_table = dplyr::mutate(
      score_table,
      selected_for_model = predictor %in% selected_predictors
    ),
    max_df = max_df
  )
}

# Build a fold plan that adapts to class imbalance and extreme rarity so the
# phenotype workflow remains valid when one class has very few isolates.
make_fold_plan <- function(y, requested_k, outcome_type) {
  y <- droplevels(factor(y))
  class_counts <- table(y)

  # In the multiclass setting, score only classes that are sufficiently
  # represented while still fitting the final model on the full cleaned dataset.
  if (outcome_type == "multiclass" && any(class_counts < 2)) {
    common_classes <- names(class_counts)[class_counts >= 2]
    singleton_classes <- names(class_counts)[class_counts < 2]

    if (length(common_classes) == 0) {
      folds <- lapply(seq_along(y), function(i) i)
      return(list(
        folds = folds,
        scoring_indices = seq_along(y),
        method_label = "Leave-one-out CV",
        validation_note = "All classes were too sparse for stratified k-fold CV."
      ))
    }

    scoring_indices <- which(y %in% common_classes)
    scoring_y <- droplevels(y[scoring_indices])
    k <- min(requested_k, min(table(scoring_y)))
    base_folds <- caret::createFolds(y = scoring_y, k = k, returnTrain = FALSE)

    folds <- lapply(base_folds, function(idx) scoring_indices[idx])

    return(list(
      folds = folds,
      scoring_indices = scoring_indices,
      method_label = sprintf("%d-fold CV on classes with >=2 samples", k),
      validation_note = paste(
        "Rare classes were retained for model fitting on the full dataset but were excluded from fold-based scoring:",
        paste(singleton_classes, collapse = ", ")
      )
    ))
  }

  min_count <- min(class_counts)
  if (min_count < 2) {
    folds <- lapply(seq_along(y), function(i) i)
    return(list(
      folds = folds,
      scoring_indices = seq_along(y),
      method_label = "Leave-one-out CV",
      validation_note = "Class counts were too small for stratified k-fold CV."
    ))
  }

  k <- min(requested_k, min_count)
  folds <- caret::createFolds(y = y, k = k, returnTrain = FALSE)

  list(
    folds = folds,
    scoring_indices = seq_along(y),
    method_label = sprintf("%d-fold stratified CV", k),
    validation_note = NA_character_
  )
}
