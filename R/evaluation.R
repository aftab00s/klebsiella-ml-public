# Performance metrics, plots and the results workbook.

# Prevent exact 0 and 1 probabilities from destabilising downstream ROC work.
clip_probabilities <- function(x, lower = 1e-06, upper = 1 - 1e-06) {
  pmin(pmax(as.numeric(x), lower), upper)
}

# Calculate a p-value against the no-information rate so classification
# performance can be described in a more conventional reporting style.
compute_nir_p_value <- function(actual, predicted) {
  actual <- factor(actual)
  predicted <- factor(predicted, levels = levels(actual))

  if (length(actual) == 0) {
    return(NA_real_)
  }

  no_information_rate <- max(prop.table(table(actual)))

  tryCatch(
    stats::binom.test(
      x = sum(actual == predicted, na.rm = TRUE),
      n = length(actual),
      p = no_information_rate,
      alternative = "greater"
    )$p.value,
    error = function(e) NA_real_
  )
}

# Compute binary and multiclass metrics through one shared interface so later
# performance summaries can be combined cleanly across outcomes.
compute_classification_metrics <- function(actual, predicted, probability = NULL, model_name, outcome_type, positive_class = NULL) {
  actual <- factor(actual)
  predicted <- factor(predicted, levels = levels(actual))

  confusion <- if (outcome_type == "binary") {
    caret::confusionMatrix(predicted, actual, positive = positive_class)
  } else {
    caret::confusionMatrix(predicted, actual)
  }

  confusion_table <- as.data.frame.matrix(confusion$table) |>
    tibble::rownames_to_column("predicted_class")

  # Use the standard binary confusion-matrix path when the endpoint is binary,
  # including ROC/AUC when probabilities are available.
  if (outcome_type == "binary") {
    sensitivity <- unname(confusion$byClass["Sensitivity"])
    specificity <- unname(confusion$byClass["Specificity"])
    balanced_accuracy <- unname(confusion$byClass["Balanced Accuracy"])
    false_positive_rate <- 1 - specificity

    roc_object <- NULL
    auc_value <- NA_real_
    if (!is.null(probability)) {
      # No `direction` is set, so pROC chooses it by comparing the median prediction
      # in each class. A near-chance model can therefore have its direction flipped.
      # This is left as it was in the original analysis.
      roc_object <- pROC::roc(
        response = actual,
        predictor = clip_probabilities(probability),
        levels = levels(actual),
        quiet = TRUE
      )
      auc_value <- as.numeric(pROC::auc(roc_object))
    }

    classwise <- tibble::tibble(
      class = positive_class,
      sensitivity = sensitivity,
      specificity = specificity,
      balanced_accuracy = balanced_accuracy,
      false_positive_rate = false_positive_rate
    )

    summary_table <- tibble::tibble(
      model = model_name,
      accuracy = unname(confusion$overall["Accuracy"]),
      kappa = unname(confusion$overall["Kappa"]),
      sensitivity = sensitivity,
      specificity = specificity,
      balanced_accuracy = balanced_accuracy,
      false_positive_rate = false_positive_rate,
      auc = auc_value,
      p_value_vs_nir = compute_nir_p_value(actual, predicted)
    )

    return(list(
      summary = summary_table,
      confusion = confusion_table,
      classwise = classwise,
      roc = roc_object
    ))
  }

  # For multiclass outcomes, calculate one-vs-rest summaries manually because
  # confusionMatrix output becomes less stable when scored folds are sparse.
  cm <- as.matrix(confusion$table)
  all_classes <- colnames(cm)

  by_class <- dplyr::bind_rows(lapply(all_classes, function(class_name) {
    tp <- cm[class_name, class_name]
    fn <- sum(cm[, class_name]) - tp
    fp <- sum(cm[class_name, ]) - tp
    tn <- sum(cm) - tp - fn - fp

    sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)

    tibble::tibble(
      class = class_name,
      sensitivity = sensitivity,
      specificity = specificity,
      balanced_accuracy = balanced_accuracy,
      false_positive_rate = 1 - specificity
    )
  }))

  summary_table <- tibble::tibble(
    model = model_name,
    accuracy = unname(confusion$overall["Accuracy"]),
    kappa = unname(confusion$overall["Kappa"]),
    sensitivity = mean(by_class$sensitivity, na.rm = TRUE),
    specificity = mean(by_class$specificity, na.rm = TRUE),
    balanced_accuracy = mean(by_class$balanced_accuracy, na.rm = TRUE),
    false_positive_rate = mean(by_class$false_positive_rate, na.rm = TRUE),
    auc = NA_real_,
    p_value_vs_nir = compute_nir_p_value(actual, predicted)
  )

  list(
    summary = summary_table,
    confusion = confusion_table,
    classwise = by_class,
    roc = NULL
  )
}

# Convert ROC objects to plotting data frames for consistent ggplot output.
roc_to_df <- function(roc_object, model_name) {
  if (is.null(roc_object)) {
    return(tibble::tibble())
  }

  coords <- pROC::coords(
    roc_object,
    x = "all",
    ret = c("specificity", "sensitivity"),
    transpose = FALSE
  )

  tibble::tibble(
    model = model_name,
    specificity = coords$specificity,
    sensitivity = coords$sensitivity,
    false_positive_rate = 1 - coords$specificity
  )
}

# Bar chart of the top predictors by variable importance.
create_rf_importance_plot <- function(importance_df, title_text) {
  importance_df |>
    dplyr::slice_head(n = 20) |>
    ggplot2::ggplot(ggplot2::aes(x = reorder(predictor, importance_score), y = importance_score)) +
    ggplot2::geom_col(fill = "#2C7FB8") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title_text, x = NULL, y = "Random forest importance")
}

# Use the proximity matrix to create an MDS plot, which gives a visual summary
# of sample similarity under the random-forest model.
create_rf_mds_plot <- function(model, full_processed_df, outcome_name, title_text) {
  if (is.null(model) || is.null(model$proximity)) {
    return(NULL)
  }

  mds_points <- stats::cmdscale(1 - model$proximity, k = 2)
  mds_df <- tibble::tibble(
    sample_id = rownames(full_processed_df),
    dim1 = mds_points[, 1],
    dim2 = mds_points[, 2],
    outcome = full_processed_df[[outcome_name]]
  )

  plot_object <- ggplot2::ggplot(
    mds_df,
    ggplot2::aes(x = dim1, y = dim2, color = outcome, label = sample_id)
  ) +
    # Add a light class-level confidence ellipse as a visual grouping boundary.
    # Ellipses are drawn only for classes with enough points for stable shape
    # estimation; sparse classes remain visible as labelled points.
    ggplot2::stat_ellipse(
      ggplot2::aes(fill = outcome, group = outcome),
      geom = "polygon",
      type = "norm",
      level = 0.68,
      alpha = 0.12,
      color = NA,
      show.legend = FALSE
    ) +
    ggplot2::stat_ellipse(
      ggplot2::aes(group = outcome),
      type = "norm",
      level = 0.68,
      linewidth = 0.7,
      linetype = "dashed",
      show.legend = FALSE
    ) +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    # Label isolates by the cleaned Sample ID row names; ggrepel keeps labels
    # readable when isolates cluster tightly in the MDS space.
    ggrepel::geom_text_repel(size = 3, max.overlaps = Inf, show.legend = FALSE) +
    ggplot2::labs(
      title = title_text,
      x = "Dimension 1",
      y = "Dimension 2",
      color = "Outcome"
    )

  list(data = mds_df, plot = plot_object)
}

# Write each result table to a workbook sheet and optionally add a distant note
# linking the sheet to a related plot in the PDF.
write_sheet_with_note <- function(wb, sheet_name, data, note = NULL) {
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet = sheet_name, x = data)
  openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)
  openxlsx::setColWidths(wb, sheet = sheet_name, cols = 1:min(20, max(1, ncol(data))), widths = "auto")

  if (!is.null(note) && !is.na(note) && nzchar(note)) {
    openxlsx::writeData(wb, sheet = sheet_name, x = note, startCol = 52, startRow = 1, colNames = FALSE)
  }
}

# Save all plots into one PDF so the final output is easy to archive and cite.
save_plot_pdf <- function(plot_list, output_path) {
  grDevices::pdf(output_path, width = 10, height = 7)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (plot_name in names(plot_list)) {
    plot_object <- plot_list[[plot_name]]
    if (!is.null(plot_object)) {
      print(plot_object)
    }
  }
}

# Build separate ROC plots for the two binary outcomes and a comparison plot
# that places GLM and random forest on the same axes.
create_binary_roc_plots <- function(glm_metrics, rf_metrics, title_prefix) {
  glm_roc_df <- roc_to_df(glm_metrics$roc, "GLM")
  rf_roc_df <- roc_to_df(rf_metrics$roc, "Random forest")
  glm_auc <- glm_metrics$summary$auc[[1]]
  rf_auc <- rf_metrics$summary$auc[[1]]

  glm_plot <- ggplot2::ggplot(glm_roc_df, ggplot2::aes(x = false_positive_rate, y = sensitivity)) +
    ggplot2::geom_path(color = "#C0392B", linewidth = 1.1) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey55") +
    ggplot2::annotate(
      "text",
      x = 0.65,
      y = 0.15,
      label = sprintf("AUC = %.3f", glm_auc),
      hjust = 0,
      size = 4
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste(title_prefix, "GLM ROC"),
      x = "False positive rate",
      y = "True positive rate"
    )

  rf_plot <- ggplot2::ggplot(rf_roc_df, ggplot2::aes(x = false_positive_rate, y = sensitivity)) +
    ggplot2::geom_path(color = "#2C7FB8", linewidth = 1.1) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey55") +
    ggplot2::annotate(
      "text",
      x = 0.65,
      y = 0.15,
      label = sprintf("AUC = %.3f", rf_auc),
      hjust = 0,
      size = 4
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste(title_prefix, "random forest ROC"),
      x = "False positive rate",
      y = "True positive rate"
    )

  comparison_df <- dplyr::bind_rows(glm_roc_df, rf_roc_df)
  comparison_plot <- ggplot2::ggplot(
    comparison_df,
    ggplot2::aes(x = false_positive_rate, y = sensitivity, color = model)
  ) +
    ggplot2::geom_path(linewidth = 1.1) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey55") +
    ggplot2::annotate(
      "text",
      x = 0.55,
      y = 0.15,
      label = sprintf("GLM AUC = %.3f\nRF AUC = %.3f", glm_auc, rf_auc),
      hjust = 0,
      size = 4
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste(title_prefix, "ROC comparison"),
      x = "False positive rate",
      y = "True positive rate",
      color = "Model"
    )

  list(glm_plot = glm_plot, rf_plot = rf_plot, comparison_plot = comparison_plot)
}
