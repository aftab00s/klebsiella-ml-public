# Step 3: summarise performance, draw the plots and write the results workbook.

source(here::here("R", "config.R"))
source(here::here("R", "preprocessing.R"))
source(here::here("R", "model_fitting.R"))
source(here::here("R", "evaluation.R"))

prepared <- readRDS(file.path(config$output_dir, "prepared_datasets.rds"))
prepared_datasets <- prepared$prepared_datasets
fits <- readRDS(file.path(config$output_dir, "fitted_models.rds"))

fold_plans <- fits$fold_plans
mv_glm <- fits$mv_glm
mv_rf <- fits$mv_rf
mv_glm_model_comparison <- fits$mv_glm_model_comparison
gent_glm <- fits$gent_glm
gent_rf <- fits$gent_rf
gent_glm_model_comparison <- fits$gent_glm_model_comparison
phenotype_model <- fits$phenotype_model
phenotype_rf <- fits$phenotype_rf

# -----------------------------
# Performance tables
# -----------------------------

mv_predictions <- dplyr::full_join(
  mv_glm$cv_predictions,
  mv_rf$cv_predictions,
  by = c("sample_id", "fold", "actual")
) |>
  dplyr::arrange(sample_id)

mv_performance_summary <- dplyr::bind_rows(
  mv_glm$cv_metrics$summary,
  mv_rf$cv_metrics$summary
) |>
  dplyr::mutate(
    outcome = outcome_specs$mv$display_name,
    n_samples = prepared_datasets$mv$n_rows,
    resampling = fold_plans$mv$method_label,
    note = fold_plans$mv$validation_note %||% ""
  )

gent_predictions <- dplyr::full_join(
  gent_glm$cv_predictions,
  gent_rf$cv_predictions,
  by = c("sample_id", "fold", "actual")
) |>
  dplyr::arrange(sample_id)

gent_performance_summary <- dplyr::bind_rows(
  gent_glm$cv_metrics$summary,
  gent_rf$cv_metrics$summary
) |>
  dplyr::mutate(
    outcome = outcome_specs$gent$display_name,
    n_samples = prepared_datasets$gent$n_rows,
    resampling = fold_plans$gent$method_label,
    note = fold_plans$gent$validation_note %||% ""
  )

all_glm_model_comparison <- dplyr::bind_rows(
  mv_glm_model_comparison,
  gent_glm_model_comparison
)

# Keep a clearly named duplicate of the detailed GLM comparison so it is easy
# to find in the workbook without confusing it with the main GLM-vs-RF summary.
glm_variant_performance_summary <- all_glm_model_comparison |>
  dplyr::arrange(outcome, model_name)

phenotype_predictions <- dplyr::full_join(
  phenotype_model$cv_predictions,
  phenotype_rf$cv_predictions,
  by = c("sample_id", "fold", "actual")
) |>
  dplyr::arrange(sample_id)

phenotype_performance_summary <- dplyr::bind_rows(
  phenotype_model$cv_metrics$summary,
  phenotype_rf$cv_metrics$summary
) |>
  dplyr::mutate(
    outcome = outcome_specs$phenotype$display_name,
    n_samples = prepared_datasets$phenotype$n_rows,
    resampling = fold_plans$phenotype$method_label,
    note = dplyr::coalesce(fold_plans$phenotype$validation_note, "")
  )

# Combine the outcome-specific summaries into one table for quick comparison
# across endpoints and model families.
combined_model_summary <- dplyr::bind_rows(
  mv_performance_summary,
  gent_performance_summary,
  phenotype_performance_summary
) |>
  dplyr::select(
    outcome,
    model,
    accuracy,
    sensitivity,
    specificity,
    balanced_accuracy,
    false_positive_rate,
    kappa,
    auc,
    p_value_vs_nir,
    n_samples,
    resampling,
    note
  )

# Record the most important methodological notes in a plain table so the final
# workbook documents key assumptions and exclusions.
data_notes <- tibble::tibble(
  item = c(
    "Workbook used",
    "Worksheet used",
    "Mechanical ventilation valid rows",
    "Gentamicin valid rows",
    "Resistance phenotype valid rows",
    "Resistance phenotype invalid codes excluded",
    "Predictor missingness strategy",
    "Binary model families",
    "Phenotype model families",
    "Random forest tuning",
    "Final test-data status",
    "Validation note for phenotype"
  ),
  value = c(
    basename(config$input_path),
    prepared$sheet_used,
    prepared_datasets$mv$n_rows,
    prepared_datasets$gent$n_rows,
    prepared_datasets$phenotype$n_rows,
    prepared_datasets$phenotype$invalid_count,
    "Numeric predictors: median imputation; categorical predictors: explicit missing level plus novel-level handling.",
    "Screened GLM with binomial family and bidirectional stepwise refinement; random forest classifier.",
    "Screened multinomial logistic regression; multiclass random forest classifier.",
    "Expanded mtry, ntree, and nodesize grid search with inverse-frequency class weights and shared cross-validation folds within each outcome.",
    "No separate held-out final test set is defined; only cross-validation is available.",
    fold_plans$phenotype$validation_note %||% "None"
  )
)

# -----------------------------
# Plots
# -----------------------------

# Generate all plots only after the model objects and summaries are finalised,
# which keeps the plotting code separate from the modelling code.
mv_roc_plots <- create_binary_roc_plots(mv_glm$cv_metrics, mv_rf$cv_metrics, "Mechanical ventilation")
gent_roc_plots <- create_binary_roc_plots(gent_glm$cv_metrics, gent_rf$cv_metrics, "Gentamicin resistance")

phenotype_model_comparison_plot <- ggplot2::ggplot(
  phenotype_performance_summary,
  ggplot2::aes(x = model, y = accuracy, fill = model)
) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::labs(
    title = "Resistance phenotype model comparison",
    x = NULL,
    y = "Accuracy"
  ) +
  ggplot2::theme(legend.position = "none")

combined_performance_plot <- ggplot2::ggplot(
  combined_model_summary,
  ggplot2::aes(x = outcome, y = balanced_accuracy, fill = model)
) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::labs(
    title = "Balanced accuracy across outcomes",
    x = NULL,
    y = "Balanced accuracy / macro-balanced accuracy"
  )

# Keep all plots in one named list so one export call can write a complete PDF
# with one page per figure.
plot_list <- list(
  mv_glm_roc = mv_roc_plots$glm_plot,
  mv_rf_roc = mv_roc_plots$rf_plot,
  mv_roc_comparison = mv_roc_plots$comparison_plot,
  mv_rf_variable_importance = create_rf_importance_plot(mv_rf$importance, "Mechanical ventilation random forest importance"),
  mv_rf_mds = mv_rf$mds$plot %||% NULL,
  gent_glm_roc = gent_roc_plots$glm_plot,
  gent_rf_roc = gent_roc_plots$rf_plot,
  gent_roc_comparison = gent_roc_plots$comparison_plot,
  gent_rf_variable_importance = create_rf_importance_plot(gent_rf$importance, "Gentamicin random forest importance"),
  gent_rf_mds = gent_rf$mds$plot %||% NULL,
  phenotype_model_comparison = phenotype_model_comparison_plot,
  phenotype_rf_variable_importance = create_rf_importance_plot(phenotype_rf$importance, "Resistance phenotype random forest importance"),
  phenotype_rf_mds = phenotype_rf$mds$plot %||% NULL,
  combined_balanced_accuracy = combined_performance_plot
)

# -----------------------------
# Export
# -----------------------------

workbook_path <- file.path(config$output_dir, config$workbook_name)
plot_pdf_path <- file.path(config$output_dir, config$plot_pdf_name)

# Write outcome-specific sheets first, then cross-outcome summaries and notes,
# so the workbook mirrors the structure of the script itself.
wb <- openxlsx::createWorkbook()

write_sheet_with_note(wb, "glm_variant_performance_summary", glm_variant_performance_summary, NULL)

write_sheet_with_note(wb, "mv_performance_summary", mv_performance_summary, "Related plot: mv_roc_comparison")
write_sheet_with_note(wb, "mv_predictions", mv_predictions, NULL)
write_sheet_with_note(wb, "mv_glm_coefficients", mv_glm$coefficients, NULL)
write_sheet_with_note(wb, "mv_glm_model_comparison", mv_glm_model_comparison, NULL)
write_sheet_with_note(wb, "mv_rf_grid", mv_rf$grid, NULL)
write_sheet_with_note(wb, "mv_rf_importance", mv_rf$importance, "Related plot: mv_rf_variable_importance")

write_sheet_with_note(wb, "gent_performance_summary", gent_performance_summary, "Related plot: gent_roc_comparison")
write_sheet_with_note(wb, "gent_predictions", gent_predictions, NULL)
write_sheet_with_note(wb, "gent_glm_coefficients", gent_glm$coefficients, NULL)
write_sheet_with_note(wb, "gent_glm_model_comparison", gent_glm_model_comparison, NULL)
write_sheet_with_note(wb, "gent_rf_grid", gent_rf$grid, NULL)
write_sheet_with_note(wb, "gent_rf_importance", gent_rf$importance, "Related plot: gent_rf_variable_importance")

write_sheet_with_note(wb, "phenotype_performance_summary", phenotype_performance_summary, "Related plot: phenotype_model_comparison")
write_sheet_with_note(wb, "phenotype_predictions", phenotype_predictions, NULL)
write_sheet_with_note(wb, "phenotype_model_coefficients", phenotype_model$coefficients, NULL)
write_sheet_with_note(wb, "phenotype_rf_grid", phenotype_rf$grid, NULL)
write_sheet_with_note(wb, "phenotype_rf_importance", phenotype_rf$importance, "Related plot: phenotype_rf_variable_importance")

write_sheet_with_note(wb, "combined_model_summary", combined_model_summary, "Related plot: combined_balanced_accuracy")
write_sheet_with_note(wb, "all_glm_model_comparison", all_glm_model_comparison, NULL)
write_sheet_with_note(wb, "data_dictionary_or_notes", data_notes, NULL)

openxlsx::saveWorkbook(wb, workbook_path, overwrite = TRUE)

# Export every plot into a single PDF to simplify review, sharing, and archive.
save_plot_pdf(plot_list, plot_pdf_path)

# Print concise summaries at the end so the run produces an immediate overview
# without needing to open the workbook first.
message("\nCombined pipeline complete.")
message("Workbook output: ", basename(workbook_path))
message("Plot PDF output: ", basename(plot_pdf_path))

message("\nMechanical ventilation performance:")
print(mv_performance_summary)

message("\nGentamicin performance:")
print(gent_performance_summary)

message("\nResistance phenotype performance:")
print(phenotype_performance_summary)

message("\nTop mechanical ventilation RF predictors:")
print(utils::head(mv_rf$importance, 10))

message("\nTop gentamicin RF predictors:")
print(utils::head(gent_rf$importance, 10))

message("\nTop phenotype RF predictors:")
print(utils::head(phenotype_rf$importance, 10))
