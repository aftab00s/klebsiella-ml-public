# Step 2: build the folds and fit every model.

source(here::here("R", "config.R"))
source(here::here("R", "preprocessing.R"))
source(here::here("R", "model_fitting.R"))
source(here::here("R", "evaluation.R"))

prepared <- readRDS(file.path(config$output_dir, "prepared_datasets.rds"))
prepared_datasets <- prepared$prepared_datasets

# Set once, before the folds. The fits below must stay in this order: moving or
# removing one changes the random numbers every later step receives.
set.seed(config$seed)

# Create a validation plan for each outcome separately because class balance and
# usable sample size differ across ventilation, gentamicin, and phenotype.
fold_plans <- lapply(names(outcome_specs), function(key) {
  spec <- outcome_specs[[key]]
  dataset_info <- prepared_datasets[[key]]
  make_fold_plan(dataset_info$data[[spec$cleaned_outcome]], config$requested_k, spec$type)
})
names(fold_plans) <- names(outcome_specs)

# Run ventilation as the first binary endpoint because it is the original
# working workflow that the combined script generalises.
mv_glm <- run_binary_glm_pipeline(prepared_datasets$mv$data, outcome_specs$mv, fold_plans$mv)
mv_rf <- run_binary_rf_pipeline(prepared_datasets$mv$data, outcome_specs$mv, fold_plans$mv)

# Compare five GLM variants for the binary ventilation endpoint so the
# dissertation can report how null, selected, and screened-full models differ.
mv_glm_model_comparison <- run_binary_glm_model_comparison(
  prepared_datasets$mv$data,
  outcome_specs$mv,
  fold_plans$mv
)

# Apply the same binary framework to gentamicin so both binary analyses are
# directly comparable in preprocessing, validation, and reporting.
gent_glm <- run_binary_glm_pipeline(prepared_datasets$gent$data, outcome_specs$gent, fold_plans$gent)
gent_rf <- run_binary_rf_pipeline(prepared_datasets$gent$data, outcome_specs$gent, fold_plans$gent)

# Repeat the five-GLM comparison for gentamicin resistance using the same
# binary validation strategy as mechanical ventilation.
gent_glm_model_comparison <- run_binary_glm_model_comparison(
  prepared_datasets$gent$data,
  outcome_specs$gent,
  fold_plans$gent
)

# Use multiclass modelling for resistance phenotype rather than forcing binary
# logic onto an endpoint with more than two categories.
phenotype_model <- run_multiclass_model_pipeline(prepared_datasets$phenotype$data, outcome_specs$phenotype, fold_plans$phenotype)
phenotype_rf <- run_multiclass_rf_pipeline(prepared_datasets$phenotype$data, outcome_specs$phenotype, fold_plans$phenotype)

saveRDS(
  list(
    fold_plans = fold_plans,
    mv_glm = mv_glm, mv_rf = mv_rf, mv_glm_model_comparison = mv_glm_model_comparison,
    gent_glm = gent_glm, gent_rf = gent_rf, gent_glm_model_comparison = gent_glm_model_comparison,
    phenotype_model = phenotype_model, phenotype_rf = phenotype_rf
  ),
  file.path(config$output_dir, "fitted_models.rds")
)
