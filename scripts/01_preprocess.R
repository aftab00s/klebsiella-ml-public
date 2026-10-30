# Step 1: read the data and build one cleaned dataset per outcome.

source(here::here("R", "config.R"))
source(here::here("R", "preprocessing.R"))
source(here::here("R", "model_fitting.R"))
source(here::here("R", "evaluation.R"))

# Import the workbook once and print basic structure information to the console
# so the analysis context is clear at run time.
raw_loaded <- read_raw_data(config$input_path, config$preferred_sheet)
raw_df <- raw_loaded$data

message("Workbook: ", basename(config$input_path))
message("Worksheet used: ", raw_loaded$sheet)
message("Rows in raw worksheet: ", nrow(raw_df))
message("Columns in raw worksheet: ", ncol(raw_df))

# Prepare one cleaned analysis dataset per endpoint before any models are fit.
prepared_datasets <- lapply(outcome_specs, function(spec) {
  prepare_outcome_dataset(raw_df, spec, all_outcome_names)
})

saveRDS(
  list(prepared_datasets = prepared_datasets, sheet_used = raw_loaded$sheet),
  file.path(config$output_dir, "prepared_datasets.rds")
)
