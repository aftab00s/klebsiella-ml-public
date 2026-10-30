# Packages, settings and outcome definitions used by every script.

library(caret)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(MASS)
library(nnet)
library(openxlsx)
library(pROC)
library(randomForest)
library(readxl)
library(tibble)
library(tidyr)

options(stringsAsFactors = FALSE, scipen = 999)
theme_set(ggplot2::theme_minimal(base_size = 12))

# All paths are relative to the project root. The clinical data is private, so
# data/ is git-ignored and must be supplied locally.
config <- list(
  seed = 20260402L,
  input_path = here::here("data", "analysis_data.xlsx"),
  preferred_sheet = "original_data",
  output_dir = here::here("results", "output"),
  workbook_name = "results.xlsx",
  plot_pdf_name = "plots.pdf",
  sample_id = "sample_id",
  requested_k = 5L,
  max_binary_df = 6L,
  max_multiclass_df = 3L
)

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

# Define each outcome once so the same preprocessing and modelling framework
# can be reused for binary and multiclass endpoints without duplicating logic.
outcome_specs <- list(
  mv = list(
    key = "mv",
    display_name = "Mechanical ventilation",
    raw_column_candidates = c("mechanical_ventilation", "mechanical ventilation"),
    cleaned_outcome = "mechanical_ventilation",
    type = "binary",
    positive_class = "Vent",
    negative_class = "NoVent",
    label_map = c("0" = "NoVent", "1" = "Vent")
  ),
  gent = list(
    key = "gent",
    display_name = "Gentamicin resistance",
    raw_column_candidates = c("gentamicin", "Gentamicin"),
    cleaned_outcome = "gentamicin",
    type = "binary",
    positive_class = "Resistant",
    negative_class = "NotResistant",
    label_map = c("0" = "NotResistant", "1" = "Resistant")
  ),
  phenotype = list(
    key = "phenotype",
    display_name = "Resistance phenotype",
    raw_column_candidates = c("resistance_phenotype", "resistance phenotype"),
    cleaned_outcome = "resistance_phenotype",
    type = "multiclass",
    positive_class = NA_character_,
    negative_class = NA_character_,
    label_map = c("0" = "MDR", "1" = "XDR", "2" = "PDR")
  )
)

all_outcome_names <- sapply(outcome_specs, function(spec) spec$cleaned_outcome)
