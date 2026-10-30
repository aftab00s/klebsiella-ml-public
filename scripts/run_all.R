# Runs everything, in order. Needs the private data file at data/analysis_data.xlsx.
# From the project root: Rscript scripts/run_all.R

source(here::here("scripts", "01_preprocess.R"))
source(here::here("scripts", "02_fit_models.R"))
source(here::here("scripts", "03_evaluate.R"))
