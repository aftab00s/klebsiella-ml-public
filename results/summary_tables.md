# Summary tables

Internal cross-validation results from the dissertation, on a small sample (n = 30) with no separate test set. Intervals are 95% bootstrap intervals that do not account for tuning, so they are likely too narrow.

## Model performance

| Outcome | Model | Balanced accuracy | 95% CI (bootstrap) | Accuracy | Sensitivity | Specificity | Kappa | AUC | p vs majority class | n used | n scored |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Mechanical ventilation | GLM (stepwise) | 0.458 | 0.29-0.64 | 0.481 | 0.250 | 0.667 | -0.086 | 0.506 | 0.834 | 27 | 27 |
| Mechanical ventilation | Random forest | 0.558 | 0.37-0.74 | 0.556 | 0.583 | 0.533 | 0.115 | 0.472 | 0.579 | 27 | 27 |
| Gentamicin resistance | GLM (stepwise) | 0.723 | 0.56-0.87 | 0.733 | 0.571 | 0.875 | 0.455 | 0.864 | 0.020 | 30 | 30 |
| Gentamicin resistance | Random forest | 0.808 | 0.67-0.94 | 0.800 | 0.929 | 0.688 | 0.605 | 0.875 | 0.002 | 30 | 30 |
| Resistance phenotype | Multinomial model | 0.750 | 0.57-0.89 | 0.750 | 0.750 | 0.750 | 0.500 | n/a | 0.006 | 29 | 28 |
| Resistance phenotype | Random forest | 0.821 | 0.68-0.96 | 0.821 | 0.821 | 0.821 | 0.643 | n/a | <0.001 | 29 | 28 |

## Predictor groups

| Predictor group | Columns |
|---|---|
| Typing (sequence, capsule and O-antigen type) | 3 |
| Resistance genes | 92 |
| Virulence genes | 60 |
| Disinfectant and heavy-metal genes | 33 |
| Mutations | 39 |
| Plasmid markers | 21 |
| Antibiotic test results | 14 |
| Biofilm phenotype | 1 |
| Total predictors | 263 |

## Samples per outcome

| Outcome | Type | Samples used | Samples scored | Note |
|---|---|---|---|---|
| Gentamicin resistance | Binary | 30 | 30 | All isolates |
| Resistance phenotype | Multiclass (3 classes) | 29 | 28 | 1 invalid code excluded; 1 sample in a rare class not scored |
| Mechanical ventilation | Binary | 27 | 27 | 3 outcome values missing |

## The 91 configurations

| Outcome | Regression configurations | Regression type | mtry values | ntree values | nodesize values | Forest configurations | Total |
|---|---|---|---|---|---|---|---|
| Mechanical ventilation | 5 | Null, forward, backward, stepwise, full | 5 | 2 | 3 | 30 | 35 |
| Gentamicin resistance | 5 | Null, forward, backward, stepwise, full | 5 | 2 | 3 | 30 | 35 |
| Resistance phenotype | 1 | Multinomial | 5 | 2 | 2 | 20 | 21 |
| Total | 11 |  |  |  |  | 80 | 91 |
