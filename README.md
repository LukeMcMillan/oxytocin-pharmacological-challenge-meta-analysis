# oxytocin-pharmacological-challenge-meta-analysis
# Peripheral Oxytocin Response to Pharmacological Challenge: 
# A Multilevel Meta-Analysis

Pre-registered systematic review and multilevel meta-analysis 
(PROSPERO: CRD42024625856)

## Software
R version 4.5.1
- metafor 4.8.0
- clubSandwich 0.6.2

## Reproducibility
All analysis scripts are provided. The extraction dataset 
(OT_meta_analysis_clean_v5.xlsx) is available on request 
from the corresponding author.

# Peripheral Oxytocin Response to Pharmacological Challenge: A Multilevel Meta-Analysis

This repository contains the analysis code for the project:

**Peripheral Oxytocin Response to Pharmacological Challenge: A Multilevel Meta-Analysis**

The analysis estimates peripheral oxytocin responses to pharmacological challenge using multilevel meta-analytic models. The script includes data cleaning, placebo/control exclusion, effect-size calculation, primary multilevel modeling, sensitivity analyses, drug-class moderator analyses, study-collapsed bias checks, sex moderator analyses, and methodological moderator analyses for extraction status and assay type.

## Repository structure

```text
README.md
analysis/
  oxytocin_meta_analysis.R
data/
  OT_meta_analysis_clean_v5.xlsx
results/
figures/
Main analysis script

The main R script is:

analysis/oxytocin_meta_analysis.R
Data file

The analysis expects the Excel dataset to be located at:

data/OT_meta_analysis_clean_v5.xlsx

In the R script, the file path should be:

file_path <- "data/OT_meta_analysis_clean_v5.xlsx"

If the Excel file is renamed, the file_path line in the R script must be updated to match the new filename exactly.

Expected data columns

The Excel file should include the following columns:

study
group
n
mean_pre
sd_pre
mean_post
sd_post
sample_id
drug_class
percent_female
extraction_status
assay_type
Required R packages

Install the required R packages before running the analysis:

install.packages(c(
  "readxl",
  "janitor",
  "dplyr",
  "stringr",
  "metafor",
  "clubSandwich",
  "ggplot2"
))
How to run the analysis

From the main repository folder, run:

source("analysis/oxytocin_meta_analysis.R")

The script reads the Excel file from the data/ folder and writes output files to the results/ folder.

Outputs

The script generates the following files:

results/dat_active_clean.csv
results/class_counts.csv
results/drug_class_results.csv
results/sensitivity_results.csv
results/study_collapsed_es.csv
results/pairwise_contrasts_BH.csv
results/extraction_means.csv
results/assay_means.csv
results/OT_metaanalysis_clean_workspace.RData
