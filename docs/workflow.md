# Analysis workflow

## 1. Data preparation

`01_clean_data.do` labels the SEER*Stat export fields, derives sex, race, age, cancer
subsite, stage, vital status, cause of death, and survival time, and applies the
primary-analysis restrictions. It saves a tidy file before exclusions and a clean
analysis file after exclusions. `01b_prepare_population_lifetable.do` prepares the
expected-mortality rates used for relative survival.

## 2. Descriptive outputs

`02_table_1.do` constructs the 2021-2022 period window and calculates aggregate stage
counts, row percentages, and age summaries. `tables/table_4.do` assigns each case to
the first applicable exclusion and calculates sequential exclusion percentages.

## 3. Main analysis

`03_main_analysis.do` fits separate flexible parametric excess-hazard models by cancer
site and sex. It models stage distributions separately for White and Black patients,
standardises over Black patients, estimates total, direct, and indirect effects, and
translates risk differences into deaths postponed for aggregate 2022 cohort sizes.

The script saves one private working result file for each effect, site, sex, and age
group. It does not create figures or manuscript tables.

## 4. Model fit and sensitivity analysis

`04_model_fit.do` compares fitted standardised survival with the non-parametric Pohar
Perme estimator. `05_sensitivity_analysis.do` repeats the mediation analysis after
assigning unknown/unstaged cases to distant stage.

## 5. Result assembly and reporting

`02b_black_counts.do` creates aggregate 2022 Black incident-cohort counts.
`06_combine_results.do` converts the individual working files to consistent long
format, exports per-site CSVs, and assembles `results/combined/Shiny.csv`.

The scripts under `code/tables/` export Tables 2-6 directly to Word. The scripts under
`code/figures/` recreate Figures 1 and 2 from the combined aggregate dataset.

