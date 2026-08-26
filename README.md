# Stage at diagnosis and Black-White cancer survival disparities

This repository contains the analysis code and disclosure-safe aggregate results for:

> **Quantifying the role of stage-at-diagnosis in racial disparities in US cancer
> survival: a population-based SEER study for 8 common cancers**

Authors: Aisha Vayani, Mark Rutherford, and Sarah Booth, University of Leicester.

## What is included

- Stata code for data cleaning, cohort summaries, causal mediation analyses, model-fit
  assessment, sensitivity analyses, result assembly, Tables 1-6, and Figures 1-2.
- Long-format aggregate estimates used in the manuscript and Shiny application.
- The R/Shiny source code and its aggregate input file.
- A variable dictionary and a description of the reproducibility workflow.

No patient-level SEER data are included. SEER data access and redistribution are
governed by the SEER Research Data Agreement. Researchers wishing to rerun the
analysis from source records must obtain their own SEER data access and create the
specified SEER*Stat exports.

## Repository map

| Location | Contents |
|---|---|
| `code/` | Stata cleaning, modelling, post-processing, table, and figure scripts |
| `results/main/` | Main-analysis aggregate estimates, including TIE, TDE, and TCE |
| `results/sensitivity/` | Aggregate sensitivity-analysis estimates |
| `results/model_fit/` | Flexible parametric model and Pohar Perme estimates |
| `results/manuscript_tables/` | Numeric table sources and manuscript table grids |
| `shiny/` | Shiny application and disclosure-safe aggregate input |
| `docs/` | Workflow, software, variable, and disclosure-control documentation |

## Reproducing the analysis

The analysis was conducted in Stata 19. The Shiny application was developed in R
4.5.1. Run Stata from the repository root so that relative paths resolve correctly.

1. Obtain authorised access to the required SEER data and create the eight site files
   described in `data/README.md`.
2. Edit `PRIVATE_DATA_ROOT` in `code/00_config.do`. This is the only path that users
   should need to change.
3. Run `code/01_clean_data.do` and `code/01b_prepare_population_lifetable.do`.
4. Run `code/02_table_1.do` and `code/tables/table_4.do` for the descriptive and
   exclusion tables.
5. Run `code/03_main_analysis.do`, `code/04_model_fit.do`, and
   `code/05_sensitivity_analysis.do`. These models are computationally intensive.
6. Run `code/02b_black_counts.do` and `code/06_combine_results.do`.
7. Run `code/run_postprocessing.do` to export Tables 2, 3, 5, and 6 and Figures 1-2.

`code/run_all.do` provides the same sequence in one file. The seed is reset to
5693454 for each cancer site, reflecting the site-by-site execution of the submitted
analysis. Confidence intervals use 200 parametric draws of the mediator-model
coefficients, conditional on the fitted survival model.

## Using the released results

The CSV files in `results/` are aggregate model estimates or aggregate counts. They
can be inspected without SEER access and are sufficient to recreate the manuscript
figures, the model-result tables, and the Shiny application.

The direct effect (`tde`) is deliberately retained in the main result files and in
Table 5, matching the current manuscript structure.

## Shiny application

GitHub stores the application source but does not run a Shiny server. The existing
deployed application is available at
[g9v1qs-aisha-vayani.shinyapps.io/Cancer](https://g9v1qs-aisha-vayani.shinyapps.io/Cancer/).
See `shiny/README.md` for local use and redeployment.

## Citation

Please cite the associated article and this repository. Machine-readable citation
metadata are provided in `CITATION.cff`.

## Licence

Code is released under the MIT License. The licence does not apply to SEER data, which
remain subject to the SEER Research Data Agreement.

