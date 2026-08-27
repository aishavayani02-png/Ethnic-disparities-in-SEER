# Stage at diagnosis and Black-White cancer survival disparities

This repository contains the code and aggregate results for:

> **Quantifying the role of stage-at-diagnosis in racial disparities in US cancer
> survival: a population-based SEER study for 8 common cancers**

Authors: Aisha Vayani, Mark Rutherford, and Sarah Booth, University of Leicester.

## Files

- `code/` contains the Stata code for cleaning, analysis, sensitivity checks,
  model-fit assessment, tables, and figures.
- `results/` contains the aggregate estimates and counts used in the paper.
- `shiny/` contains the Shiny application and its aggregate input file.
- `docs/variable_dictionary.csv` explains the variables used in the analysis.
- `data/README.md` lists the SEER files needed to rerun the work.

No patient-level SEER data are included. Researchers who wish to rerun the analysis
from the original records must obtain their own access through the SEER Data Request
System and comply with the SEER Research Data Agreement.

## Running the analysis

The analysis was conducted in Stata 19. The Shiny application was developed in R
4.5.1. Run the Stata files from the repository root.

1. Create the private data folders described in `data/README.md`.
2. Set `PRIVATE_DATA_ROOT` in `code/00_config.do`.
3. Run `code/01_clean_data.do` and `code/01b_prepare_population_lifetable.do`.
4. Run `code/02_table_1.do`, `code/02b_black_counts.do`,
   `code/03_main_analysis.do`, `code/04_model_fit.do`, and
   `code/05_sensitivity_analysis.do`.
5. Run `code/06_combine_results.do` to create the combined result files.
6. Run `code/run_postprocessing.do` to create the manuscript tables and figures.

`code/run_all.do` runs the same steps in order. The model-fitting stages take the
longest to complete.

## Results

The CSV files contain aggregate estimates or counts, not individual patient records.
The main results include the total indirect effect (TIE), total direct effect (TDE),
total causal effect (TCE), and proportion mediated. The table source files contain
the values used in Tables 1-6.

## Shiny application

The application code can be run locally from the `shiny/` folder. The existing
version is also available at
[g9v1qs-aisha-vayani.shinyapps.io/Cancer](https://g9v1qs-aisha-vayani.shinyapps.io/Cancer/).

## Licence

The code is released under the MIT License. The licence does not apply to SEER data.
