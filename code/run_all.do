version 19.0
clear all
set more off

* Full rerun. The modelling stages are computationally intensive.
do "code/01_clean_data.do"
do "code/01b_prepare_population_lifetable.do"
do "code/02_table_1.do"
do "code/tables/table_4.do"
do "code/03_main_analysis.do"
do "code/04_model_fit.do"
do "code/05_sensitivity_analysis.do"
do "code/02b_black_counts.do"
do "code/06_combine_results.do"
do "code/run_postprocessing.do"

