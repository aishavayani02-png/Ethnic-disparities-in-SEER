/****************************************************************************************
Project configuration

Run every do-file from the repository root. Edit only PRIVATE_DATA_ROOT below.
The private folder must never be committed to GitHub.

Expected private folder structure:
  raw/    SEER-exported site files named Breast.dta, Colon.dta, ... and popmort.dta
  tidy/   created by code/01_clean_data.do
  clean/  created by code/01_clean_data.do and code/01b_prepare_lifetable.do
****************************************************************************************/

version 19.0

* ----------------------------- USER SETTING -------------------------------------------
global PRIVATE_DATA_ROOT "EDIT_TO_YOUR_PRIVATE_DATA_FOLDER"
* --------------------------------------------------------------------------------------

global REPO_ROOT         "`c(pwd)'"
global RAW_DATA_DIR      "$PRIVATE_DATA_ROOT/raw"
global TIDY_DATA_DIR     "$PRIVATE_DATA_ROOT/tidy"
global CLEAN_DATA_DIR    "$PRIVATE_DATA_ROOT/clean"

global WORK_DIR          "$REPO_ROOT/work"
global MAIN_RESULT_DIR   "$WORK_DIR/main"
global SENS_RESULT_DIR   "$WORK_DIR/sensitivity"
global MODEL_FIT_DIR     "$WORK_DIR/model_fit"
global TABLE_OUTPUT_DIR  "$REPO_ROOT/outputs/tables"
global FIGURE_OUTPUT_DIR "$REPO_ROOT/outputs/figures"
global RELEASE_DIR       "$REPO_ROOT/results"

global ANALYSIS_SITES "Breast Colon Corpus_and_uterus Lung_and_bronchus Oral_cavity_and_pharynx Ovary Pancreas Rectum"
global ANALYSIS_SEED 5693454
global PARAMETRIC_DRAWS 200

global PATHS_CONFIGURED 1
if "$PRIVATE_DATA_ROOT" == "EDIT_TO_YOUR_PRIVATE_DATA_FOLDER" {
    global PATHS_CONFIGURED 0
}

* Create output folders if they are absent. mkdir errors are harmless when folders exist.
capture mkdir "$WORK_DIR"
capture mkdir "$MAIN_RESULT_DIR"
capture mkdir "$SENS_RESULT_DIR"
capture mkdir "$MODEL_FIT_DIR"
capture mkdir "$REPO_ROOT/outputs"
capture mkdir "$TABLE_OUTPUT_DIR"
capture mkdir "$FIGURE_OUTPUT_DIR"
