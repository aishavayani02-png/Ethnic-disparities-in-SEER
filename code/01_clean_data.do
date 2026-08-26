/****************************************************************************************
Create tidy and analysis-ready files from SEER*Stat exports.

Input:  one private .dta file per cancer site, with exported variables v1-v12.
Output: a tidy file retaining exclusion variables, and a clean analysis file.

No patient-level data are written inside the repository.
****************************************************************************************/

version 19.0
clear all
set more off
do "code/00_config.do"

if $PATHS_CONFIGURED != 1 {
    display as error "Edit PRIVATE_DATA_ROOT in code/00_config.do before running."
    exit 198
}

capture mkdir "$TIDY_DATA_DIR"
capture mkdir "$CLEAN_DATA_DIR"

foreach site of global ANALYSIS_SITES {
    display as text "Cleaning `site'..."
    use "$RAW_DATA_DIR/`site'.dta", clear

    rename v1 year
    label variable year "Year of diagnosis"

    rename v2 sex_text
    generate byte sex = 1 if sex_text == "Male"
    replace sex = 2 if sex_text == "Female"
    label define sex_label 1 "Male" 2 "Female", replace
    label values sex sex_label
    label variable sex "Sex"
    drop sex_text

    generate byte race = .
    replace race = 1 if v3 == "White"
    replace race = 2 if v3 == "Black"
    replace race = 3 if v3 == "Other (American Indian/AK Native, Asian/Pacific Islander)"
    replace race = 4 if v3 == "Unknown"
    label define race_label 1 "White" 2 "Black" 3 "Other" 4 "Unknown", replace
    label values race race_label
    label variable race "Race"
    drop v3

    replace v4 = "90 years" if v4 == "90+ years"
    generate str10 age_text = substr(v4, 1, strpos(v4, "years") - 1)
    destring age_text, generate(age)
    label variable age "Age in years"
    drop v4 age_text

    if "`site'" == "Colon" {
        generate byte site = .
        replace site = 0 if v5 == "Appendix"
        replace site = 1 if inlist(v5, "Cecum", "Ascending Colon", "Hepatic Flexure", "Transverse Colon")
        replace site = 2 if inlist(v5, "Splenic Flexure", "Descending Colon", "Sigmoid Colon")
        replace site = 3 if v5 == "Large Intestine, NOS"
        label define site_label 0 "Appendix" 1 "Proximal" 2 "Distal" 3 "Large intestine, NOS", replace
        label values site site_label
        drop v5
    }
    else if "`site'" == "Oral_cavity_and_pharynx" {
        generate byte site = .
        replace site = 1 if inlist(v5, "Nasopharynx", "Oropharynx", "Hypopharynx", "Tonsil")
        replace site = 2 if inlist(v5, "Floor of Mouth", "Gum and Other Mouth", "Salivary Gland")
        replace site = 3 if v5 == "Tongue"
        replace site = 4 if inlist(v5, "Lip", "Other Oral Cavity and Pharynx")
        label define site_label 1 "Pharynx" 2 "Mouth" 3 "Tongue" 4 "Other oral cavity", replace
        label values site site_label
        drop v5
    }
    else {
        encode v5, generate(site)
        label variable site "Location of tumour"
        drop v5
    }

    generate byte stage = .
    replace stage = 0 if v6 == "In situ"
    replace stage = 1 if v6 == "Localized only"
    replace stage = 2 if inlist(v6, ///
        "Regional by both direct extension and lymph node involvement", ///
        "Regional by direct extension only", ///
        "Regional lymph nodes involved only")
    replace stage = 3 if v6 == "Distant site(s)/node(s) involved"
    replace stage = 4 if v6 == "Unknown/unstaged/unspecified/DCO"
    label define stage_label 0 "In situ" 1 "Localised" 2 "Regional" 3 "Distant" 4 "Unknown/unstaged", replace
    label values stage stage_label
    label variable stage "Stage at diagnosis"

    generate byte stage_expanded = .
    replace stage_expanded = 0 if stage == 0
    replace stage_expanded = 1 if stage == 1
    replace stage_expanded = 2 if v6 == "Regional by direct extension only"
    replace stage_expanded = 3 if v6 == "Regional lymph nodes involved only"
    replace stage_expanded = 4 if v6 == "Regional by both direct extension and lymph node involvement"
    replace stage_expanded = 5 if stage == 3
    replace stage_expanded = 6 if stage == 4
    label define expanded_label 0 "In situ" 1 "Localised" 2 "Regional: direct extension" ///
        3 "Regional: lymph node involvement" 4 "Regional: both" 5 "Distant" ///
        6 "Unknown/unstaged", replace
    label values stage_expanded expanded_label
    drop v6

    generate byte cod = 0 if v7 == "Alive or dead of other cause" & v8 == "Alive or dead due to cancer"
    replace cod = 1 if v7 == "Dead (attributable to this cancer dx)"
    replace cod = 2 if v8 == "Dead (attributable to causes other than this cancer dx)"
    replace cod = 3 if v7 == "Dead (missing/unknown COD)"
    label define cod_label 0 "Alive" 1 "Dead due to cancer" 2 "Dead due to other causes" 3 "Unknown", replace
    label values cod cod_label

    generate byte dead = 1 if v11 == "Dead"
    replace dead = 0 if v11 == "Alive"
    label define dead_label 0 "Alive" 1 "Dead", replace
    label values dead dead_label
    drop v7 v8 v11

    replace v9 = "" if v9 == "Unknown"
    destring v9, generate(surv_mm)
    label variable surv_mm "Survival time in months"
    drop v9

    rename v10 SurvMonthsFlag
    rename v12 PrimaryByInternationalRules
    compress
    save "$TIDY_DATA_DIR/`site'_tidy.dta", replace

    if "`site'" == "Breast" keep if sex == 2
    if "`site'" == "Colon" drop if site == 0
    keep if PrimaryByInternationalRules == "Yes"
    drop PrimaryByInternationalRules
    drop if SurvMonthsFlag == "Not calculated because a Death Certificate Only or Autopsy Only case"
    replace surv_mm = surv_mm + 0.5
    replace surv_mm = 0.5 / 30.5 if SurvMonthsFlag == ///
        "Complete dates are available and there are 0 days of survival"
    drop SurvMonthsFlag
    keep if inrange(age, 15, 89)
    keep if inlist(race, 1, 2)
    keep if inlist(stage, 1, 2, 3)
    compress
    save "$CLEAN_DATA_DIR/`site'.dta", replace
}

display as result "Tidy and clean site files created outside the repository."

