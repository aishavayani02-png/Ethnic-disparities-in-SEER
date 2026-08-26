/****************************************************************************************
Create aggregate 2022 Black incident-cohort counts used to express indirect effects
as numbers of deaths postponed. Unknown/unstaged cases are retained by assigning them
to distant stage, as described in the manuscript.
****************************************************************************************/

version 19.0
clear all
set more off
do "code/00_config.do"

if $PATHS_CONFIGURED != 1 {
    display as error "Edit PRIVATE_DATA_ROOT in code/00_config.do before running."
    exit 198
}

tempname counts_post
tempfile counts
postfile `counts_post' str40 site str10 sex byte agegroup long n_black using `counts', replace

foreach site of global ANALYSIS_SITES {
    use "$TIDY_DATA_DIR/`site'_tidy.dta", clear
    if "`site'" == "Breast" keep if sex == 2
    if "`site'" == "Colon" drop if site == 0
    keep if PrimaryByInternationalRules == "Yes"
    drop if SurvMonthsFlag == "Not calculated because a Death Certificate Only or Autopsy Only case"
    keep if inrange(age, 15, 89)
    keep if inlist(race, 1, 2)
    drop if stage == 0 | missing(stage)
    replace stage = 3 if stage == 4
    generate byte agegroup = cond(age < 75, 1, 2)

    levelsof sex, local(sexes)
    foreach sex_value of local sexes {
        local sex_label : label (sex) `sex_value'
        quietly count if sex == `sex_value' & race == 2 & year == 2022
        post `counts_post' ("`site'") ("`sex_label'") (0) (r(N))
        forvalues age_value = 1/2 {
            quietly count if sex == `sex_value' & race == 2 & year == 2022 & agegroup == `age_value'
            post `counts_post' ("`site'") ("`sex_label'") (`age_value') (r(N))
        }
    }
}

postclose `counts_post'
use `counts', clear
sort site sex agegroup
save "$WORK_DIR/black_counts.dta", replace
export delimited using "$RELEASE_DIR/combined/black_counts.csv", replace

display as result "Aggregate Black incident-cohort counts created."

