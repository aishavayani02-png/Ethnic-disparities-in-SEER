/****************************************************************************************
Prepare the SEER population mortality file used for expected survival.

Input:  PRIVATE_DATA_ROOT/raw/popmort.dta with SEER*Stat fields v1-v5.
Output: PRIVATE_DATA_ROOT/clean/popmort.dta.
****************************************************************************************/

version 19.0
clear all
set more off
do "code/00_config.do"

if $PATHS_CONFIGURED != 1 {
    display as error "Edit PRIVATE_DATA_ROOT in code/00_config.do before running."
    exit 198
}

use "$RAW_DATA_DIR/popmort.dta", clear

generate double prob = v5 / 1000000
generate double rate = -ln(prob)
drop v5

rename v4 sex_text
generate byte sex = 1 if sex_text == "Male"
replace sex = 2 if sex_text == "Female"
label define sex_label 1 "Male" 2 "Female", replace
label values sex sex_label
drop sex_text

generate byte race = .
replace race = 1 if v3 == "White"
replace race = 2 if v3 == "Black"
replace race = 3 if v3 == "Other (American Indian/AK Native, Asian/Pacific Islander)"
replace race = 7 if v3 == "Other unspecified (1991+)"
replace race = 9 if v3 == "Unknown"
keep if inlist(race, 1, 2)
label define race_label 1 "White" 2 "Black", replace
label values race race_label
drop v3

rename v1 age
rename v2 year
order year age sex race prob rate
isid year age sex race
compress
save "$CLEAN_DATA_DIR/popmort.dta", replace

display as result "Population mortality file created outside the repository."

