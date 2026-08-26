/****************************************************************************************
Combine private working result files into disclosure-safe long-format result datasets.

Prerequisites:
  code/03_main_analysis.do
  code/04_model_fit.do
  code/05_sensitivity_analysis.do
  code/02b_black_counts.do

The released CSVs contain model estimates and aggregate counts only. They contain no
patient-level SEER records.
****************************************************************************************/

version 19.0
clear all
set more off
do "code/00_config.do"

capture program drop append_effect_file
program define append_effect_file
    syntax, File(string) Site(string) Sex(string) Estimate(string) ///
        Agegroup(integer) Master(string)

    use "`file'", clear
    capture confirm variable PE
    if _rc {
        quietly ds *_PE
        local candidate : word 1 of `r(varlist)'
        rename `candidate' PE
    }
    capture confirm variable lci
    if _rc {
        quietly ds *_lci
        local candidate : word 1 of `r(varlist)'
        rename `candidate' lci
    }
    capture confirm variable uci
    if _rc {
        quietly ds *_uci
        local candidate : word 1 of `r(varlist)'
        rename `candidate' uci
    }

    confirm variable tt PE lci uci
    keep tt PE lci uci
    generate str40 site = "`site'"
    generate str10 sex = "`sex'"
    generate str12 estimate = "`estimate'"
    generate byte agegroup = `agegroup'
    generate str10 race = ""
    generate str12 stage = ""
    order site sex estimate race stage agegroup tt PE lci uci
    append using "`master'"
    save "`master'", replace
end

capture program drop append_rs_file
program define append_rs_file
    syntax, File(string) Site(string) Sex(string) Agegroup(integer) Master(string)

    use "`file'", clear
    foreach group in lw rw dw lb rb db {
        rename `group' PE_`group'
        rename `group'_lci lci_`group'
        rename `group'_uci uci_`group'
    }
    reshape long PE_ lci_ uci_, i(tt) j(group) string
    rename PE_ PE
    rename lci_ lci
    rename uci_ uci
    generate str40 site = "`site'"
    generate str10 sex = "`sex'"
    generate str12 estimate = "acs"
    generate byte agegroup = `agegroup'
    generate str10 race = cond(substr(group, 2, 1) == "w", "White", "Black")
    generate str12 stage = cond(substr(group, 1, 1) == "l", "Localised", ///
        cond(substr(group, 1, 1) == "r", "Regional", "Distant"))
    keep site sex estimate race stage agegroup tt PE lci uci
    append using "`master'"
    save "`master'", replace
end

* Main causal-effect result files --------------------------------------------------------
tempfile main_long
clear
save `main_long', emptyok replace

foreach site of global ANALYSIS_SITES {
    foreach sex in Male Female {
        foreach estimate in AD ADa ADb tce tde tie pm {
            local all_file "$MAIN_RESULT_DIR/`estimate'_`site'_`sex'.dta"
            capture confirm file "`all_file'"
            if !_rc {
                quietly append_effect_file, file("`all_file'") site("`site'") ///
                    sex("`sex'") estimate("`estimate'") agegroup(0) master("`main_long'")
            }
            forvalues agegroup = 1/2 {
                local age_file "$MAIN_RESULT_DIR/`estimate'`agegroup'_`site'_`sex'.dta"
                capture confirm file "`age_file'"
                if !_rc {
                    quietly append_effect_file, file("`age_file'") site("`site'") ///
                        sex("`sex'") estimate("`estimate'") agegroup(`agegroup') master("`main_long'")
                }
            }
        }
    }
}

use `main_long', clear
sort site sex estimate agegroup tt
save "$MAIN_RESULT_DIR/results_long.dta", replace
export delimited using "$RELEASE_DIR/main/results_long.csv", replace

foreach site of global ANALYSIS_SITES {
    foreach sex in Male Female {
        preserve
        keep if site == "`site'" & sex == "`sex'"
        quietly count
        if r(N) > 0 export delimited using ///
            "$RELEASE_DIR/main/by_site/`site'_`sex'.csv", replace
        restore
    }
}

* Sensitivity-analysis result files ----------------------------------------------------
tempfile sensitivity_long
clear
save `sensitivity_long', emptyok replace

foreach site of global ANALYSIS_SITES {
    foreach sex in Male Female {
        foreach estimate in ADb tce tie pm {
            local all_file "$SENS_RESULT_DIR/`estimate'_`site'_`sex'.dta"
            capture confirm file "`all_file'"
            if !_rc {
                quietly append_effect_file, file("`all_file'") site("`site'") ///
                    sex("`sex'") estimate("`estimate'") agegroup(0) master("`sensitivity_long'")
            }
            forvalues agegroup = 1/2 {
                local age_file "$SENS_RESULT_DIR/`estimate'`agegroup'_`site'_`sex'.dta"
                capture confirm file "`age_file'"
                if !_rc {
                    quietly append_effect_file, file("`age_file'") site("`site'") ///
                        sex("`sex'") estimate("`estimate'") agegroup(`agegroup') master("`sensitivity_long'")
                }
            }
        }
    }
}

use `sensitivity_long', clear
sort site sex estimate agegroup tt
save "$SENS_RESULT_DIR/results_long.dta", replace
export delimited using "$RELEASE_DIR/sensitivity/results_long.csv", replace

foreach site of global ANALYSIS_SITES {
    foreach sex in Male Female {
        preserve
        keep if site == "`site'" & sex == "`sex'"
        quietly count
        if r(N) > 0 export delimited using ///
            "$RELEASE_DIR/sensitivity/by_site/`site'_`sex'.csv", replace
        restore
    }
}

* Standardised survival by race and stage ---------------------------------------------
tempfile survival_long
clear
save `survival_long', emptyok replace

foreach site of global ANALYSIS_SITES {
    foreach sex in Male Female {
        forvalues agegroup = 1/2 {
            local rs_file "$MAIN_RESULT_DIR/rs`agegroup'_`site'_`sex'.dta"
            capture confirm file "`rs_file'"
            if !_rc {
                quietly append_rs_file, file("`rs_file'") site("`site'") ///
                    sex("`sex'") agegroup(`agegroup') master("`survival_long'")
            }
        }
    }
}
use `survival_long', clear
sort site sex race stage agegroup tt
save "$MAIN_RESULT_DIR/standardised_survival_long.dta", replace
export delimited using "$RELEASE_DIR/main/standardised_survival_long.csv", replace

* Model-fit output ---------------------------------------------------------------------
use "$MODEL_FIT_DIR/model_fit.dta", clear
sort site sex race stage agegroup estimate tt
export delimited using "$RELEASE_DIR/model_fit/model_fit.csv", replace

* Combined Shiny input -----------------------------------------------------------------
use "$MAIN_RESULT_DIR/results_long.dta", clear
append using "$MAIN_RESULT_DIR/standardised_survival_long.dta"
append using "$MODEL_FIT_DIR/model_fit.dta"
keep if inlist(agegroup, 1, 2)

merge m:1 site sex agegroup using "$WORK_DIR/black_counts.dta", keep(master match) nogen
replace n_black = . if !inlist(estimate, "AD", "ADa", "ADb")

replace site = "Colon (excluding Appendix)" if site == "Colon"
replace site = "Uterine corpus" if site == "Corpus_and_uterus"
replace site = "Lung and bronchus" if site == "Lung_and_bronchus"
replace site = "Oral cavity and pharynx" if site == "Oral_cavity_and_pharynx"
generate str4 agegroup_label = cond(agegroup == 1, "<75", "75+")
drop agegroup
rename agegroup_label agegroup
order site estimate agegroup stage sex race tt PE lci uci n_black
sort site estimate agegroup stage sex race tt

save "$WORK_DIR/Shiny.dta", replace
export delimited using "$RELEASE_DIR/combined/Shiny.csv", replace
copy "$RELEASE_DIR/combined/Shiny.csv" "$REPO_ROOT/shiny/data/Shiny.csv", replace

display as result "Main, sensitivity, model-fit, and Shiny release files created."

