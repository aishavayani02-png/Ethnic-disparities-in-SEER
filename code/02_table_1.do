/****************************************************************************************
Table 1: stage distribution and age at diagnosis among period-analysis-eligible cases.

The random seed is reset by cancer site to reproduce the site-by-site analyses.
This script writes an aggregate CSV and a formatted Word table; no individual records
are exported.
****************************************************************************************/

version 19.0
clear all
set more off
do "code/00_config.do"

if $PATHS_CONFIGURED != 1 {
    display as error "Edit PRIVATE_DATA_ROOT in code/00_config.do before running."
    exit 198
}

tempfile master
clear
set obs 0
generate str40 cancer_site = ""
generate byte sex = .
generate byte age_group = .
generate byte race = .
generate byte stage = .
generate long count = .
generate double mean_age = .
generate double sd_age = .
generate long denominator = .
generate double percentage = .
save `master', replace

foreach site of global ANALYSIS_SITES {
    set seed $ANALYSIS_SEED
    use "$CLEAN_DATA_DIR/`site'.dta", clear
    levelsof sex if !missing(sex), local(sexes)

    foreach sex_value of local sexes {
        use "$CLEAN_DATA_DIR/`site'.dta", clear
        keep if sex == `sex_value'
        generate long id = _n

        generate byte month_diag = runiformint(1, 12)
        generate byte day_diag = runiformint(1, 28)
        generate double random_diagnosis = mdy(month_diag, day_diag, year)
        generate double date_exit = dofm(mofd(random_diagnosis) + surv_mm)
        generate double censored_diagnosis = dofm(mofd(date_exit) - surv_mm)
        generate double date_diagnosis = cond(dead == 1, random_diagnosis, censored_diagnosis)

        stset date_exit, enter(time mdy(1, 1, 2021)) ///
            exit(time min(date_diagnosis + 5.1 * 365.24, mdy(12, 31, 2022))) ///
            origin(date_diagnosis) failure(dead == 1) scale(365.24) id(id)

        rename age age_diagnosis
        rename year year_diagnosis
        generate byte age = min(int(age_diagnosis + _t), 99)
        generate int year = min(floor(year_diagnosis + _t), 2021)
        sort sex year age race
        merge m:1 sex year age race using "$CLEAN_DATA_DIR/popmort.dta", ///
            keep(match master) keepusing(rate)
        drop _merge

        keep if _st == 1 & _t > _t0
        generate byte age_group = cond(age_diagnosis < 75, 1, 2)
        bysort sex age_group race: egen double stratum_mean = mean(age_diagnosis)
        bysort sex age_group race: egen double stratum_sd = sd(age_diagnosis)

        collapse (count) count=age_diagnosis ///
            (firstnm) mean_age=stratum_mean sd_age=stratum_sd, ///
            by(sex age_group race stage)
        bysort sex age_group race: egen long denominator = total(count)
        generate double percentage = 100 * count / denominator
        generate str40 cancer_site = "`site'"

        append using `master'
        save `master', replace
    }
}

use `master', clear
drop if missing(count)
sort cancer_site sex age_group race stage
export delimited using "$RELEASE_DIR/manuscript_tables/table_1_source.csv", replace

quietly summarize count
assert r(sum) == 1299386

generate str32 n_percent = strtrim(string(count, "%12.0fc")) + " (" + ///
    strtrim(string(percentage, "%5.2f")) + ")"
generate str24 mean_sd = strtrim(string(mean_age, "%5.2f")) + " (" + ///
    strtrim(string(sd_age, "%5.2f")) + ")"
keep cancer_site sex age_group race stage n_percent mean_sd
reshape wide n_percent, i(cancer_site sex age_group race mean_sd) j(stage)
rename n_percent1 localised
rename n_percent2 regional
rename n_percent3 distant
decode sex, generate(sex_label)
decode race, generate(race_label)
generate str4 age_label = cond(age_group == 1, "<75", "75+")
replace cancer_site = "Colon (excluding Appendix)" if cancer_site == "Colon"
replace cancer_site = "Uterine corpus" if cancer_site == "Corpus_and_uterus"
replace cancer_site = "Lung and bronchus" if cancer_site == "Lung_and_bronchus"
replace cancer_site = "Oral cavity and pharynx" if cancer_site == "Oral_cavity_and_pharynx"
sort cancer_site sex age_group race

putdocx clear
putdocx begin, pagesize(A4) landscape margin(top, 0.5) margin(bottom, 0.5) ///
    margin(left, 0.5) margin(right, 0.5)
putdocx paragraph, style(Title)
putdocx text ("Table 1. Stage at diagnosis and age among period-analysis-eligible patients")
putdocx table table1 = data(cancer_site sex_label age_label race_label ///
    localised regional distant mean_sd), varnames layout(autofitcontents)
putdocx table table1(1,.), bold
putdocx table table1(.,.), font("Calibri", 8)
putdocx save "$TABLE_OUTPUT_DIR/Table_1.docx", replace

display as result "Table 1 source CSV and Word table created."

