/****************************************************************************************
Table 4: sequential exclusions and final analysis sample.

Each case is assigned to its first applicable exclusion. Exclusion percentages use
the number remaining immediately before that step; the included percentage uses all
eligible White/Black cases before sequential exclusions.
****************************************************************************************/

version 19.0
clear all
set more off
do "code/00_config.do"

if $PATHS_CONFIGURED != 1 {
    display as error "Edit PRIVATE_DATA_ROOT in code/00_config.do before running."
    exit 198
}

tempfile all_results
clear
save `all_results', emptyok replace

foreach site of global ANALYSIS_SITES {
    use "$TIDY_DATA_DIR/`site'_tidy.dta", clear
    if "`site'" == "Breast" keep if sex == 2
    if "`site'" == "Colon" drop if site == 0
    keep if inlist(race, 1, 2)

    generate str40 cancer_site = "`site'"
    generate byte age_group = 1 if age < 75
    replace age_group = 2 if age >= 75 & age < .
    drop if missing(age_group)

    generate byte exclusion = .
    replace exclusion = 1 if PrimaryByInternationalRules != "Yes" & missing(exclusion)
    replace exclusion = 2 if SurvMonthsFlag == ///
        "Not calculated because a Death Certificate Only or Autopsy Only case" & missing(exclusion)
    replace exclusion = 3 if age < 15 & missing(exclusion)
    replace exclusion = 4 if age > 89 & missing(exclusion)
    replace exclusion = 5 if stage == 0 & missing(exclusion)
    replace exclusion = 6 if (stage == 4 | missing(stage)) & missing(exclusion)

    generate long cases = 1
    forvalues reason = 1/6 {
        generate long excluded_`reason' = exclusion == `reason'
    }
    generate long included = missing(exclusion)

    collapse (sum) cases excluded_1-excluded_6 included, ///
        by(cancer_site sex age_group race)

    generate long denominator_1 = cases
    forvalues reason = 2/6 {
        local previous = `reason' - 1
        generate long denominator_`reason' = denominator_`previous' - excluded_`previous'
    }
    forvalues reason = 1/6 {
        generate double percent_`reason' = 100 * excluded_`reason' / denominator_`reason'
    }
    generate double percent_included = 100 * included / cases

    append using `all_results'
    save `all_results', replace
}

use `all_results', clear
sort cancer_site sex age_group race
export delimited using "$RELEASE_DIR/manuscript_tables/table_4_source.csv", replace

generate str24 cases_text = strtrim(string(cases, "%12.0fc")) + " (100)"
forvalues reason = 1/6 {
    generate str28 excluded_text_`reason' = strtrim(string(excluded_`reason', "%12.0fc")) + ///
        " (" + strtrim(string(percent_`reason', "%5.2f")) + ")"
}
generate str28 included_text = strtrim(string(included, "%12.0fc")) + ///
    " (" + strtrim(string(percent_included, "%5.2f")) + ")"
decode sex, generate(sex_label)
decode race, generate(race_label)
generate str4 age_label = cond(age_group == 1, "<75", "75+")
replace cancer_site = "Colon (excluding Appendix)" if cancer_site == "Colon"
replace cancer_site = "Uterine corpus" if cancer_site == "Corpus_and_uterus"
replace cancer_site = "Lung and bronchus" if cancer_site == "Lung_and_bronchus"
replace cancer_site = "Oral cavity and pharynx" if cancer_site == "Oral_cavity_and_pharynx"

putdocx clear
putdocx begin, pagesize(A4) landscape margin(top, 0.4) margin(bottom, 0.4) ///
    margin(left, 0.35) margin(right, 0.35)
putdocx paragraph, style(Title)
putdocx text ("Table 4. Sequential exclusions and final analysis sample")
putdocx table table4 = data(cancer_site sex_label age_label race_label cases_text ///
    excluded_text_1 excluded_text_2 excluded_text_3 excluded_text_4 ///
    excluded_text_5 excluded_text_6 included_text), varnames layout(autofitcontents)
putdocx table table4(1,.), bold
putdocx table table4(.,.), font("Calibri", 7)
putdocx save "$TABLE_OUTPUT_DIR/Table_4.docx", replace

display as result "Table 4 source CSV and Word table created."
