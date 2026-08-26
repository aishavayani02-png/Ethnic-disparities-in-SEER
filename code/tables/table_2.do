/****************************************************************************************
Table 2: five-year age-standardised all-cause survival by race, sex, age group,
and stage at diagnosis.
****************************************************************************************/

version 19.0
clear all
set more off
do "code/00_config.do"

import delimited using "$RELEASE_DIR/combined/Shiny.csv", clear varnames(1) case(preserve)
keep if estimate == "acs" & abs(tt - 5) < 0.000001
generate double estimate_percent = 100 * PE
generate double lower_percent = 100 * lci
generate double upper_percent = 100 * uci
keep site sex agegroup race stage estimate_percent lower_percent upper_percent
sort site agegroup race sex stage
export delimited using "$RELEASE_DIR/manuscript_tables/table_2_source.csv", replace

generate str32 result = strtrim(string(estimate_percent, "%4.1f")) + " (" + ///
    strtrim(string(lower_percent, "%4.1f")) + ", " + ///
    strtrim(string(upper_percent, "%4.1f")) + ")"
keep site sex agegroup race stage result
reshape wide result, i(site sex agegroup race) j(stage) string
rename resultLocalised localised
rename resultRegional regional
rename resultDistant distant

putdocx clear
putdocx begin, pagesize(A4) landscape margin(top, 0.5) margin(bottom, 0.5) ///
    margin(left, 0.5) margin(right, 0.5)
putdocx paragraph, style(Title)
putdocx text ("Table 2. Five-year age-standardised all-cause survival (%)")
putdocx table table2 = data(site sex agegroup race localised regional distant), ///
    varnames layout(autofitcontents)
putdocx table table2(1,.), bold
putdocx table table2(.,.), font("Calibri", 8)
putdocx save "$TABLE_OUTPUT_DIR/Table_2.docx", replace
