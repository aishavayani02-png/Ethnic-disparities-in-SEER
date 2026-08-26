/****************************************************************************************
Table 3: deaths postponed at 1, 3, and 5 years in a typical 2022 Black incident cohort.
The ADb result files already contain the indirect effect multiplied by the aggregate
2022 Black incident-cohort size; they are not multiplied a second time here.
****************************************************************************************/

version 19.0
clear all
set more off
do "code/00_config.do"

import delimited using "$RELEASE_DIR/combined/Shiny.csv", clear varnames(1) case(preserve)
keep if estimate == "ADb" & inlist(tt, 1, 3, 5)
generate double deaths = PE
generate double deaths_lower = lci
generate double deaths_upper = uci
keep site sex agegroup tt n_black deaths deaths_lower deaths_upper
sort site sex agegroup tt
export delimited using "$RELEASE_DIR/manuscript_tables/table_3_source.csv", replace

generate str36 result = strtrim(string(deaths, "%6.1f")) + " (" + ///
    strtrim(string(deaths_lower, "%6.1f")) + ", " + ///
    strtrim(string(deaths_upper, "%6.1f")) + ")"
generate byte follow_up = round(tt)
keep site sex agegroup n_black follow_up result
reshape wide result, i(site sex agegroup n_black) j(follow_up)
rename result1 one_year
rename result3 three_years
rename result5 five_years

putdocx clear
putdocx begin, pagesize(A4) landscape margin(top, 0.5) margin(bottom, 0.5) ///
    margin(left, 0.5) margin(right, 0.5)
putdocx paragraph, style(Title)
putdocx text ("Table 3. Estimated deaths postponed in a typical 2022 Black incident cohort")
putdocx table table3 = data(site sex agegroup n_black one_year three_years five_years), ///
    varnames layout(autofitcontents)
putdocx table table3(1,.), bold
putdocx table table3(.,.), font("Calibri", 8)
putdocx save "$TABLE_OUTPUT_DIR/Table_3.docx", replace
