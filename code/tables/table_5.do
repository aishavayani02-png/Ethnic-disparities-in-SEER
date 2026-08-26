/****************************************************************************************
Table 5: five-year total indirect, total direct, and total causal effects, and
proportion mediated. Direct effects and their confidence intervals are retained to
match the current manuscript.
****************************************************************************************/

version 19.0
clear all
set more off
do "code/00_config.do"

import delimited using "$RELEASE_DIR/main/results_long.csv", clear varnames(1) case(preserve)
keep if inlist(estimate, "tie", "tde", "tce", "pm") & ///
    inlist(agegroup, 1, 2) & abs(tt - 5) < 0.000001

generate double estimate_percent = cond(estimate == "pm", PE, 100 * PE)
generate double lower_percent = cond(estimate == "pm", lci, 100 * lci)
generate double upper_percent = cond(estimate == "pm", uci, 100 * uci)
keep site sex agegroup estimate estimate_percent lower_percent upper_percent
sort site sex agegroup estimate
export delimited using "$RELEASE_DIR/manuscript_tables/table_5_source.csv", replace

generate str36 result = strtrim(string(estimate_percent, "%5.2f")) + " (" + ///
    strtrim(string(lower_percent, "%5.2f")) + ", " + ///
    strtrim(string(upper_percent, "%5.2f")) + ")"
keep site sex agegroup estimate result
reshape wide result, i(site sex agegroup) j(estimate) string
rename resulttie TIE
rename resulttde TDE
rename resulttce TCE
rename resultpm proportion_mediated
generate str4 age_label = cond(agegroup == 1, "<75", "75+")
replace site = "Colon (excluding Appendix)" if site == "Colon"
replace site = "Uterine corpus" if site == "Corpus_and_uterus"
replace site = "Lung and bronchus" if site == "Lung_and_bronchus"
replace site = "Oral cavity and pharynx" if site == "Oral_cavity_and_pharynx"

putdocx clear
putdocx begin, pagesize(A4) landscape margin(top, 0.5) margin(bottom, 0.5) ///
    margin(left, 0.5) margin(right, 0.5)
putdocx paragraph, style(Title)
putdocx text ("Table 5. Five-year causal contrasts and proportion mediated")
putdocx table table5 = data(site sex age_label TIE TDE TCE proportion_mediated), ///
    varnames layout(autofitcontents)
putdocx table table5(1,.), bold
putdocx table table5(.,.), font("Calibri", 8)
putdocx save "$TABLE_OUTPUT_DIR/Table_5.docx", replace
