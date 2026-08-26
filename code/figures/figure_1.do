version 19.0
clear all
set more off
do "code/00_config.do"

import delimited using "$RELEASE_DIR/combined/Shiny.csv", clear varnames(1) case(preserve)

* Recreate the numeric codes used by the submitted plotting script.
generate byte site_code = .
replace site_code = 1 if site == "Breast"
replace site_code = 2 if site == "Colon (excluding Appendix)"
replace site_code = 3 if site == "Uterine corpus"
replace site_code = 6 if site == "Lung and bronchus"
replace site_code = 7 if site == "Oral cavity and pharynx"
replace site_code = 8 if site == "Ovary"
replace site_code = 9 if site == "Pancreas"
replace site_code = 10 if site == "Rectum"
drop site
rename site_code site

generate byte estimate_code = .
replace estimate_code = 1 if estimate == "AD"
replace estimate_code = 2 if estimate == "ADa"
replace estimate_code = 3 if estimate == "ADb"
replace estimate_code = 4 if estimate == "Model fit"
replace estimate_code = 5 if estimate == "PP"
replace estimate_code = 6 if estimate == "acs"
replace estimate_code = 7 if estimate == "pm"
replace estimate_code = 8 if estimate == "tce"
replace estimate_code = 9 if estimate == "tde"
replace estimate_code = 10 if estimate == "tie"
drop estimate
rename estimate_code estimate

generate byte agegroup_code = cond(agegroup == "<75", 1, 2)
drop agegroup
rename agegroup_code agegroup

generate byte sex_code = cond(sex == "Male", 1, 2)
drop sex
rename sex_code sex
label define sex_label 1 "Male" 2 "Female", replace
label values sex sex_label

label define site 1 "Breast" 2 "Colon (excl. Appendix)" 3 "Corpus uteri" ///
    6 "Lung & bronchus" 7 "Oral cavity & pharynx" 8 "Ovary" 9 "Pancreas" 10 "Rectum", replace
label values site site

* Proportion mediated at 5 years since diagnosis
drop if agegroup == 0
keep if tt == 5
replace PE = PE*100 if estimate > 7 
replace lci = lci*100 if estimate > 7
replace uci = uci*100 if estimate > 7  

label drop site
label define site 1 "Breast" 2  "Colon (excl. Appendix)" 3 "Corpus uteri" ///
 4 "Kidney & renal pelvis" 5 "Liver & IBD" 6 "Lung & bronchus" ///
 7 "Oral cavity & pharynx" 8 "Ovary" 9 "Pancreas" 10 "Rectum"
 
label values site site
drop if site == 4 | site == 5
drop if site == 10 & sex == 2

* Males
preserve
keep if estimate==7 & sex==1

decode site, gen(site_str)
encode site_str, gen(site_id)

gen site_y = site_id
replace site_y = site_id + 0.1 if agegroup==1
replace site_y = site_id - 0.1 if agegroup==2

* confirm the value-label name used by site_id
local vl : value label site_id
*display "`vl'"   // should print something like: site_id

* attach that value label to site_y too
label values site_y `vl'

* Plot
twoway (scatter site_y PE if agegroup == 1, mcol("28 40 94")) ///
 (rcap lci uci site_y if agegroup == 1, hor lcol("28 40 94"%50)) ///
 (scatter site_y PE if agegroup == 2, mcol("252 93 93")) ///
 (rcap lci uci site_y if agegroup == 2, hor lcol("252 93 93"%50)), ///
 yscale(range(0.8 5)) xlabel(0 10(5)50 85 90 100, labsize(medium)) ///
 ylabel(1/5, valuelabel angle(0) noticks nogrid labsize(medium)) ///
 ytitle("") xtitle("Proportion mediated by stage at diagnosis (%)", size(medium)) ///
 legend(order(1 "<75" 3 "75+") ring(0) pos(1) size(medium)) ///
 title("", size(medium)) ///
 plotregion(margin(b=0 r=1 l=0)) graphregion(color(white)) name(Prop_M, replace)
restore


* Females
preserve
keep if estimate==7 & sex==2

decode site, gen(site_str)
encode site_str, gen(site_id)

gen site_y = site_id
replace site_y = site_id + 0.1 if agegroup==1
replace site_y = site_id - 0.1 if agegroup==2

* confirm the value-label name used by site_id
local vl : value label site_id

* attach that value label to site_y too
label values site_y `vl'

* Plot
twoway (scatter site_y PE if agegroup == 1, mcol("28 40 94")) ///
 (rcap lci uci site_y if agegroup == 1, hor lcol("28 40 94"%50)) ///
 (scatter site_y PE if agegroup == 2, mcol("252 93 93")) ///
 (rcap lci uci site_y if agegroup == 2, hor lcol("252 93 93"%50)), ///
 yscale(range(0 7)) xlabel(-10 0(5)20 40(5)75 100, labsize(medium) format(%9.0f)) ///
 ylabel(1/7, valuelabel angle(0) noticks nogrid labsize(medium)) ///
 ytitle("") xtitle("Proportion mediated by stage at diagnosis (%)", size(medium)) ///
 legend(order(1 "<75" 3 "75+") ring(0) pos(1) size(medium)) title("", size(medium)) ///
 plotregion(margin(b=0 r=1 l=0)) graphregion(color(white)) name(Prop_F, replace)
restore


* Scatter plot for TIE at 5 years since diagnosis
* Males
preserve
keep if estimate==10 & sex==1

decode site, gen(site_str)
encode site_str, gen(site_id)

gen site_y = site_id
replace site_y = site_id + 0.1 if agegroup==1
replace site_y = site_id - 0.1 if agegroup==2

* confirm the value-label name used by site_id
local vl : value label site_id
*display "`vl'"   // should print something like: site_id

* attach that value label to site_y too
label values site_y `vl'

* Plot
twoway (scatter site_y PE if agegroup == 1, mcol("28 40 94")) ///
 (rcap lci uci site_y if agegroup == 1, hor lcol("28 40 94"%50)) ///
 (scatter site_y PE if agegroup == 2, mcol("252 93 93")) ///
 (rcap lci uci site_y if agegroup == 2, hor lcol("252 93 93"%50)), ///
 yscale(range(0 5)) xlabel(0(1)4, labsize(medium) format(%9.1f)) ///
 ylabel(1/5, valuelabel angle(0) noticks nogrid labsize(medium)) ///
 ytitle("") xtitle("Difference in standardised all-cause probability of death (%)", size(medium)) ///
 legend(off) title("Indirect effect", size(medium)) ///
 plotregion(margin(b=0 r=1 l=0)) graphregion(color(white)) name(tie_M, replace)
restore


* Females
preserve
keep if estimate==10 & sex==2

decode site, gen(site_str)
encode site_str, gen(site_id)

gen site_y = site_id
replace site_y = site_id + 0.1 if agegroup==1
replace site_y = site_id - 0.1 if agegroup==2

* confirm the value-label name used by site_id
local vl : value label site_id

* attach that value label to site_y too
label values site_y `vl'

* Plot
twoway (scatter site_y PE if agegroup == 1, mcol("28 40 94")) ///
 (rcap lci uci site_y if agegroup == 1, hor lcol("28 40 94"%50)) ///
 (scatter site_y PE if agegroup == 2, mcol("252 93 93")) ///
 (rcap lci uci site_y if agegroup == 2, hor lcol("252 93 93"%50)), ///
 yscale(range(0.8 6.8)) xlabel(-1 0(2)8, labsize(medium) format(%9.1f)) ///
 ylabel(1/7, valuelabel angle(0) noticks nogrid labsize(medium)) ///
 ytitle("") xtitle("Difference in standardised all-cause probability of death (%)", size(medium)) ///
 legend(off) title("Indirect effect", size(medium)) ///
 plotregion(margin(b=0 r=1 l=0)) graphregion(color(white)) name(tie_F, replace)
restore


* Scatter plot for TCE at 5 years since diagnosis
* Males
preserve
keep if estimate==8 & sex==1

decode site, gen(site_str)
encode site_str, gen(site_id)

gen site_y = site_id
replace site_y = site_id + 0.1 if agegroup==1
replace site_y = site_id - 0.1 if agegroup==2

* confirm the value-label name used by site_id
local vl : value label site_id
*display "`vl'"   // should print something like: site_id

* attach that value label to site_y too
label values site_y `vl'

* Plot
twoway (scatter site_y PE if agegroup == 1, mcol("28 40 94")) ///
 (rcap lci uci site_y if agegroup == 1, hor lcol("28 40 94"%50)) ///
 (scatter site_y PE if agegroup == 2, mcol("252 93 93")) ///
 (rcap lci uci site_y if agegroup == 2, hor lcol("252 93 93"%50)), ///
 yscale(range(0 5)) xlabel(0(2)14, labsize(medium) format(%9.1f)) ///
 ylabel(none) ///
 ytitle("") xtitle("Difference in standardised all-cause probability of death (%)", size(medium)) ///
 legend(off) title("Total causal effect", size(medium)) ///
 plotregion(margin(b=0 r=1 l=0)) graphregion(color(white)) name(tce_M, replace)
restore


* Females
preserve
keep if estimate==8 & sex==2

decode site, gen(site_str)
encode site_str, gen(site_id)

gen site_y = site_id
replace site_y = site_id + 0.1 if agegroup==1
replace site_y = site_id - 0.1 if agegroup==2

* confirm the value-label name used by site_id
local vl : value label site_id

* attach that value label to site_y too
label values site_y `vl'

* Plot
twoway (scatter site_y PE if agegroup == 1, mcol("28 40 94")) ///
 (rcap lci uci site_y if agegroup == 1, hor lcol("28 40 94"%50)) ///
 (scatter site_y PE if agegroup == 2, mcol("252 93 93")) ///
 (rcap lci uci site_y if agegroup == 2, hor lcol("252 93 93"%50)), ///
 yscale(range(0.8 6.8)) xlabel(0(2)8 18 , labsize(medium) format(%9.1f)) ///
 ylabel(none) ytitle("") ///
 xtitle("Difference in standardised all-cause probability of death (%)", size(medium)) ///
 legend(off) title("Total causal effect", size(medium)) ///
 plotregion(margin(b=0 r=1 l=0)) graphregion(color(white)) name(tce_F, replace)
restore

graph combine tie_M tce_M, imargin(0 0 0 0) cols(2) name(te_M,replace) ///
plotregion(margin(zero)) graphregion(col(white)) title("") 

graph combine tie_F tce_F, imargin(0 1 0 0) cols(2) name(te_F,replace) ///
plotregion(margin(zero)) graphregion(col(white)) title("")

graph combine te_M Prop_M, imargin(b=0 r=0 l=0) rows(2) name(causal_M,replace) ///
plotregion(margin(zero)) graphregion(col(white)) title("Males", size(medium) span) 

graph export "$FIGURE_OUTPUT_DIR/Figure_1_males.png", replace

graph combine te_F Prop_F, imargin(b=0 r=0 l=0) rows(2) name(causal_F,replace) ///
plotregion(margin(zero)) graphregion(col(white)) title("Females", size(medium) span) ///
note("Plots exclude estimates of female rectum due to wide confidence intervals and/or negative proportions, presented in Table 5", size(vsmall) span color("60 60 60"))

graph export "$FIGURE_OUTPUT_DIR/Figure_1_females.png", replace
