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

* The original combined plotting file inherited these transformations from Figure 1.
* State them explicitly now that Figure 2 is a standalone script.
keep if abs(tt - 5) < 0.000001
replace PE = PE * 100 if estimate > 7
replace lci = lci * 100 if estimate > 7
replace uci = uci * 100 if estimate > 7


* Plot for deaths postponed

*-----------------------------------------------
* Black people count by agegroup/site (Females) 
*-----------------------------------------------
forvalues i = 1/2 {
	preserve
	keep if sex == `i'
	local sexlabel : label (sex) `i'
	drop if missing(n_black)
	keep if estimate==1 

	* Make it one row per site x agegroup (to mimic graph hbar (sum))
	collapse (sum) n_black, by(site agegroup)

	* Ensure consistent ordering
	sort site agegroup

	* 1=<75, 2=75+  (adjust if needed)
	bysort site (agegroup): gen byte within_site = _n
	egen site_id = group(site), label

	local GAP = 1
	gen double y = (site_id-1)*(2+`GAP') + within_site

	* ---- Build ylabel() macro from site value labels ----
	local ylbl ""
	local sitelab : value label site

	quietly levelsof site_id, local(sites)
	foreach sid of local sites {

		quietly summarize y if site_id==`sid' & within_site==1, meanonly
		local y0 = r(mean)

		quietly levelsof site if site_id==`sid' & within_site==1, local(sitecode)
		local sc : word 1 of `sitecode'

		if "`sitelab'" != "" {
			local nm : label `sitelab' `sc'
		}
		else {
			local nm "`sc'"
		}

		local ylbl `"`ylbl' `y0' "`nm'""'
	}

	if `i' == 1 {
		twoway (bar n_black y if agegroup==1, horizontal barw(0.75) ///
				fcolor("28 40 94") lcolor("28 40 94") lwidth(medthick)) ///
				(bar n_black y if agegroup==2, horizontal barw(0.75) ///
				fcolor("252 93 93") lcolor("252 93 93") lwidth(medthick)), ///
				yscale(reverse) ///
				ylabel(`ylbl', labsize(medium) noticks nogrid) ///
				ytitle("") ///
				title("Number of Black diagnoses in 2022", size(medium)) ///
				xlabel(0(1000)3000, format(%9.0f) labsize(medium) grid) ///
				xtitle("Counts (N)", size(medium)) legend(off) ///
				plotregion(margin(0 r+0.5 0 0)) ///
				graphregion(color(white)) ///
				name(Black_`sexlabel', replace)
	}

	if `i' == 2{
		twoway (bar n_black y if agegroup==1, horizontal barw(0.75) ///
				fcolor("28 40 94") lcolor("28 40 94") lwidth(medthick)) ///
				(bar n_black y if agegroup==2, horizontal barw(0.75) ///
				fcolor("252 93 93") lcolor("252 93 93") lwidth(medthick)), ///
				yscale(reverse) ///
				ylabel(`ylbl', labsize(medium) noticks nogrid) ///
				ytitle("") ///
				title("Number of Black diagnoses in 2022", size(medium)) ///
				xlabel(0(4000)12000, format(%9.0f) labsize(medium) grid) ///
				xtitle("Counts (N)", size(medium)) legend(off) ///
				plotregion(margin(0 0 0 0)) ///
				graphregion(color(white)) ///
				name(Black_`sexlabel', replace)
	}
	restore
}


*---------------------------------------------
* Total indirect effect (TIE) by agegroup/site
*---------------------------------------------
forvalues i = 1/2 {
	preserve
	keep if sex == `i'
	local sexlabel : label (sex) `i'
	keep if estimate==10 
	
	* Ensure consistent ordering
	sort site agegroup

	* 1=<75, 2=75+ (adjust if needed)
	bysort site (agegroup): gen byte within_site = _n

	* numeric id per site (keeps current order)
	egen site_id = group(site), label

	* gap size between sites (tune this)
	local GAP = 1
	gen double y = (site_id-1)*(2+`GAP') + within_site

	* get readable site names
	capture confirm value label site
	if _rc==0 {
		decode site, gen(site_name)
	}
	else {
		tostring site, gen(site_name) usedisplayformat
	}

	* Build ylabel macro: label only the first row of each site (within_site==1)
	* ---- Build ylabel() macro with *value labels* for site ----
	local ylbl ""

	* get the name of the value label attached to site (if any)
	local sitelab : value label site

	quietly levelsof site_id, local(sites)
	foreach sid of local sites {

		* y position for the first row of this site (within_site==1)
		quietly summarize y if site_id==`sid' & within_site==1, meanonly
		local y0 = r(mean)

		* get the numeric site code for this site_id
		quietly levelsof site if site_id==`sid' & within_site==1, local(sitecode)
		local sc : word 1 of `sitecode'

		* get the display name
		if "`sitelab'" != "" {
			local nm : label `sitelab' `sc'
		}
		else {
			local nm "`sc'"
		}

		* append to ylabel macro (quotes important!)
		local ylbl `"`ylbl' `y0' "`nm'""'
	}

	* (optional) display to check
	di `"`ylbl'"'

	if `i' == 1 {
		twoway (bar PE y if agegroup==1, horizontal barw(0.75) ///
				fcolor(white) lcolor("28 40 94") lwidth(medthick)) ///
				(bar PE y if agegroup==2, horizontal barw(0.75) ///
				fcolor(white) lcolor("252 93 93") lwidth(medthick)) ///
				(rcap lci uci y if agegroup==1, horizontal lcolor("28 40 94") lwidth(medium)) ///
				(rcap lci uci y if agegroup==2, horizontal lcolor("252 93 93") ///
				lwidth(medium)), ///
				yscale(reverse) ytitle("") xlabel(,format(%9.2f) labsize(medium)) ///
				title("Total indirect effect", size(medium)) ///
				xtitle("Difference in standardised all-cause probability of death (%)", ///
				size(medium)) ///
				ylabel(none) ///
				legend(off) ///
				plotregion(margin(0 0 0 0)) ///
				graphregion(color(white)) ///
				name(tie_`sexlabel', replace)
	}
	if `i' == 2 {
		twoway (bar PE y if agegroup==1, horizontal barw(0.75) ///
				fcolor(white) lcolor("28 40 94") lwidth(medthick)) ///
				(bar PE y if agegroup==2, horizontal barw(0.75) ///
				fcolor(white) lcolor("252 93 93") lwidth(medthick)) ///
				(rcap lci uci y if agegroup==1, horizontal lcolor("28 40 94") ///
				lwidth(medium)) ///
				(rcap lci uci y if agegroup==2, horizontal lcolor("252 93 93") ///
				lwidth(medium)), ///
				yscale(reverse) ytitle("") ///
				ylabel(none) xlabel(, format(%9.2f) labsize(medium)) ///
				title("Total indirect effect", size(medium)) ///
				xtitle("Difference in standardised all-cause probability of death (%)", ///
				size(medium)) ///
				legend(off) ///
				plotregion(margin(l+1 0 0 0)) ///
				graphregion(color(white)) ///
				name(tie_`sexlabel', replace)
	}
			
		restore
}

graph combine Black_Male tie_Male, ycommon name(Male,replace) ///
    imargin(0 r+1 0 0) graphregion(margin(0 0 0 0))
graph combine Black_Female tie_Female, ycommon name(Female,replace) ///
    imargin(0 r+2 0 0)

*-----------------------------------------------------------------
* Avoidable deaths due to stage differences (ADb) by agegroup/site
*-----------------------------------------------------------------
forvalues i = 1/2 {
	preserve
	keep if sex == `i'
	local sexlabel : label (sex) `i'
	keep if estimate==3 
	
	* Ensure consistent ordering
	sort site agegroup

	* 1=<75, 2=75+ (adjust if needed)
	bysort site (agegroup): gen byte within_site = _n

	* numeric id per site (keeps current order)
	egen site_id = group(site), label

	* gap size between sites (tune this)
	local GAP = 0.5
	gen double y = (site_id-1)*(2+`GAP') + within_site

	* get readable site names
	capture confirm value label site
	if _rc==0 {
		decode site, gen(site_name)
	}
	else {
		tostring site, gen(site_name) usedisplayformat
	}

	* Build ylabel macro: label only the first row of each site (within_site==1)
	* ---- Build ylabel() macro with *value labels* for site ----
	local ylbl ""

	* get the name of the value label attached to site (if any)
	local sitelab : value label site

	quietly levelsof site_id, local(sites)
	foreach sid of local sites {

		* y position for the first row of this site (within_site==1)
		quietly summarize y if site_id==`sid' & within_site==1, meanonly
		local y0 = r(mean)

		* get the numeric site code for this site_id
		quietly levelsof site if site_id==`sid' & within_site==1, local(sitecode)
		local sc : word 1 of `sitecode'

		* get the display name
		if "`sitelab'" != "" {
			local nm : label `sitelab' `sc'
		}
		else {
			local nm "`sc'"
		}

		* append to ylabel macro (quotes important!)
		local ylbl `"`ylbl' `y0' "`nm'""'
	}

	* (optional) display to check
	di `"`ylbl'"'

	if `i' == 1 {
		twoway (bar PE y if agegroup==1, horizontal barw(0.75) ///
				fcolor(white) lcolor("28 40 94") lwidth(medthick)) ///
				(bar PE y if agegroup==2, horizontal barw(0.75) ///
				fcolor(white) lcolor("252 93 93") lwidth(medthick)) ///
				(rcap lci uci y if agegroup==1, horizontal lcolor("28 40 94") lwidth(medium)) ///
				(rcap lci uci y if agegroup==2, horizontal lcolor("252 93 93") ///
				lwidth(medium)), ///
				yscale(reverse) ///
				ylabel(`ylbl', labsize(medium) noticks nogrid) ///
				ytitle("") ///
				xlabel(0(5)15 30(5)60,format(%9.0f) labsize(medium)) ///
				xtitle("Postponable deaths at 5 years since diagnosis", size(medium)) ///
				legend(order(1 "<75" 2 "75+") pos(5) ring(0) rows(2) size(medium)) ///
				plotregion(margin(0 r+1 0 0)) ///
				graphregion(color(white)) ///
				name(AD_`sexlabel', replace)
	}
	if `i' == 2 {
		twoway (bar PE y if agegroup==1, horizontal barw(0.75) ///
				fcolor(white) lcolor("28 40 94") lwidth(medthick)) ///
				(bar PE y if agegroup==2, horizontal barw(0.75) ///
				fcolor(white) lcolor("252 93 93") lwidth(medthick)) ///
				(rcap lci uci y if agegroup==1, horizontal lcolor("28 40 94") lwidth(medium)) ///
				(rcap lci uci y if agegroup==2, horizontal lcolor("252 93 93") ///
				lwidth(medium)), ///
				yscale(reverse) ///
				ylabel(`ylbl', labsize(medium) noticks nogrid) ///
				ytitle("") ///
				xlabel(0(25)50 200 225 280, format(%9.0f) labsize(medium)) ///
				xtitle("Postponable deaths at 5 years since diagnosis", size(medium)) ///
				legend(order(1 "<75" 2 "75+") pos(5) ring(0) rows(2) size(medium)) ///
				plotregion(margin(l+1 0 0 0)) ///
				graphregion(color(white)) ///
				name(AD_`sexlabel', replace)
	}		
		restore
}

graph combine Male AD_Male, ycommon title("Males", size(medium)) ///
rows(2) name(Male,replace) plotregion(margin(zero)) graphregion(col(white)) 

*graph export Male.png, replace
graph combine Female AD_Female, ycommon title("Females", size(medium)) ///
rows(2) name(Female,replace) plotregion(margin(zero)) graphregion(col(white)) 

*graph export Female.png, replace

graph combine Male Female, ycommon name(AD,replace) cols(2) ///
    imargin(0 0 0 0) graphregion(col(white)) plotregion(margin(zero)) ///
	title("Avoidable deaths due to differences in stage at diagnosis at 5 years since diagnosis", size(medium) span color("60 60 60")) ///
	note("Plots exclude estimates of female rectum due to wide confidence intervals and/or negative proportions, presented in Table 5", size(vsmall) span color("60 60 60")) 
	
graph export "$FIGURE_OUTPUT_DIR/Figure_2.png", replace
