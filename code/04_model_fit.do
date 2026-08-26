/****************************************************************************************
Model-fit assessment: compare fitted flexible parametric models with the
non-parametric Pohar Perme estimator, by stage, race, age group, and sex.
****************************************************************************************/

version 19.0
clear all
set more off
do "code/00_config.do"
set seed $ANALYSIS_SEED

if $PATHS_CONFIGURED != 1 {
	display as error "Edit PRIVATE_DATA_ROOT in code/00_config.do before running."
	exit 198
}

* master dataset for RS by stage and age group across all sites - stratified by sex
tempfile Shiny
save `Shiny', emptyok replace


** for all sites 
foreach site in  "Breast" "Colon" "Corpus_and_uterus" "Lung_and_bronchus" ///
"Oral_cavity_and_pharynx" "Ovary" "Pancreas" "Rectum" { 
	
	frame change default
	* load dataset
	cd "$CLEAN_DATA_DIR"
	use `site'.dta, clear 
	levelsof sex if !missing(sex), local(sexlist)
	
	foreach j of local sexlist { 
		
		cd "$CLEAN_DATA_DIR"
		use `site'.dta, clear 
		
		keep if sex == `j'
		local sexlabel : label (sex) `j'
		
		gen id = _n 
		
		* Randomly assigned dates for all 
		gen month_diag = runiformint(1,12) 
		gen day_diag = runiformint(1,28) 
		gen datediag_random = mdy(month_diag, day_diag, year) 
		
		* Exit date based on surv_mm 
		gen dateexit = dofm(mofd(datediag_random)+surv_mm) 
		
		* Systematic diagnosis dates for censored cases 
		gen datediag_cens = dofm(mofd(dateexit)-surv_mm) 
		
		* Combine systematic and random 
		gen datediag = cond(dead==1, datediag_random, datediag_cens) 
		
		* Clean up 
		drop month_diag day_diag datediag_random datediag_cens 
		
		stset dateexit, enter(time mdy(1,1,2021)) exit(time min(datediag + ///
		5.1*365.24, mdy(12, 31, 2022))) origin(datediag) failure(dead==1) ///
		scale(365.24) id(id) 
		
		* Merge in expected survival using population lifetable (popmort.dta) 
		rename age agediag 
		rename year yeardiag 
		
		* Create categorical age groups
		gen agegrp = .
		replace agegrp=1 if agediag<75 
		replace agegrp=2 if agediag>74
		replace agegrp=4 if agegrp==.
		 
		egen groups = group(race stage agegrp), label 
		tempfile gmap 
			preserve
				keep groups race stage agegrp sex 
				duplicates drop 
				save `gmap', replace 
			restore 
			
		// PP estimator
		stpp R_pp using popmort, by(groups) ///
		agediag(agediag) datediag(datediag) /// 
		pmother(sex race) pmage(age) /// 
		pmyear(year) pmmaxyear(2021) pmmaxage(99) ///
		list(1 (1) 5) graphname(R_pp, replace) frame(stpp,replace)  
		
		gen age = min(int(agediag+_t),99) 
		gen year = min(floor(yeardiag + _t),2021) 
		sort sex year age race 
		
		merge m:1 sex year age race using popmort.dta, ///
		keep(match master) keepusing(rate) 
		
		drop _merge 
		
		* Fit the survival model using race, stage and age 
		stpm3 i.race i.stage @rcs(agediag, df(3) winsor(2.5 97.5)) i.race#i.stage ///
		i.race#@rcs(agediag, df(3) winsor(2.5 97.5)) ///
		i.stage#@rcs(agediag, df(3) winsor(2.5 97.5)),  df(3) ///
		tvc(i.race @rcs(agediag, df(3) winsor(2.5 97.5)) i.stage i.race#i.stage ///
		i.race#@rcs(agediag, df(3) winsor(2.5 97.5))) dftvc(2) ///
		scale(lncumhazard) bhaz(rate) eform neq(1) initvaluesloop iter(100) 
		
		* Generate time variable for predictions 
		keep if _t0 == 0 
		range tt 0 5 101 
		
		// Model fit check (RS scale) - age groups
		standsurv , surv timevar(tt) over(groups) ci ///
		frame(agegrp,replace) atvars(RS_groups_*) 
		
		frame agegrp{
			
			* ---- Model fit (RS scale) ----
			keep tt RS_groups_*

			* reshape long: groups index + PE/lci/uci
			reshape long RS_groups_ RS_groups_@_lci RS_groups_@_uci, i(tt) j(groups)

			rename RS_groups_        PE
			rename RS_groups__lci    lci
			rename RS_groups__uci    uci

			* attach race, stage, agegrp, sex (numeric + labels likely)
			merge m:1 groups using `gmap', nogen

			* agegroup like rs_long_all
			rename agegrp agegroup
			gen site = "`site'"
			gen estimate = "Model fit"
			
			tempfile rs_this
			save `rs_this', replace
		}
		
		frame stpp{
			
			keep time PP PP_lci PP_uci groups
			rename time   tt
			rename PP     PE
			rename PP_lci lci
			rename PP_uci uci
			gen site = "`site'"
			gen estimate = "PP"

			* attach race, stage, agegrp, sex
			merge m:1 groups using `gmap', nogen
			rename agegrp agegroup
			
			tempfile pp
			save `pp', replace
			
		}
		
		// Model fit check (RS scale) - marginal
		cd "$CLEAN_DATA_DIR"
		use `site'.dta, clear 
		
		keep if sex == `j'
		local sexlabel : label (sex) `j'
		
		gen id = _n 
		
		* Randomly assigned dates for all 
		gen month_diag = runiformint(1,12) 
		gen day_diag = runiformint(1,28) 
		gen datediag_random = mdy(month_diag, day_diag, year) 
		
		* Exit date based on surv_mm 
		gen dateexit = dofm(mofd(datediag_random)+surv_mm) 
		
		* Systematic diagnosis dates for censored cases 
		gen datediag_cens = dofm(mofd(dateexit)-surv_mm) 
		
		* Combine systematic and random 
		gen datediag = cond(dead==1, datediag_random, datediag_cens) 
		
		* Clean up 
		drop month_diag day_diag datediag_random datediag_cens 
		
		stset dateexit, enter(time mdy(1,1,2021)) exit(time min(datediag + ///
		5.1*365.24, mdy(12, 31, 2022))) origin(datediag) failure(dead==1) ///
		scale(365.24) id(id) 
		
		* Merge in expected survival using population lifetable (popmort.dta) 
		rename age agediag 
		rename year yeardiag 
		
		* Create categorical age groups
		gen agegrp = 0
		 
		egen groups = group(race stage agegrp), label 
		tempfile gmap 
			preserve
				keep groups race stage agegrp sex 
				duplicates drop 
				save `gmap', replace 
			restore 
			
		// PP estimator
		stpp R_pp using popmort, by(groups) ///
		agediag(agediag) datediag(datediag) /// 
		pmother(sex race) pmage(age) /// 
		pmyear(year) pmmaxyear(2021) pmmaxage(99) ///
		list(1 (1) 5) graphname(R_pp, replace) frame(stpp_all,replace)  
		
		gen age = min(int(agediag+_t),99) 
		gen year = min(floor(yeardiag + _t),2021) 
		sort sex year age race 
		
		merge m:1 sex year age race using popmort.dta, ///
		keep(match master) keepusing(rate) 
		
		drop _merge 
		
		* Fit the survival model using race, stage and age 
		stpm3 i.race i.stage @rcs(agediag, df(3) winsor(2.5 97.5)) i.race#i.stage ///
		i.race#@rcs(agediag, df(3) winsor(2.5 97.5)) ///
		i.stage#@rcs(agediag, df(3) winsor(2.5 97.5)),  df(3) ///
		tvc(i.race @rcs(agediag, df(3) winsor(2.5 97.5)) i.stage i.race#i.stage ///
		i.race#@rcs(agediag, df(3) winsor(2.5 97.5))) dftvc(2) ///
		scale(lncumhazard) bhaz(rate) eform neq(1) initvaluesloop iter(100) 
		
		* Generate time variable for predictions 
		keep if _t0 == 0 
		range tt 0 5 101 
		
		// Model fit check (RS scale) - age groups
		standsurv , surv timevar(tt) over(groups) ci ///
		frame(all,replace) atvars(RS_groups_*) 
		
		frame all{
			
			* ---- Model fit (RS scale) ----
			keep tt RS_groups_*

			* reshape long: groups index + PE/lci/uci
			reshape long RS_groups_ RS_groups_@_lci RS_groups_@_uci, i(tt) j(groups)

			rename RS_groups_        PE
			rename RS_groups__lci    lci
			rename RS_groups__uci    uci

			* attach race, stage, agegrp, sex (numeric + labels likely)
			merge m:1 groups using `gmap', nogen

			* agegroup like rs_long_all
			rename agegrp agegroup
			gen site = "`site'"
			gen estimate = "Model fit"
			
			tempfile rs_all
			save `rs_all', replace
		}
		
		frame stpp_all{
			
			keep time PP PP_lci PP_uci groups
			rename time   tt
			rename PP     PE
			rename PP_lci lci
			rename PP_uci uci
			gen site = "`site'"
			gen estimate = "PP"

			* attach race, stage, agegrp, sex
			merge m:1 groups using `gmap', nogen
			rename agegrp agegroup
			
			tempfile pp_all
			save `pp_all', replace
			
		}


* Now append both into your running tempfile `Shiny`
use `Shiny', clear
append using `rs_this'
append using `pp'
append using `rs_all'
append using `pp_all'
save `Shiny', replace

		}
		
}

use `Shiny', clear
order site sex estimate race stage agegroup tt PE lci uci
cd "$MODEL_FIT_DIR"
save model_fit.dta, replace

gen str6 sex_str = ""
replace sex_str = "Male"   if sex == 1
replace sex_str = "Female" if sex == 2
drop sex
rename sex_str sex

gen str6 race_str = ""
replace race_str = "White"   if race == 1
replace race_str = "Black" if race == 2
drop race
rename race_str race

gen str6 stage_str = ""
replace stage_str = "Localised"   if stage == 1
replace stage_str = "Regional" if stage == 2
replace stage_str = "Distant" if stage == 3
drop stage
rename stage_str stage

replace agegroup = . if agegroup == 0

drop groups

order site sex estimate race stage agegroup tt PE lci uci
cd "$MODEL_FIT_DIR"
save model_fit.dta, replace
