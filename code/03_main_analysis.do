/****************************************************************************************
Main analysis: standardised survival, causal mediation effects, and deaths postponed.

Run this file from the repository root after editing code/00_config.do. The script
saves one private working .dta file per estimate/site/sex/age group. Use
code/06_combine_results.do to assemble those files and create disclosure-safe CSVs.

The uncertainty calculation reproduces the submitted analysis: mediator-model
parameters are simulated while the fitted survival model is held fixed.
****************************************************************************************/

version 19.0
clear all
set more off
do "code/00_config.do"

if $PATHS_CONFIGURED != 1 {
	display as error "Edit PRIVATE_DATA_ROOT in code/00_config.do before running."
	exit 198
}

capture pr drop sampbeta
program sampbeta, 
  syntax name
  tempname b V
  matrix `b' = e(b)
  matrix `V' = e(V)
  local Np = colsof(`b')
  local cnames:  colfullnames `b'
  local rnames:  rowfullnames `b'  
  forvalues i = 1/`Np' {
    tempvar p`i'
    local plist `plist' `p`i''
  }
  tempname m
  frame create `m'
  frame `m' {
    drawnorm `plist', mean(`b') cov(`V') n(1)
    mkmat `plist', matrix(`namelist')
  }
  matrix colnames `namelist' = `cnames'
  matrix rownames `namelist' = `rnames'  
  erepost b = `namelist'  V=`V', noesample 
end

* ALL-CAUSE SURVIVAL MODEL
***************************
* Reset the seed for each site because the original analyses were run site by site.
foreach site of global ANALYSIS_SITES {
	set seed $ANALYSIS_SEED
	
	frame change default
	* load dataset
	cd "$CLEAN_DATA_DIR"
	use `site'.dta, clear 
	levelsof sex if !missing(sex), local(sexlist) 
		
	foreach j of local sexlist {
		frame reset
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
		
		stset dateexit, enter(time mdy(1,1,2021)) exit(time min(datediag+5.1*365.24, ///
		mdy(12, 31, 2022))) origin(datediag) failure(dead==1) scale(365.24) id(id)
		
		* Merge in expected survival using population lifetable (popmort.dta)
		rename age agediag
		rename year yeardiag
		gen age = min(int(agediag+_t),99)
		gen year = min(floor(yeardiag + _t),2021)
		sort sex year age race
		
		merge  m:1 sex year age race using  popmort.dta, keep(match master) keepusing(rate)
		drop _merge
		
		
		// STEP 1
		**********
		* Fit the survival model using race, stage and age 
		stpm3 i.race i.stage @rcs(agediag, df(3) winsor(2.5 97.5)) i.race#i.stage ///
		i.race#@rcs(agediag, df(3) winsor(2.5 97.5)) ///
		i.stage#@rcs(agediag, df(3) winsor(2.5 97.5)),  df(3) ///
		tvc(i.race @rcs(agediag, df(3) winsor(2.5 97.5)) i.stage i.race#i.stage ///
		i.race#@rcs(agediag, df(3) winsor(2.5 97.5))) dftvc(2) ///
		scale(lncumhazard) bhaz(rate) eform neq(1) initvaluesloop iter(100) 
		
		estimates store surv

		
		// STEP 2
		*********
		* Fit a model for the mediator including the exposure and confounders.
		
		//Generate splines for age
		rcsgen agediag, df(3) gen(rcsa) orthog
		* Generates orthogonal restricted cubic spline variables for agediag with 3df
		
		local ageknots `r(knots)' // store knot locations into local macro ageknots
		matrix R = r(R) // save orthogonalisation matrix into R
		matrix list R // display content of matrix R

		//Fit a multinomial regression model for White ethnicity
		quietly: mlogit stage rcsa* if race==1 
		estimates store ph1
			
		//Fit a multinomial regression model for Black ethnicity
		quietly: mlogit stage rcsa* if race==2 
		estimates store ph0
		
		
		// STEP 3
		*********
		* For each individual, obtain predictions for the prob of being in a specific level
		** of the mediator (stage) at each level of the exposure (ethnicity).

		// Point estimates
		// For Black ethnicity
		estimates restore ph0
		gen tmprace = race // new variable with same observations as race column
		replace race = 2 // make every obs in race column of black ethnicity
		predict p11 p12 p13 // predicted probabilities for each stage (3 categories)
		replace race = tmprace // return race column to how it was
		drop tmprace

		// For White ethnicity
		estimates restore ph1
		capture drop tmprace
		gen tmprace = race
		replace race = 1
		predict p01 p02 p03 
		replace race = tmprace
		drop tmprace
		
		
		// STEP 4
		*********		
		* Generate time variable for predictions
		keep if _t0==0
		range tt 0 5 101
		
		estimates restore surv
		frame change default
		
		//Standardised relative survival by stage and race
		capture drop _at*
		standsurv if race==2, surv timevar(tt) ci ///
		at1(race 1 stage 1) ///
		at2(race 1 stage 2) ///
		at3(race 1 stage 3) ///
		at4(race 2 stage 1) ///
		at5(race 2 stage 2) ///
		at6(race 2 stage 3) ///
		atvars(lw rw dw lb rb db) ///
		frame(rs_`sexlabel',mergecreate) ///
		expsurv(using (popmort.dta) ///
				datediag(datediag)      ///
				agediag(agediag)        ///
				pmrate(rate)            ///
				pmage(age)              ///
				pmyear(year)            ///
				pmother(race sex)       ///
				pmmaxyear(2021)         ///
				at1(.)                  ///
				at2(.)                  ///
				at3(.)                  ///
				at4(.)                  ///
				at5(.)                  ///
				at6(.))
		
		// TCE
		capture drop _at* 
		standsurv if race==2, timevar(tt) failure  ///
		at1(race 2 stage 1, atindweights(p11)) ///
		at2(race 2 stage 2, atindweights(p12)) ///
		at3(race 2 stage 3, atindweights(p13)) ///
		at4(race 1 stage 1, atindweights(p01)) ///
		at5(race 1 stage 2, atindweights(p02)) ///
		at6(race 1 stage 3, atindweights(p03)) ///
		lincom(1 1 1 -1 -1 -1) lincomvar(PE)	///
		frame(tce_`sexlabel',mergecreate)		///
		expsurv(using (popmort.dta) ///
			datediag(datediag)      ///
			agediag(agediag)        ///
			pmrate(rate)            ///
			pmage(age)              ///
			pmyear(year)            ///
			pmother(race sex)       ///
			pmmaxyear(2021)         ///
			at1(.)                  ///
			at2(.)                  ///
			at3(.)                  ///
			at4(.)                  ///
			at5(.)                  ///
			at6(.)) 

		// NDE
		capture drop _at* 
		standsurv if race==2, failure timevar(tt)  ///
		at1(race 2 stage 1, atindweights(p01)) ///
		at2(race 2 stage 2, atindweights(p02)) ///
		at3(race 2 stage 3, atindweights(p03)) ///
		at4(race 1 stage 1, atindweights(p01)) ///
		at5(race 1 stage 2, atindweights(p02)) ///
		at6(race 1 stage 3, atindweights(p03)) ///
		lincom(1 1 1 -1 -1 -1) lincomvar(PE) ///
		frame(tde_`sexlabel',mergecreate) ///
		expsurv(using (popmort.dta) ///
			datediag(datediag)      ///
			agediag(agediag)        ///
			pmrate(rate)            ///
			pmage(age)              ///
			pmyear(year)            ///
			pmother(race sex)       ///
			pmmaxyear(2021)         ///
			at1(.)                  ///
			at2(.)                  ///
			at3(.)                  ///
			at4(.)                  ///
			at5(.)                  ///
			at6(.)) 
				
		// NIE
		capture drop _at*
		standsurv if race==2, failure timevar(tt)  ///
		at1(race 2 stage 1, atindweights(p11)) ///
		at2(race 2 stage 2, atindweights(p12)) ///
		at3(race 2 stage 3, atindweights(p13)) ///
		at4(race 2 stage 1, atindweights(p01)) ///
		at5(race 2 stage 2, atindweights(p02)) ///
		at6(race 2 stage 3, atindweights(p03)) ///
		lincom(1 1 1 -1 -1 -1) lincomvar(PE) ///
		frame(tie_`sexlabel',mergecreate) ///
		expsurv(using (popmort.dta) ///
			datediag(datediag)      ///
			agediag(agediag)        ///
			pmrate(rate)            ///
			pmage(age)              ///
			pmyear(year)            ///
			pmother(race sex)       ///
			pmmaxyear(2021)         ///
			at1(.)                  ///
			at2(.)                  ///
			at3(.)                  ///
			at4(.)                  ///
			at5(.)                  ///
			at6(.))
				
		* N* = total number of Black patients diagnosed in 2022
		** (inclusive of those with missing stage)
		local Black = .
		if "`site'" == "Breast" local Black = 12054
		else if "`site'" == "Colon" {
			if `j' == 1 local Black = 2347
			else if `j' == 2 local Black = 2383
		}
		else if "`site'" == "Corpus_and_uterus" local Black = 3326
		else if "`site'" == "Lung_and_bronchus" {
			if `j' == 1 local Black = 4207
			else if `j' == 2 local Black = 3870
		}
		else if "`site'" == "Oral_cavity_and_pharynx" {
			if `j' == 1 local Black = 1034
			else if `j' == 2 local Black = 475 
		}
		else if "`site'" == "Ovary" local Black = 948
		else if "`site'" == "Pancreas" {
			if `j' == 1 local Black = 1287 
			else if `j' == 2 local Black = 1489 
		}
		else if "`site'" == "Rectum" {
			if `j' == 1 local Black = 1122 
			else if `j' == 2 local Black = 911 
		}
		else if "`site'" == "Stomach" {
			if `j' == 1 local Black = 942 
			else if `j' == 2 local Black = 929 
		}
		else if "`site'" == "Bladder" {
			if `j' == 1 local Black = 736 
			else if `j' == 2 local Black = 356 
		}
		
		// Total avoidable deaths
		capture drop _at*
		standsurv if race==2, timevar(tt) failure verbose per(`Black')         ///
		at1(race 2 stage 1, atindweights(p11)) ///
		at2(race 2 stage 2, atindweights(p12)) ///
		at3(race 2 stage 3, atindweights(p13)) ///
		at4(race 1 stage 1, atindweights(p01)) ///
		at5(race 1 stage 2, atindweights(p02)) ///
		at6(race 1 stage 3, atindweights(p03)) ///
		lincom(1 1 1 -1 -1 -1) lincomvar(PE) ///
		frame(AD_`sexlabel',mergecreate) ///
		expsurv(using (popmort.dta) ///
			datediag(datediag)      ///
			agediag(agediag)        ///
			pmrate(rate)            ///
			pmage(age)              ///
			pmyear(year)            ///
			pmother(race sex)       ///
			pmmaxyear(2021)         ///
			at1(.)                  ///
			at2(.)                  ///
			at3(.)                  ///
			at4(.)                  ///
			at5(.)                  ///
			at6(.)) 
		
		// Avoidable deaths by eliminating non-stage differences
		capture drop _at* 
		standsurv if race==2, failure timevar(tt) verbose per(`Black')  ///
		at1(race 2 stage 1, atindweights(p01)) ///
		at2(race 2 stage 2, atindweights(p02)) ///
		at3(race 2 stage 3, atindweights(p03)) ///
		at4(race 1 stage 1, atindweights(p01)) ///
		at5(race 1 stage 2, atindweights(p02)) ///
		at6(race 1 stage 3, atindweights(p03)) ///
		lincom(1 1 1 -1 -1 -1) lincomvar(PE) ///
		frame(ADa_`sexlabel',mergecreate) ///
		expsurv(using (popmort.dta) ///
			datediag(datediag)      ///
			agediag(agediag)        ///
			pmrate(rate)            ///
			pmage(age)              ///
			pmyear(year)            ///
			pmother(race sex)       ///
			pmmaxyear(2021)         ///
			at1(.)                  ///
			at2(.)                  ///
			at3(.)                  ///
			at4(.)                  ///
			at5(.)                  ///
			at6(.)) 
				
		// Avoidable deaths by eliminating stage differences
		capture drop _at*
		standsurv if race==2, failure timevar(tt) verbose per(`Black')  ///
		at1(race 2 stage 1, atindweights(p11)) ///
		at2(race 2 stage 2, atindweights(p12)) ///
		at3(race 2 stage 3, atindweights(p13)) ///
		at4(race 2 stage 1, atindweights(p01)) ///
		at5(race 2 stage 2, atindweights(p02)) ///
		at6(race 2 stage 3, atindweights(p03)) ///
		lincom(1 1 1 -1 -1 -1) lincomvar(PE) ///
		frame(ADb_`sexlabel',mergecreate) ///
		expsurv(using (popmort.dta) ///
			datediag(datediag)      ///
			agediag(agediag)        ///
			pmrate(rate)            ///
			pmage(age)              ///
			pmyear(year)            ///
			pmother(race sex)       ///
			pmmaxyear(2021)         ///
			at1(.)                  ///
			at2(.)                  ///
			at3(.)                  ///
			at4(.)                  ///
			at5(.)                  ///
			at6(.)) 
				
			
		* Repeat above by age group
		
		* Create categorical age groups
		gen agegrp = .
		replace agegrp=1 if agediag<75 
		replace agegrp=2 if agediag>74
		
		forvalues k = 1/2 {
			
			//Standardised relative survival by stage and race
			capture drop _at*
			standsurv if race==2 & agegrp == `k', surv ci timevar(tt) ///
			at1(race 1 stage 1) ///
			at2(race 1 stage 2) ///
			at3(race 1 stage 3) ///
			at4(race 2 stage 1) ///
			at5(race 2 stage 2) ///
			at6(race 2 stage 3) ///
			atvars(lw rw dw lb rb db) ///
			frame(rs_`k'_`sexlabel',mergecreate) ///
			expsurv(using (popmort.dta) ///
					datediag(datediag)      ///
					agediag(agediag)        ///
					pmrate(rate)            ///
					pmage(age)              ///
					pmyear(year)            ///
					pmother(race sex)       ///
					pmmaxyear(2021)         ///
					at1(.)                  ///
					at2(.)                  ///
					at3(.)                  ///
					at4(.)                  ///
					at5(.)                  ///
					at6(.)) 
			
			// TCE
			capture drop _at*
			standsurv if race==2 & agegrp==`k', timevar(tt) failure            ///
			at1(race 2 stage 1, atindweights(p11)) ///
			at2(race 2 stage 2, atindweights(p12)) ///
			at3(race 2 stage 3, atindweights(p13)) ///
			at4(race 1 stage 1, atindweights(p01)) ///
			at5(race 1 stage 2, atindweights(p02)) ///
			at6(race 1 stage 3, atindweights(p03)) ///
			lincom(1 1 1 -1 -1 -1) lincomvar(PE) ///
			frame(tce_`k'_`sexlabel',mergecreate) ///
			expsurv(using (popmort.dta) ///
					datediag(datediag)      ///
					agediag(agediag)        ///
					pmrate(rate)            ///
					pmage(age)              ///
					pmyear(year)            ///
					pmother(race sex)       ///
					pmmaxyear(2021)         ///
					at1(.)                  ///
					at2(.)                  ///
					at3(.)                  ///
					at4(.)                  ///
					at5(.)                  ///
					at6(.)) 
			
			// NDE
			capture drop _at* 
			standsurv if race==2 & agegrp==`k', failure timevar(tt)  ///
			at1(race 2 stage 1, atindweights(p01)) ///
			at2(race 2 stage 2, atindweights(p02)) ///
			at3(race 2 stage 3, atindweights(p03)) ///
			at4(race 1 stage 1, atindweights(p01)) ///
			at5(race 1 stage 2, atindweights(p02)) ///
			at6(race 1 stage 3, atindweights(p03)) ///
			lincom(1 1 1 -1 -1 -1) lincomvar(PE) ///
			frame(tde_`k'_`sexlabel',mergecreate) ///
			expsurv(using (popmort.dta) ///
					datediag(datediag)      ///
					agediag(agediag)        ///
					pmrate(rate)            ///
					pmage(age)              ///
					pmyear(year)            ///
					pmother(race sex)       ///
					pmmaxyear(2021)         ///
					at1(.)                  ///
					at2(.)                  ///
					at3(.)                  ///
					at4(.)                  ///
					at5(.)                  ///
					at6(.)) 
				
			// NIE
			capture drop _at*
			standsurv if race==2 & agegrp==`k', failure timevar(tt)  ///
			at1(race 2 stage 1, atindweights(p11)) ///
			at2(race 2 stage 2, atindweights(p12)) ///
			at3(race 2 stage 3, atindweights(p13)) ///
			at4(race 2 stage 1, atindweights(p01)) ///
			at5(race 2 stage 2, atindweights(p02)) ///
			at6(race 2 stage 3, atindweights(p03)) ///
			lincom(1 1 1 -1 -1 -1) lincomvar(PE) ///
			frame(tie_`k'_`sexlabel',mergecreate) ///
			expsurv(using (popmort.dta) ///
					datediag(datediag)      ///
					agediag(agediag)        ///
					pmrate(rate)            ///
					pmage(age)              ///
					pmyear(year)            ///
					pmother(race sex)       ///
					pmmaxyear(2021)         ///
					at1(.)                  ///
					at2(.)                  ///
					at3(.)                  ///
					at4(.)                  ///
					at5(.)                  ///
					at6(.)) 
				
			* N* = total number of Black patients diagnosed in 2022 within age group
			** (inclusive of those with missing stage)
			local Black = .
			if "`site'" == "Breast" {
				if `k' == 1 local Black = 10309
				else if `k' == 2 local Black = 1745
			}
			else if "`site'" == "Colon" {
				if `j' == 1 & `k' == 1 local Black = 1918
				else if `j' == 1 & `k' == 2 local Black = 429
				else if `j' == 2 & `k' == 1 local Black = 1802
				else if `j' == 2 & `k' == 2 local Black = 581
			}
			else if "`site'" == "Corpus_and_uterus" {
				if `k' == 1 local Black = 2843
				else if `k' == 2 local Black = 483
			}
			else if "`site'" == "Lung_and_bronchus" {
				if `j' == 1 & `k' == 1 local Black = 3188
				else if `j' == 1 & `k' == 2 local Black = 1019
				else if `j' == 2 & `k' == 1 local Black = 2757
				else if `j' == 2 & `k' == 2 local Black = 1113
			}
			else if "`site'" == "Oral_cavity_and_pharynx" {
				if `j' == 1 & `k' == 1 local Black = 902
				else if `j' == 1 & `k' == 2 local Black = 132
				else if `j' == 2 & `k' == 1 local Black = 408
				else if `j' == 2 & `k' == 2 local Black = 67
			}
			else if "`site'" == "Ovary" {
				if `k' == 1 local Black = 808 
				else if `k' == 2 local Black = 140
			}
			else if "`site'" == "Pancreas" {
				if `j' == 1 & `k' == 1 local Black = 980
				else if `j' == 1 & `k' == 2 local Black = 307
				else if `j' == 2 & `k' == 1 local Black = 1051
				else if `j' == 2 & `k' == 2 local Black = 438
			}
			else if "`site'" == "Rectum" {
				if `j' == 1 & `k' == 1 local Black = 991
				else if `j' == 1 & `k' == 2 local Black = 131
				else if `j' == 2 & `k' == 1 local Black = 788
				else if `j' == 2 & `k' == 2 local Black = 123
			}
			else if "`site'" == "Stomach" {
				if `j' == 1 & `k' == 1 local Black = 746
				else if `j' == 1 & `k' == 2 local Black = 196
				else if `j' == 2 & `k' == 1 local Black = 668
				else if `j' == 2 & `k' == 2 local Black = 261
			}
			else if "`site'" == "Bladder" {
				if `j' == 1 & `k' == 1 local Black = 480
				else if `j' == 1 & `k' == 2 local Black = 256
				else if `j' == 2 & `k' == 1 local Black = 211
				else if `j' == 2 & `k' == 2 local Black = 145
			}
		
			// Total avoidable deaths
			capture drop _at*
			standsurv if race==2 & agegrp==`k', timevar(tt) failure ///
			verbose per(`Black')         ///
			at1(race 2 stage 1, atindweights(p11)) ///
			at2(race 2 stage 2, atindweights(p12)) ///
			at3(race 2 stage 3, atindweights(p13)) ///
			at4(race 1 stage 1, atindweights(p01)) ///
			at5(race 1 stage 2, atindweights(p02)) ///
			at6(race 1 stage 3, atindweights(p03)) ///
			lincom(1 1 1 -1 -1 -1) lincomvar(PE) ///
			frame(AD_`k'_`sexlabel',mergecreate) ///
			expsurv(using (popmort.dta) ///
					datediag(datediag)      ///
					agediag(agediag)        ///
					pmrate(rate)            ///
					pmage(age)              ///
					pmyear(year)            ///
					pmother(race sex)       ///
					pmmaxyear(2021)         ///
					at1(.)                  ///
					at2(.)                  ///
					at3(.)                  ///
					at4(.)                  ///
					at5(.)                  ///
					at6(.)) 

		
			// Avoidable deaths by eliminating non-stage differences
			capture drop _at* 
			standsurv if race==2 & agegrp==`k', failure timevar(tt) verbose per(`Black')  ///
			at1(race 2 stage 1, atindweights(p01)) ///
			at2(race 2 stage 2, atindweights(p02)) ///
			at3(race 2 stage 3, atindweights(p03)) ///
			at4(race 1 stage 1, atindweights(p01)) ///
			at5(race 1 stage 2, atindweights(p02)) ///
			at6(race 1 stage 3, atindweights(p03)) ///
			lincom(1 1 1 -1 -1 -1) lincomvar(PE) ///
			frame(ADa_`k'_`sexlabel',mergecreate) ///
			expsurv(using (popmort.dta) ///
					datediag(datediag)      ///
					agediag(agediag)        ///
					pmrate(rate)            ///
					pmage(age)              ///
					pmyear(year)            ///
					pmother(race sex)       ///
					pmmaxyear(2021)         ///
					at1(.)                  ///
					at2(.)                  ///
					at3(.)                  ///
					at4(.)                  ///
					at5(.)                  ///
					at6(.)) 

				
			// Avoidable deaths by eliminating stage differences
			capture drop _at*
			standsurv if race==2 & agegrp==`k', failure timevar(tt) verbose per(`Black')  ///
			at1(race 2 stage 1, atindweights(p11)) ///
			at2(race 2 stage 2, atindweights(p12)) ///
			at3(race 2 stage 3, atindweights(p13)) ///
			at4(race 2 stage 1, atindweights(p01)) ///
			at5(race 2 stage 2, atindweights(p02)) ///
			at6(race 2 stage 3, atindweights(p03)) ///
			lincom(1 1 1 -1 -1 -1) lincomvar(PE) ///
			frame(ADb_`k'_`sexlabel',mergecreate) ///
			expsurv(using (popmort.dta) ///
					datediag(datediag)      ///
					agediag(agediag)        ///
					pmrate(rate)            ///
					pmage(age)              ///
					pmyear(year)            ///
					pmother(race sex)       ///
					pmmaxyear(2021)         ///
					at1(.)                  ///
					at2(.)                  ///
					at3(.)                  ///
					at4(.)                  ///
					at5(.)                  ///
					at6(.)) 
								
		}


		// STEP 5
		*********
		* Parametric bootstrapping for CIs

		// Repeat the parametric simulation 200 times, as reported in the manuscript.
		global m $PARAMETRIC_DRAWS
		forvalues i = 1/$m {

		// Weights for Black ethnicity 
		preserve
		estimates restore ph0
		sampbeta b 
		restore

		capture drop tmprace
		capture drop p11_`i' p12_`i' p13_`i' 
		gen tmprace = race
		replace race = 2 
		predict p11_`i' p12_`i' p13_`i' 
		replace race = tmprace
		drop tmprace

		// Weights for White ethnicity 
		preserve
		estimates restore ph1
		sampbeta b
		restore
		
		capture drop tmprace
		capture drop p01_`i' p02_`i' p03_`i'
		gen tmprace = race
		replace race = 1 
		predict p01_`i' p02_`i' p03_`i'
		replace race = tmprace
		drop tmprace


		// Same as before for the survival model
		preserve
		estimates restore surv
		matrix bsurv = e(b)
		matrix VVVsurv= e(V) 
		sampbeta b	
		restore
		

		// Obtain predictions for i bootstrap sample
		frame tce_`sexlabel':  capture drop _at*
		frame tie_`sexlabel': capture drop _at*
		frame tde_`sexlabel': capture drop _at*
		frame AD_`sexlabel': capture drop _at*
		frame ADa_`sexlabel': capture drop _at*
		frame ADb_`sexlabel': capture drop _at*
		
		estimates restore surv
		
		// TCE
		capture drop _at*
		standsurv if race==2, timevar(tt) failure verbose     ///
		at1(race 2 stage 1, atindweights(p11_`i')) ///
		at2(race 2 stage 2, atindweights(p12_`i')) ///
		at3(race 2 stage 3, atindweights(p13_`i')) ///
		at4(race 1 stage 1, atindweights(p01_`i')) ///
		at5(race 1 stage 2, atindweights(p02_`i')) ///
		at6(race 1 stage 3, atindweights(p03_`i')) ///
		lincom(1 1 1 -1 -1 -1) lincomvar(tce_`i') ///
		frame(tce_`sexlabel',merge)  ///
		expsurv(using (popmort.dta) ///
				datediag(datediag)      ///
				agediag(agediag)        ///
				pmrate(rate)            ///
				pmage(age)              ///
				pmyear(year)            ///
				pmother(race sex)       ///
				pmmaxyear(2021)         ///
				at1(.)                  ///
				at2(.)                  ///
				at3(.)                  ///
				at4(.)                  ///
				at5(.)                  ///
				at6(.))  
				
		// NDE
		capture drop _at*
		standsurv if race==2, failure timevar(tt) verbose            ///
		at1(race 2 stage 1, atindweights(p01_`i')) ///
		at2(race 2 stage 2, atindweights(p02_`i')) ///
		at3(race 2 stage 3, atindweights(p03_`i')) ///
		at4(race 1 stage 1, atindweights(p01_`i')) ///
		at5(race 1 stage 2, atindweights(p02_`i')) ///
		at6(race 1 stage 3, atindweights(p03_`i')) ///
		lincom(1 1 1 -1 -1 -1) lincomvar(tde_`i') ///
		frame(tde_`sexlabel',merge) ///
		expsurv(using (popmort.dta) ///
				datediag(datediag)      ///
				agediag(agediag)        ///
				pmrate(rate)            ///
				pmage(age)              ///
				pmyear(year)            ///
				pmother(race sex)       ///
				pmmaxyear(2021)         ///
				at1(.)                  ///
				at2(.)                  ///
				at3(.)                  ///
				at4(.)                  ///
				at5(.)                  ///
				at6(.)) 
				
		// NIE
		capture drop _at* 
		standsurv if race==2, failure timevar(tt) verbose           ///
		at1(race 2 stage 1, atindweights(p11_`i')) ///
		at2(race 2 stage 2, atindweights(p12_`i')) ///
		at3(race 2 stage 3, atindweights(p13_`i')) ///
		at4(race 2 stage 1, atindweights(p01_`i')) ///
		at5(race 2 stage 2, atindweights(p02_`i')) ///
		at6(race 2 stage 3, atindweights(p03_`i')) ///
		lincom(1 1 1 -1 -1 -1) lincomvar(tie_`i') ///
		frame(tie_`sexlabel',merge) ///
		expsurv(using (popmort.dta) ///
				datediag(datediag)      ///
				agediag(agediag)        ///
				pmrate(rate)            ///
				pmage(age)              ///
				pmyear(year)            ///
				pmother(race sex)       ///
				pmmaxyear(2021)         ///
				at1(.)                  ///
				at2(.)                  ///
				at3(.)                  ///
				at4(.)                  ///
				at5(.)                  ///
				at6(.)) 

		// AD
		* N* = total number of Black patients diagnosed in 2022
		** (inclusive of those with missing stage)
		local Black = .
		if "`site'" == "Breast" local Black = 12054
		else if "`site'" == "Colon" {
			if `j' == 1 local Black = 2347
			else if `j' == 2 local Black = 2383
		}
		else if "`site'" == "Corpus_and_uterus" local Black = 3326
		else if "`site'" == "Lung_and_bronchus" {
			if `j' == 1 local Black = 4207
			else if `j' == 2 local Black = 3870
		}
		else if "`site'" == "Oral_cavity_and_pharynx" {
			if `j' == 1 local Black = 1034
			else if `j' == 2 local Black = 475 
		}
		else if "`site'" == "Ovary" local Black = 948
		else if "`site'" == "Pancreas" {
			if `j' == 1 local Black = 1287 
			else if `j' == 2 local Black = 1489 
		}
		else if "`site'" == "Rectum" {
			if `j' == 1 local Black = 1122 
			else if `j' == 2 local Black = 911 
		}
		else if "`site'" == "Stomach" {
			if `j' == 1 local Black = 942 
			else if `j' == 2 local Black = 929 
		}
		else if "`site'" == "Bladder" {
			if `j' == 1 local Black = 736 
			else if `j' == 2 local Black = 356 
		}
		
		// Total avoidable deaths
		capture drop _at*
		standsurv if race==2, timevar(tt) failure verbose per(`Black')   ///
		at1(race 2 stage 1, atindweights(p11_`i')) ///
		at2(race 2 stage 2, atindweights(p12_`i')) ///
		at3(race 2 stage 3, atindweights(p13_`i')) ///
		at4(race 1 stage 1, atindweights(p01_`i')) ///
		at5(race 1 stage 2, atindweights(p02_`i')) ///
		at6(race 1 stage 3, atindweights(p03_`i')) ///
		lincom(1 1 1 -1 -1 -1) lincomvar(AD_`i') ///
		frame(AD_`sexlabel',merge)  ///
		expsurv(using (popmort.dta) ///
				datediag(datediag)      ///
				agediag(agediag)        ///
				pmrate(rate)            ///
				pmage(age)              ///
				pmyear(year)            ///
				pmother(race sex)       ///
				pmmaxyear(2021)         ///
				at1(.)                  ///
				at2(.)                  ///
				at3(.)                  ///
				at4(.)                  ///
				at5(.)                  ///
				at6(.))  
		
		// Avoidable deaths by removing other relative survival differences
		capture drop _at*
		standsurv if race==2, failure timevar(tt) verbose per(`Black')            ///
		at1(race 2 stage 1, atindweights(p01_`i')) ///
		at2(race 2 stage 2, atindweights(p02_`i')) ///
		at3(race 2 stage 3, atindweights(p03_`i')) ///
		at4(race 1 stage 1, atindweights(p01_`i')) ///
		at5(race 1 stage 2, atindweights(p02_`i')) ///
		at6(race 1 stage 3, atindweights(p03_`i')) ///
		lincom(1 1 1 -1 -1 -1) lincomvar(ADa_`i') ///
		frame(ADa_`sexlabel',merge)  ///
		expsurv(using (popmort.dta) ///
				datediag(datediag)      ///
				agediag(agediag)        ///
				pmrate(rate)            ///
				pmage(age)              ///
				pmyear(year)            ///
				pmother(race sex)       ///
				pmmaxyear(2021)         ///
				at1(.)                  ///
				at2(.)                  ///
				at3(.)                  ///
				at4(.)                  ///
				at5(.)                  ///
				at6(.)) 
			
		// Avoidable deaths by removing stage differences
		standsurv if race==2, failure timevar(tt) verbose per(`Black')          ///
		at1(race 2 stage 1, atindweights(p11_`i')) ///
		at2(race 2 stage 2, atindweights(p12_`i')) ///
		at3(race 2 stage 3, atindweights(p13_`i')) ///
		at4(race 2 stage 1, atindweights(p01_`i')) ///
		at5(race 2 stage 2, atindweights(p02_`i')) ///
		at6(race 2 stage 3, atindweights(p03_`i')) ///
		lincom(1 1 1 -1 -1 -1) lincomvar(ADb_`i') ///
		frame(ADb_`sexlabel',merge)  ///
		expsurv(using (popmort.dta) ///
				datediag(datediag)      ///
				agediag(agediag)        ///
				pmrate(rate)            ///
				pmage(age)              ///
				pmyear(year)            ///
				pmother(race sex)       ///
				pmmaxyear(2021)         ///
				at1(.)                  ///
				at2(.)                  ///
				at3(.)                  ///
				at4(.)                  ///
				at5(.)                  ///
				at6(.)) 
		
		
		* by age group
		
		forvalues k = 1/2 {
						
			frame tce_`k'_`sexlabel':  capture drop _at*
			frame tie_`k'_`sexlabel': capture drop _at*
			frame tde_`k'_`sexlabel': capture drop _at*

			frame AD_`k'_`sexlabel':  capture drop _at*
			frame ADa_`k'_`sexlabel': capture drop _at*
			frame ADb_`k'_`sexlabel': capture drop _at*
		
			estimates restore surv
			
			// TCE
			capture drop _at*
			standsurv if race==2 & agegrp==`k', timevar(tt) failure verbose     ///
			at1(race 2 stage 1, atindweights(p11_`i')) ///
			at2(race 2 stage 2, atindweights(p12_`i')) ///
			at3(race 2 stage 3, atindweights(p13_`i')) ///
			at4(race 1 stage 1, atindweights(p01_`i')) ///
			at5(race 1 stage 2, atindweights(p02_`i')) ///
			at6(race 1 stage 3, atindweights(p03_`i')) ///
			lincom(1 1 1 -1 -1 -1) lincomvar(tce_`i') ///
			frame(tce_`k'_`sexlabel',merge) ///
			expsurv(using (popmort.dta) ///
					datediag(datediag)      ///
					agediag(agediag)        ///
					pmrate(rate)            ///
					pmage(age)              ///
					pmyear(year)            ///
					pmother(race sex)       ///
					pmmaxyear(2021)         ///
					at1(.)                  ///
					at2(.)                  ///
					at3(.)                  ///
					at4(.)                  ///
					at5(.)                  ///
					at6(.)) 
					
			// NDE
			capture drop _at*
			standsurv if race==2 & agegrp==`k', failure timevar(tt) verbose            ///
			at1(race 2 stage 1, atindweights(p01_`i')) ///
			at2(race 2 stage 2, atindweights(p02_`i')) ///
			at3(race 2 stage 3, atindweights(p03_`i')) ///
			at4(race 1 stage 1, atindweights(p01_`i')) ///
			at5(race 1 stage 2, atindweights(p02_`i')) ///
			at6(race 1 stage 3, atindweights(p03_`i')) ///
			lincom(1 1 1 -1 -1 -1) lincomvar(tde_`i') ///
			frame(tde_`k'_`sexlabel',merge) ///
			expsurv(using (popmort.dta) 	///
					datediag(datediag)      ///
					agediag(agediag)        ///
					pmrate(rate)            ///
					pmage(age)              ///
					pmyear(year)            ///
					pmother(race sex)       ///
					pmmaxyear(2021)         ///
					at1(.)                  ///
					at2(.)                  ///
					at3(.)                  ///
					at4(.)                  ///
					at5(.)                  ///
					at6(.)) 
			
			// NDE
			capture drop _at* 
			standsurv if race==2 & agegrp==`k', failure timevar(tt) verbose           ///
			at1(race 2 stage 1, atindweights(p11_`i')) ///
			at2(race 2 stage 2, atindweights(p12_`i')) ///
			at3(race 2 stage 3, atindweights(p13_`i')) ///
			at4(race 2 stage 1, atindweights(p01_`i')) ///
			at5(race 2 stage 2, atindweights(p02_`i')) ///
			at6(race 2 stage 3, atindweights(p03_`i')) ///
			lincom(1 1 1 -1 -1 -1) lincomvar(tie_`i') ///
			frame(tie_`k'_`sexlabel',merge) ///
			expsurv(using (popmort.dta) ///
					datediag(datediag)      ///
					agediag(agediag)        ///
					pmrate(rate)            ///
					pmage(age)              ///
					pmyear(year)            ///
					pmother(race sex)       ///
					pmmaxyear(2021)         ///
					at1(.)                  ///
					at2(.)                  ///
					at3(.)                  ///
					at4(.)                  ///
					at5(.)                  ///
					at6(.)) 
					

			* N* = total number of Black patients diagnosed in 2022 within age group
			** (inclusive of those with missing stage)
			local Black = .
			if "`site'" == "Breast" {
				if `k' == 1 local Black = 10309
				else if `k' == 2 local Black = 1745
			}
			else if "`site'" == "Colon" {
				if `j' == 1 & `k' == 1 local Black = 1918
				else if `j' == 1 & `k' == 2 local Black = 429
				else if `j' == 2 & `k' == 1 local Black = 1802
				else if `j' == 2 & `k' == 2 local Black = 581
			}
			else if "`site'" == "Corpus_and_uterus" {
				if `k' == 1 local Black = 2843
				else if `k' == 2 local Black = 483
			}
			else if "`site'" == "Lung_and_bronchus" {
				if `j' == 1 & `k' == 1 local Black = 3188
				else if `j' == 1 & `k' == 2 local Black = 1019
				else if `j' == 2 & `k' == 1 local Black = 2757
				else if `j' == 2 & `k' == 2 local Black = 1113
			}
			else if "`site'" == "Oral_cavity_and_pharynx" {
				if `j' == 1 & `k' == 1 local Black = 902
				else if `j' == 1 & `k' == 2 local Black = 132
				else if `j' == 2 & `k' == 1 local Black = 408
				else if `j' == 2 & `k' == 2 local Black = 67
			}
			else if "`site'" == "Ovary" {
				if `k' == 1 local Black = 808 
				else if `k' == 2 local Black = 140
			}
			else if "`site'" == "Pancreas" {
				if `j' == 1 & `k' == 1 local Black = 980
				else if `j' == 1 & `k' == 2 local Black = 307
				else if `j' == 2 & `k' == 1 local Black = 1051
				else if `j' == 2 & `k' == 2 local Black = 438
			}
			else if "`site'" == "Rectum" {
				if `j' == 1 & `k' == 1 local Black = 991
				else if `j' == 1 & `k' == 2 local Black = 131
				else if `j' == 2 & `k' == 1 local Black = 788
				else if `j' == 2 & `k' == 2 local Black = 123
			}
			else if "`site'" == "Stomach" {
				if `j' == 1 & `k' == 1 local Black = 746
				else if `j' == 1 & `k' == 2 local Black = 196
				else if `j' == 2 & `k' == 1 local Black = 668
				else if `j' == 2 & `k' == 2 local Black = 261
			}
			else if "`site'" == "Bladder" {
				if `j' == 1 & `k' == 1 local Black = 480
				else if `j' == 1 & `k' == 2 local Black = 256
				else if `j' == 2 & `k' == 1 local Black = 211
				else if `j' == 2 & `k' == 2 local Black = 145
			}
		
			// Total AD
			capture drop _at*
			standsurv if race==2 & agegrp==`k', timevar(tt) failure verbose per(`Black')   ///
			at1(race 2 stage 1, atindweights(p11_`i')) ///
			at2(race 2 stage 2, atindweights(p12_`i')) ///
			at3(race 2 stage 3, atindweights(p13_`i')) ///
			at4(race 1 stage 1, atindweights(p01_`i')) ///
			at5(race 1 stage 2, atindweights(p02_`i')) ///
			at6(race 1 stage 3, atindweights(p03_`i')) ///
			lincom(1 1 1 -1 -1 -1) lincomvar(AD_`i') ///
			frame(AD_`k'_`sexlabel',merge)  ///
			expsurv(using (popmort.dta) ///
					datediag(datediag)      ///
					agediag(agediag)        ///
					pmrate(rate)            ///
					pmage(age)              ///
					pmyear(year)            ///
					pmother(race sex)       ///
					pmmaxyear(2021)         ///
					at1(.)                  ///
					at2(.)                  ///
					at3(.)                  ///
					at4(.)                  ///
					at5(.)                  ///
					at6(.))  
		
			// Avoidable deaths by removing other relative survival differences
			capture drop _at*
			standsurv if race==2 & agegrp==`k', failure timevar(tt) ///
			verbose per(`Black')            ///
			at1(race 2 stage 1, atindweights(p01_`i')) ///
			at2(race 2 stage 2, atindweights(p02_`i')) ///
			at3(race 2 stage 3, atindweights(p03_`i')) ///
			at4(race 1 stage 1, atindweights(p01_`i')) ///
			at5(race 1 stage 2, atindweights(p02_`i')) ///
			at6(race 1 stage 3, atindweights(p03_`i')) ///
			lincom(1 1 1 -1 -1 -1) lincomvar(ADa_`i') ///
			frame(ADa_`k'_`sexlabel',merge)  ///
			expsurv(using (popmort.dta) ///
					datediag(datediag)      ///
					agediag(agediag)        ///
					pmrate(rate)            ///
					pmage(age)              ///
					pmyear(year)            ///
					pmother(race sex)       ///
					pmmaxyear(2021)         ///
					at1(.)                  ///
					at2(.)                  ///
					at3(.)                  ///
					at4(.)                  ///
					at5(.)                  ///
					at6(.)) 
			
		

			// Avoidable deaths by removing stage differences
			capture drop _at* 
			standsurv if race==2 & agegrp==`k', failure timevar(tt) ///
			verbose per(`Black')          ///
			at1(race 2 stage 1, atindweights(p11_`i')) ///
			at2(race 2 stage 2, atindweights(p12_`i')) ///
			at3(race 2 stage 3, atindweights(p13_`i')) ///
			at4(race 2 stage 1, atindweights(p01_`i')) ///
			at5(race 2 stage 2, atindweights(p02_`i')) ///
			at6(race 2 stage 3, atindweights(p03_`i')) ///
			lincom(1 1 1 -1 -1 -1) lincomvar(ADb_`i') ///
			frame(ADb_`k'_`sexlabel',merge)  ///
			expsurv(using (popmort.dta) ///
					datediag(datediag)      ///
					agediag(agediag)        ///
					pmrate(rate)            ///
					pmage(age)              ///
					pmyear(year)            ///
					pmother(race sex)       ///
					pmmaxyear(2021)         ///
					at1(.)                  ///
					at2(.)                  ///
					at3(.)                  ///
					at4(.)                  ///
					at5(.)                  ///
					at6(.)) 
			
		}

		capture drop _at*
		drop p11_`i' p12_`i' p13_`i' p01_`i' p02_`i' p03_`i'
				
		}


		// Step 6
		*********
		// Standardised rs
		frame rs_`sexlabel' {
			cd "$MAIN_RESULT_DIR"
			save rs_`site'_`sexlabel'.dta, replace
		}
		
		
		* Calculate 95% CI using the sd of estimates obtained from bootstrap samples
		local z = invnormal(1-(1-c(level)/100)/2)
		
		        // Proportion mediated via stage: (TIE/TCE)*100
        // Uses bootstrap draws: pm_i = (tie_i / tce_i)*100, CI via +/- z*SD as per other estimates
        capture frame drop pm_`sexlabel'
        frame copy tce_`sexlabel' pm_`sexlabel', replace

        frame pm_`sexlabel' {
            // pull in TIE point estimate and bootstrap draws by tt
            frlink 1:1 tt, frame(tie_`sexlabel')
            frget PE, from(tie_`sexlabel') prefix(tie_)
            frget tie_*, from(tie_`sexlabel')

            // point estimate (note: PE here is TCE PE because we copied from tce frame)
            capture drop pm_PE pm_all_sd pm_all_lci pm_all_uci pm_*
            gen pm_PE = (tie_PE / PE) * 100

            // bootstrap draws of proportion mediated
            forvalues i = 1/$m {
                gen pm_`i' = (tie_`i' / tce_`i') * 100
            }

            egen pm_all_sd  = rowsd(pm_*)
            gen  pm_all_lci = pm_PE - `z'*pm_all_sd
            gen  pm_all_uci = pm_PE + `z'*pm_all_sd

            keep tt pm_PE pm_all_lci pm_all_uci
            replace pm_PE=0 if tt==0
            replace pm_all_lci=0 if tt==0
            replace pm_all_uci=0 if tt==0

            cd "$MAIN_RESULT_DIR"
            save pm_`site'_`sexlabel'.dta, replace
        }

		
		// Total causal effect
		frame tce_`sexlabel' {
			capture drop tce_all_sd tce_all_lci tce_all_uci
			
			egen tce_all_sd = rowsd(tce_*)
			gen  tce_all_lci = PE - `z'*tce_all_sd
			gen  tce_all_uci = PE + `z'*tce_all_sd
			
			keep tt PE tce_all_lci tce_all_uci
			replace PE=0 if tt==0
			replace tce_all_lci=0 if tt==0
			replace tce_all_uci=0 if tt==0
			
			cd "$MAIN_RESULT_DIR"
			save tce_`site'_`sexlabel'.dta,replace
		}

		// Direct effect
		frame tde_`sexlabel' {
			capture drop tde_all_sd tde_all_lci tde_all_uci
			
			egen tde_all_sd = rowsd(tde_*)
			gen  tde_all_lci = PE - `z'*tde_all_sd
			gen  tde_all_uci = PE + `z'*tde_all_sd
			
			keep tt PE tde_all_lci tde_all_uci
			replace PE=0 if tt==0
			replace tde_all_lci=0 if tt==0
			replace tde_all_uci=0 if tt==0
			
			cd "$MAIN_RESULT_DIR"
			save tde_`site'_`sexlabel'.dta,replace
		}
		
		// Indirect effect
		frame tie_`sexlabel' {
			capture drop tie_all_sd tie_all_lci tie_all_uci
			
			egen tie_all_sd = rowsd(tie_*)
			gen  tie_all_lci = PE - `z'*tie_all_sd
			gen  tie_all_uci = PE + `z'*tie_all_sd
			
			keep tt PE tie_all_lci tie_all_uci
			replace PE=0 if tt==0
			replace tie_all_lci=0 if tt==0
			replace tie_all_uci=0 if tt==0
			
			cd "$MAIN_RESULT_DIR"
			save tie_`site'_`sexlabel'.dta,replace
		}
		
		// Total avoidable deaths
		frame AD_`sexlabel' {
			capture drop AD_all_sd AD_all_lci AD_all_uci
			
			egen AD_all_sd = rowsd(AD_*)
			gen  AD_all_lci = PE - `z'*AD_all_sd
			gen  AD_all_uci = PE + `z'*AD_all_sd

			keep tt PE AD_all_lci AD_all_uci
			replace PE=0 if tt==0
			replace AD_all_lci=0 if tt==0
			replace AD_all_uci=0 if tt==0
		
			cd "$MAIN_RESULT_DIR"
			save AD_`site'_`sexlabel'.dta,replace		
		}
		
		// AD by removing relative survival differences
		frame ADa_`sexlabel' {
			capture drop ADa_all_sd ADa_all_lci ADa_all_uci
			
			egen ADa_all_sd = rowsd(ADa_*)
			gen  ADa_all_lci = PE - `z'*ADa_all_sd
			gen  ADa_all_uci = PE + `z'*ADa_all_sd

			keep tt PE ADa_all_lci ADa_all_uci
			replace PE=0 if tt==0
			replace ADa_all_lci=0 if tt==0
			replace ADa_all_uci=0 if tt==0
			
			cd "$MAIN_RESULT_DIR"
			save ADa_`site'_`sexlabel'.dta,replace
		}
		
		// AD by removing stage differences
		frame ADb_`sexlabel' {
			capture drop ADb_all_sd ADb_all_lci ADb_all_uci
			
			egen ADb_all_sd = rowsd(ADb_*)
			gen  ADb_all_lci = PE - `z'*ADb_all_sd
			gen  ADb_all_uci = PE + `z'*ADb_all_sd

			keep tt PE ADb_all_lci ADb_all_uci
			replace PE=0 if tt==0
			replace ADb_all_lci=0 if tt==0
			replace ADb_all_uci=0 if tt==0
			
			cd "$MAIN_RESULT_DIR"
			save ADb_`site'_`sexlabel'.dta,replace
		}

		forvalues k=1/2 {
			
			// Standardised rs
			frame rs_`k'_`sexlabel' {
				
			cd "$MAIN_RESULT_DIR"
			save rs`k'_`site'_`sexlabel'.dta,replace
			}
			
			            // Proportion mediated via stage (age group): (TIE/TCE)*100
            capture frame drop pm_`k'_`sexlabel'
            frame copy tce_`k'_`sexlabel' pm_`k'_`sexlabel', replace

            frame pm_`k'_`sexlabel' {
                frlink 1:1 tt, frame(tie_`k'_`sexlabel')
                frget PE, from(tie_`k'_`sexlabel') prefix(tie_)
                frget tie_*, from(tie_`k'_`sexlabel')

                capture drop pm_PE pm_all_sd pm_all_lci pm_all_uci pm_*
                gen pm_PE = (tie_PE / PE) * 100

                forvalues i = 1/$m {
                    gen pm_`i' = (tie_`i' / tce_`i') * 100
                }

                egen pm_all_sd  = rowsd(pm_*)
                gen  pm_all_lci = pm_PE - `z'*pm_all_sd
                gen  pm_all_uci = pm_PE + `z'*pm_all_sd

                keep tt pm_PE pm_all_lci pm_all_uci
                replace pm_PE=0 if tt==0
                replace pm_all_lci=0 if tt==0
                replace pm_all_uci=0 if tt==0

                cd "$MAIN_RESULT_DIR"
                save pm`k'_`site'_`sexlabel'.dta, replace
            }

			
			// Total causal effect
			frame tce_`k'_`sexlabel' {
				capture drop tce_all_sd tce_all_lci tce_all_uci
				
				egen tce_all_sd = rowsd(tce_*)
				gen  tce_all_lci = PE - `z'*tce_all_sd
				gen  tce_all_uci = PE + `z'*tce_all_sd
				
				keep tt PE tce_all_lci tce_all_uci
				replace PE=0 if tt==0
				replace tce_all_lci=0 if tt==0
				replace tce_all_uci=0 if tt==0
				
				cd "$MAIN_RESULT_DIR"
				save tce`k'_`site'_`sexlabel'.dta,replace
			}
		
			// Direct effect
			frame tde_`k'_`sexlabel' {
				capture drop tde_all_sd tde_all_lci tde_all_uci
				
				egen tde_all_sd = rowsd(tde_*)
				gen  tde_all_lci = PE - `z'*tde_all_sd
				gen  tde_all_uci = PE + `z'*tde_all_sd
				
				keep tt PE tde_all_lci tde_all_uci
				replace PE=0 if tt==0
				replace tde_all_lci=0 if tt==0
				replace tde_all_uci=0 if tt==0
				
				cd "$MAIN_RESULT_DIR"
				save tde`k'_`site'_`sexlabel'.dta,replace
			}
		
			// Indirect effect
			frame tie_`k'_`sexlabel' {
				capture drop tie_all_sd tie_all_lci tie_all_uci
				
				egen tie_all_sd = rowsd(tie_*)
				gen  tie_all_lci = PE - `z'*tie_all_sd
				gen  tie_all_uci = PE + `z'*tie_all_sd
				
				keep tt PE tie_all_lci tie_all_uci
				replace PE=0 if tt==0
				replace tie_all_lci=0 if tt==0
				replace tie_all_uci=0 if tt==0
				
				cd "$MAIN_RESULT_DIR"
				save tie`k'_`site'_`sexlabel'.dta,replace
			}

			// Total avoidable deaths
			frame AD_`k'_`sexlabel' {
				capture drop AD_all_sd AD_all_lci AD_all_uci
				
				egen AD_all_sd = rowsd(AD_*)
				gen  AD_all_lci = PE - `z'*AD_all_sd
				gen  AD_all_uci = PE + `z'*AD_all_sd

				keep tt PE AD_all_lci AD_all_uci
				replace PE=0 if tt==0
				replace AD_all_lci=0 if tt==0
				replace AD_all_uci=0 if tt==0
				
				cd "$MAIN_RESULT_DIR"
				save AD`k'_`site'_`sexlabel'.dta,replace
			}

			// AD by removing relative survival differences
			frame ADa_`k'_`sexlabel' {
				capture drop ADa_all_sd ADa_all_lci ADa_all_uci
				
				egen ADa_all_sd = rowsd(ADa_*)
				gen  ADa_all_lci = PE - `z'*ADa_all_sd
				gen  ADa_all_uci = PE + `z'*ADa_all_sd

				keep tt PE ADa_all_lci ADa_all_uci
				replace PE=0 if tt==0
				replace ADa_all_lci=0 if tt==0
				replace ADa_all_uci=0 if tt==0
				
				cd "$MAIN_RESULT_DIR"
				save ADa`k'_`site'_`sexlabel'.dta,replace
			}

			// AD by removing stage differences
			frame ADb_`k'_`sexlabel' {
				capture drop ADb_all_sd ADb_all_lci ADb_all_uci
				
				egen ADb_all_sd = rowsd(ADb_*)
				gen  ADb_all_lci = PE - `z'*ADb_all_sd
				gen  ADb_all_uci = PE + `z'*ADb_all_sd

				keep tt PE ADb_all_lci ADb_all_uci
				replace PE=0 if tt==0
				replace ADb_all_lci=0 if tt==0
				replace ADb_all_uci=0 if tt==0
				
				cd "$MAIN_RESULT_DIR"
				save ADb`k'_`site'_`sexlabel'.dta,replace
			}
		}
	}
}
