#delimit cr
capture clear all 
capture log close 
set more off 

global data_directory "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\data\rsei\raw\" 
global savedata_directory "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\data\rsei\"
global crosswalk_directory "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\data\crosswalks\" 
global log_directory "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\programs\_LOGS\" 
log using "${log_directory}05_create_v01_rsei.log", replace 

/******************************************************************
PROGRAM NAME: 05_create_v01_rsei.do
AUTHOR: KW, GW
PURPOSE: input, clean, and aggregate data from RSEI-GM
*******************************************************************/

/*******************
INPUT RAW RSEI DATA
********************/
use "${data_directory}2001neuro_2000tracts.dta"
gen year=.
replace year = 2001

/*********
CLEAN DATA
**********/
rename *, lower
rename fips geo2000
replace geo2000 = " " if regexm(geo2000, "nodata")
gen g = substr(geo2000,1,2)
order g
gen h = substr(geo2000,4,3)
order g h 
gen i = substr(geo2000,8,6)
order g h i
gen geo = g + h + i
order geo
drop geo2000 g h i
rename geo geo2000

destring geo2000, replace
format geo2000 %13.0f

drop *tc
drop *score

rename (*c)(*)
rename (*)(*_)
rename geo2000_ geo2000
rename year_ year
order geo2000 year
sort geo2000 year

unab vlist: chem*
di "`vlist'"
local j "year"
di "`j'"

levelsof `j', local(J) //
foreach var of varlist `vlist' {
	foreach j of local J {
		local newlist `newlist' `var'`j'
		local lablist "`lablist' `"`:variable label `var'' (`j')"'"
	}
}

reshape wide chem*, i(geo2000) j(year)

foreach new of local newlist {
	gettoken lab lablist : lablist
	quietly label var `new' "`lab'"
}

save "${savedata_directory}rseiv1.dta", replace 

/***************************************
2000 CENSUS TRACT TO 2000 ZIP CODE XWALK
****************************************/
use "${crosswalk_directory}crosswalk_lndareawt.dta", clear 

gen t = substr(tract,1,4)
gen u = substr(tract,6,2)
drop tract
gen tract = t + u 
gen geo2000 = county + tract
destring geo2000, replace
format geo2000 %15.0g

destring afact, replace 

keep geo2000 zcta5 afact
sort zcta5 geo2000

label var geo2000 "Concatenation of 2000 state, county and tract for Census 2000 Tract"
label var zcta5 "2000 ZIP Codes"
label var afact "Portion of 2000 Censust tract land area located within 2000 ZIP/ZCTA"

egen sumafact = sum(afact), by(geo2000)
sum sumafact, det
drop sumafact

merge m:1 geo2000 using "${savedata_directory}rseiv1.dta"
tab _merge
keep if _merge==3
drop _merge

sort zcta5 geo2000

egen sumwt = sum(afact), by(zcta5)
sum sumwt, det
replace afact=afact/sumwt
drop sumwt

foreach var of varlist chem6_2001-chem599_2001 {
	quietly replace `var'= `var'*afact
	} 

foreach v of var * {
	local l`v' : variable label `v'
	if `"`l`v''"' == "" {
		local l`v' "`v'"
		}
	}
		
collapse (sum) chem6_2001-chem599_2001, by(zcta5)

foreach v of var * {
	quietly label var `v' "`l`v''"
	}
	
sort zcta5

drop if regexm(zcta5, "HH") 
drop if regexm(zcta5, "XX") 
destring zcta5, replace

/********************
RECODE CONCENTRATIONS
*********************/
/*CONVERT TO NANOGRAMS PER M^3*/
foreach i of var chem* {
	quietly replace `i'=`i'*1000
	}

/*DROP O3 (REDUNDANT WITH CACES)*/
drop chem441*

/********
SAVE RSEI
*********/
save "${savedata_directory}v01_rsei.dta", replace 

erase "${savedata_directory}rseiv1.dta"

/***********
DESCRIPTIVES
************/
codebook 

clear 

log close

