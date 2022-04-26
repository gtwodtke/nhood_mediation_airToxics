#delimit cr
capture clear all 
capture log close 
set more off 

global data_directory "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\data\" 
global project_directory "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\" 
global log_directory "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\programs\_LOGS\" 
log using "${log_directory}06_create_v03_eclsb_merged.log", replace 

/******************************************************************
PROGRAM NAME: 06_create_v03_eclsb_merged.do
AUTHOR: KW, GW
PURPOSE: Merge ECLS-B with NCDB, RSEI, and CACES
*******************************************************************/

/*******************************
MERGE NCDB, RSEI, AND CACES DATA
********************************/
use "${data_directory}ncdb\v01_ncdb.dta"
keep zcta5 region nhdadvg* nhpovrt* nhpop*
sort zcta5

merge 1:m zcta5 using "${data_directory}caces\v01_caces.dta" 
tab _merge
keep if _merge==3
drop _merge 

merge 1:1 zcta5 using "${data_directory}rsei\v01_rsei.dta" 
tab _merge
keep if _merge==3
drop _merge

rename zcta5 zip01
rename region region01
keep zip01 *01 *_2001
sort zip01

save "${data_directory}_TEMP\ncdb_rsei_caces.dta", replace
clear

/****************
MERGE WITH ECLS-B
*****************/
use "${data_directory}eclsb\v02_eclsb_nvars.dta"

merge m:1 zip01 using "${data_directory}_TEMP\ncdb_rsei_caces.dta"
drop if _merge == 2 
drop _merge

/**********
RECODE VARS
***********/
/*TRUNCATE EXTREME OUTLIERS FOR RSEI TOXICS*/
foreach i of var chem* {
	quietly centile `i', c(99) 
	quietly replace `i'=. if `i'>r(c_1) & `i'<9999999
	}

/*DROP RSEI TOXICS WITH <=1% NONZERO EXPOSURE*/
drop chem189_2001 chem233_2001 chem541_2001 chem557_2001

/**********
SAVE MERGED
***********/
saveold "${data_directory}\eclsb\v03_eclsb_merged.dta", replace v(12)

erase "${data_directory}_TEMP\ncdb_rsei_caces.dta"

/***********
DESCRIPTIVES
************/
codebook 
centile nhpovrt01, c(10 20 30 40 50 60 70 80 90) 

capture clear all

log close