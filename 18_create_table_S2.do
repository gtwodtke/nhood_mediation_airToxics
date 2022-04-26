#delimit cr
capture clear all 
capture log close 
set more off 

global projdir "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\" 
log using "${projdir}programs\_LOGS\18_create_table_S2.log", replace 

/******************************************************************
PROGRAM NAME: 18_create_table_S2.do
AUTHOR: GW, BP
PURPOSE: create table of sample descriptives
*******************************************************************/

/********
LOAD DATA
*********/
use "${projdir}data\eclsb\v04_eclsb_mi.dta"
keep if minum!=0

/**********************
CREATE OUTPUT FOR TABLE
***********************/
/*COVARIATES*/
tab1 gender race twinid region01 urban01 married01 biodad01 house01 ///
	prmlang01 rbooks01 dademp01 momemp01 pared01 ///
	wic01 foodst01 medicd01 tanf01 [iw=fnlwt]
sum birthwt age01 age05 momage01 dadage01 faminc01 parocc01 hhtotal01 ///
	nhpovrt01 maththeta05 readtheta05 nhpopden01 [iw=fnlwt]

log close
