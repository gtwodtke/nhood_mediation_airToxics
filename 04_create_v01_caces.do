#delimit cr
capture clear all 
capture log close 
set more off 

global data_directory "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\data\caces\" 
global savedata_directory "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\data\caces\"
global crosswalk_directory "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\data\crosswalks\" 
global log_directory "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\programs\_LOGS\" 
log using "${log_directory}04_create_v01_caces.log", replace 

/******************************************************************
PROGRAM NAME: 04_create_v01_caces.do
AUTHOR: KW, GW
PURPOSE: input, clean, and aggregate data from CACES
*******************************************************************/

/*******************
INPUT RAW CACES DATA
********************/
import delimited "${data_directory}caces_raw.csv", encoding(ISO-8859-1)

format fips %15.0g
sort fips year pollutant

rename pred_wght poll_ 
rename fips geo2010

label var geo2010 "Concatenation of state, county and tract for 2010 Census Tracts"

keep geo2010 pollutant year poll_ 
order geo2010 year pollutant poll_  

sort geo2010 year pollutant

reshape wide poll_, i(geo2010 year) j(pollutant) string

rename poll_co co_
rename poll_no2 no2_
rename poll_o3 o3_
rename poll_pm10 pm10_
rename poll_pm25 pm25_
rename poll_so2 so2_

reshape wide co_ no2_ o3_ pm10_ pm25_ so2_, i(geo2010) j(year)

sum co_1980 - so2_2015

save "${savedata_directory}cacesv1.dta", replace 
clear 

/******************************
2010 TO 2000 CENSUS TRACT XWALK
*******************************/
insheet using "${crosswalk_directory}us2010trf.txt", clear delim(",")

keep v4 v13 v22 v24

rename v4 geo2000
rename v13 geo2010
rename v22 landpct00
rename v24 landpct10

label var geo2000 "Concatenation of 2000 state, county and tract for Census 2000 Tract"
label var geo2010 "Concatenation of state, county and tract for 2010 Census Tracts"
label var landpct00 "Percent of the total 2000 tract land area represents"
label var landpct10 "Percent of the total 2010 tract land area represents"

format geo2000 %15.0g
format geo2010 %15.0g

sort geo2010 geo2000

by geo2010: gen freq=_n
by geo2010: egen maxfreq=max(freq)

gen xwgt=landpct10/100 if maxfreq>1
replace xwgt=landpct00/100 if maxfreq==1

egen sumwt = sum(xwgt), by(geo2010)
sum sumwt, det
drop sumwt

sort geo2000 geo2010
egen sumwt = sum(xwgt), by(geo2000)
replace xwgt = xwgt/sumwt
drop sumwt

egen sumwt = sum(xwgt), by(geo2000)
sum sumwt, det
drop sumwt

keep geo* xwgt

sort geo2010 geo2000
merge m:1 geo2010 using "${savedata_directory}cacesv1.dta"
tab _merge
drop if _merge == 1 
drop _merge

sort geo2010 geo2000
foreach var of varlist ///
		co_1980-so2_2015 {
			replace `var'= `var'*xwgt
			} 

collapse (sum) co_1980-so2_2015, by(geo2000)
sort geo2000
sum co_1980-so2_2015

save "${savedata_directory}cacesv1.dta", replace 
clear 

/***************************************
2000 CENSUS TRACT TO 2000 ZIP CODE XWALK
****************************************/
use "${crosswalk_directory}crosswalk_lndareawt.dta"

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

egen sumwt = sum(afact), by(zcta5)
replace afact = afact/sumwt
drop sumwt

egen sumwt = sum(afact), by(zcta5)
sum sumwt, det
drop sumwt

merge m:1 geo2000 using "${savedata_directory}cacesv1.dta"
tab _merge
keep if _merge == 3 
drop _merge 

sort zcta5 geo2000

foreach var of varlist ///
	co_1980-so2_2015 {
		replace `var'= `var'*afact
		} 

collapse (sum) co_1980-so2_2015, by(zcta5)
sum co_1980-so2_2015

/**********
CLEAN ZCTAs
***********/
drop if regexm(zcta5, "HH") 
drop if regexm(zcta5, "XX") 
destring zcta5, replace

/*********
SAVE CACES
**********/
save "${savedata_directory}v01_caces.dta", replace 

erase "${savedata_directory}cacesv1.dta"

/***********
DESCRIPTIVES
************/
*codebook 

clear 

log close











