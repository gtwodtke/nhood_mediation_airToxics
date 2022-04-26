#delimit cr
capture clear all 
capture log close 
set more off 

global data_directory "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\data\ncdb\csvFiles\" 
global savedata_directory "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\data\ncdb\"
global log_directory "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\programs\_LOGS\" 
global crosswalk_directory "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\data\crosswalks\"
log using "${log_directory}03_create_v01_ncdb.log", replace 

/******************************************************************
PROGRAM NAME: 03_create_v01_ncdb.do
AUTHOR: KW, GW
PURPOSE: input, clean, interpolate, and aggregate data from NCDB
*******************************************************************/

/******************
INPUT RAW NCDB DATA
*******************/
import delimited "${data_directory}ncdb_raw.csv", encoding(ISO-8859-1) 

/**********
RENAME VARS
***********/
rename trctpop0 trctpop_2000 
rename trctpop1 trctpop_2010 
rename educ80 	educ8_2000 
rename educ110 	educ11_2000 
rename educ120 	educ12_2000 
rename educ150 	educ15_2000 
rename educa0 	educa_2000 
rename educ160 	educ16_2000 
rename educpp0	educpp_2000 
rename educ81a 	educ8_2010 
rename educ111a educ11_2010 
rename educ121a educ12_2010 
rename educ151a educ15_2010 
rename educa1a 	educa_2010 
rename educ161a	educ16_2010 
rename educpp1a	educpp_2010 
rename	ffh0d	ffhd_2000 
rename	ffh0n	ffhn_2000 
rename	ffh1ad	ffhd_2010 
rename	ffh1an	ffhn_2010 
rename povrat0d		povratd_2000 
rename povrat0n		povratn_2000 
rename povrat1ad	povratd_2010 
rename povrat1an	povratn_2010 
rename unempt0d		unemptd_2000
rename unempt0n		unemptn_2000
rename unempt1an	unemptn_2010
rename unempt1ad	unemptd_2010
rename shrwht0n	shrwhtn_2000 
rename shrwht1n	shrwhtn_2010 
rename shrblk0n shrblkn_2000 
rename shrblk1n shrblkn_2010 

keep geo2010 state region division county trctpop* educ* ffhd* ffhn* ///
povratd* povratn* unemptd* unemptn* shrwhtn* shrblkn* 

sort geo2010
save "${savedata_directory}ncdbv1.dta", replace 
clear 

/******************************
2010 TO 2000 CENSUS TRACT XWALK
*******************************/
insheet using "${crosswalk_directory}us2010trf.txt", clear delim(",")

keep v4 v13 v29 v30

rename v4 geo2000
rename v13 geo2010
rename v29 hupct00
rename v30 hupct10

label var geo2000 "Concatenation of 2000 state, county and tract for Census 2000 Tract"
label var geo2010 "Concatenation of state, county and tract for 2010 Census Tracts"
label var hupct00 "Percent of the HU00 this record contains"
label var hupct10 "Percent of the HU10 this record contains"

format geo2000 %15.0g
format geo2010 %15.0g

sort geo2010 geo2000
by geo2010: gen freq=_n
by geo2010: egen maxfreq=max(freq)

gen xwgt=hupct10/100 if maxfreq>1
replace xwgt=hupct00/100 if maxfreq==1

egen sumwt = sum(xwgt), by(geo2010)
sum sumwt, det
drop sumwt

keep geo* xwgt

merge m:1 geo2010 using "${savedata_directory}ncdbv1.dta"
drop if _merge == 1
drop _merge

sort geo2010 geo2000
foreach var of varlist ///
	trctpop_2000-povratd_2010 {
		replace `var'= `var'*xwgt
		} 

collapse (sum) trctpop_2000-povratd_2010 (firstnm) state division county, by(geo2000)
sort geo2000

save "${savedata_directory}ncdbv1.dta", replace 
clear

/***************************************
2000 CENSUS TRACT TO 2000 ZIP CODE XWALK
****************************************/
use "${crosswalk_directory}crosswalk_popsizewt.dta"

gen t = substr(tract,1,4)
gen u = substr(tract,6,2)
drop tract
gen tract = t + u 
gen geo2000 = county + tract
destring geo2000, replace
format geo2000 %15.0g

destring afact, replace 
destring pop2k, replace 

gen zipx = strpos(zcta5, "X") > 0 
tabulate zipx 
drop if zipx == 1
destring zcta5, replace

keep geo2000 zcta5 afact
sort geo2000 zcta5

label var geo2000 "Concatenation of 2000 state, county and tract for Census 2000 Tract"
label var zcta5 "2000 ZIP Codes"
label var afact "Portion of 2000 Censust tract population located within 2000 ZIP/ZCTA"

egen sumafact = sum(afact), by(geo2000)
sum sumafact, det
drop sumafact

merge m:1 geo2000 using "${savedata_directory}ncdbv1.dta"
tab _merge
drop if _merge == 2  

foreach var of varlist ///
	trctpop_2000-povratd_2010 {
		replace `var'= `var'*afact
		} 

collapse (sum) trctpop_2000-povratd_2010 (firstnm) state division county, by(zcta5)
sort zcta5

save "${savedata_directory}ncdbv1.dta", replace 

/********
LAND AREA
*********/
use "${crosswalk_directory}crosswalk_lndareawt.dta", clear 

gen t = substr(tract,1,4)
gen u = substr(tract,6,2)
drop tract
gen tract = t + u 
gen geo2000 = county + tract
destring geo2000, replace
format geo2000 %15.0g

destring afact landsqmi, replace 

keep geo2000 zcta5 landsqmi 
sort geo2000 zcta5

collapse (sum) landsqmi, by(zcta5)

replace zcta5 = " " if regexm(zcta5, "HH")
replace zcta5 = " " if regexm(zcta5, "XX")
destring zcta5, replace
drop if zcta5==.

merge 1:1 zcta5 using "${savedata_directory}ncdbv1.dta"
keep if _merge==3
drop _merge

/*****************************************
LINEAR INTERPOLATION FOR INTERCENSAL YEARS
******************************************/
capture macro drop vars08 
global vars08 educ8 educ11 educ12 educ15 educa educ16 educpp ///
			  ffhd ffhn povratd povratn unemptd unemptn ///
			  trctpop shrwhtn shrblkn

foreach v of global vars08 { 
	forval t=2001/2009 { 
		gen `v'_`t'=`v'_2000+(`v'_2010-`v'_2000)*((`t'-2000)/(2010-2000)) 
		} 
	} 

/***********
RESHAPE DATA
************/
capture macro drop stubs 
global stubs  educ8_ educ11_ educ12_ educ15_ educa_ educ16_ educpp_ ///
			  ffhd_ ffhn_ povratd_ povratn_ unemptd_ unemptn_ ///
			  trctpop_ shrwhtn_ shrblkn_ 

reshape long $stubs, i(zcta5) j(year) 

foreach v in educ8 educ11 educ12 educ15 educa educ16 educpp ///
			 ffhd ffhn povratd povratn unemptd unemptn ///
			 trctpop shrwhtn shrblkn { 
			 rename `v'_ `v' 
			 } 

/*******************
CREATE NEW VARIABLES
********************/
/***EDUCATIONAL COMPOSITION***/
gen nhlesshs=(educ8+educ11)/educpp 
gen nhhsgrad=(educ12)/educpp 
gen nhsomcol=(educ15+educa)/educpp 
gen nhcolgrd=(educ16)/educpp 
foreach v in nhlesshs nhhsgrad nhsomcol nhcolgrd { 
	replace `v'=. if inrange(`v',1,99) 
	} 
sum nhlesshs-nhcolgrd 

/***FEMALE-HEADED FAMILIES WITH CHILDREN***/
gen nhfemhd=ffhn/ffhd 
replace nhfemhd=. if nhfemhd>1 
sum nhfemhd 

/***POVERTY RATE***/
gen nhpovrt=povratn/povratd 
sum nhpovrt 

/***UNEMPLOYMENT RATE*/
gen nhunemprt=unemptn/unemptd
sum nhunemprt 

/***RACIAL COMPOSITION***/
gen nhshrwht=1-(shrwhtn/trctpop )
replace nhshrwht=. if nhshrwht>1 
sum nhshrwht 

gen nhshrblk=shrblkn/trctpop 
replace nhshrblk=. if nhshrblk>1 
sum nhshrblk 

/***COMPOSITE NH DADVG SCALE***/
pca nhpovrt nhunemprt nhfemhd nhlesshs nhshrwht
predict nhdadvg 

/***POPULATION DENSITY***/
gen nhpopden=trctpop/landsqmi
sum nhpopden, detail 

/***CENSUS REGION***/
recode division (1 2 = 1) (3 4 = 2) (5 6 7 = 3) (8 9 = 4), gen(region)

label def region_lbl ///
	1 "Northeast" ///
	2 "Midwest" ///
	3 "South" ///
	4 "West"
	
label values region region_lbl

/***********
RESHAPE DATA
************/
sort zcta5 year
keep zcta5 state region division county year n* 

capture macro drop stubs 
global stubs nhlesshs nhhsgrad nhsomcol nhcolgrd nhfemhd nhpovrt ///
nhunemprt nhshrwht nhshrblk nhpopden nhdadvg 
reshape wide $stubs, i(zcta5) j(year) 
foreach v of global stubs { 
	rename `v'2000 `v'00 
	rename `v'2001 `v'01 
	rename `v'2002 `v'02 
	rename `v'2003 `v'03 
	rename `v'2004 `v'04 
	rename `v'2005 `v'05 
	rename `v'2006 `v'06 
	rename `v'2007 `v'07 
	rename `v'2008 `v'08 
	rename `v'2009 `v'09 
	rename `v'2010 `v'10 
	} 

save "${savedata_directory}ncdbv1.dta", replace

/********
SAVE NCDB
*********/
save "${savedata_directory}v01_ncdb.dta", replace 

erase "${savedata_directory}ncdbv1.dta"

/***********
DESCRIPTIVES
************/
*codebook 
centile nhpovrt01, c(10 20 30 40 50 60 70 80 90) 
centile nhdadvg01, c(10 20 30 40 50 60 70 80 90) 

clear

log close 


