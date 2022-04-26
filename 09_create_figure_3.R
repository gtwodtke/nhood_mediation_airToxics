################################################
################################################
##                                            ##
## PROGRAM NAME: 09_create_figure_3.R         ##
## AUTHOR: GW                                 ##
## DESCRIPTION:                               ##
##                                            ##
##  create scatterplot matrix of nh poverty   ##
##  and selected toxics                       ##
##                                            ##
################################################
################################################

##### LOAD LIBRARIES #####
rm(list=ls())
library(haven)
library(foreign)
library(dplyr)
library(tidyr)
library(ggplot2)
library(psych)

startTime<-Sys.time()

##### LOAD ECLS-B #####
eclsb<-read.dta("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\eclsb\\v04_eclsb_mi.dta")
eclsb<-as.data.frame(eclsb)
eclsb<-eclsb[order(eclsb$caseid,eclsb$minum),]

##### DEFINE VARIABLE SETS #####
vars.chem01<-c(
	"chem6_2001",
	"chem10_2001",
	"chem12_2001",
	"chem17_2001",
	"chem36_2001",
	"chem39_2001",
	"chem40_2001",
	"chem70_2001",
	"chem75_2001",
	"chem108_2001",
	"chem109_2001",
	"chem114_2001",
	"chem141_2001",
	"chem163_2001",
	"chem180_2001",
	"chem217_2001",
	"chem225_2001",
	"chem236_2001",
	"chem290_2001",
	"chem293_2001",
	"chem323_2001",
	"chem332_2001",
	"chem346_2001",
	"chem347_2001",
	"chem351_2001",
	"chem355_2001",
	"chem356_2001",
	"chem359_2001",
	"chem360_2001",
	"chem362_2001",
	"chem364_2001",
	"chem381_2001",
	"chem447_2001",
	"chem454_2001",
	"chem474_2001",
	"chem519_2001",
	"chem522_2001",
	"chem531_2001",
	"chem549_2001",
	"chem565_2001",
	"chem566_2001",
	"chem567_2001",
	"chem572_2001",
	"chem586_2001",
	"chem592_2001",
	"chem599_2001") 

vars.mtx<-c("nhpovrt01","co_2001","pm10_2001","chem364_2001","chem359_2001","chem355_2001")

eclsb<-eclsb[which(eclsb$minum==1),c(vars.mtx)]

eclsb<-rename(eclsb,
	Nh.Poverty=nhpovrt01,
	CO=co_2001,
	PM10=pm10_2001,
	Methanol=chem364_2001,
	Hg=chem359_2001,
	Mn=chem355_2001)

##### CREATE SCATTERPLOT MATRIX #####
tiff("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\figures\\figure_3.tiff",
	width=7,
	height=7,
	units='in',
	res=600)

pairs.panels(eclsb,
	smooth=TRUE,
	ci=FALSE,
	scale=FALSE,
	density=FALSE,
	ellipses=FALSE,
	method="pearson",
	pch=1,
	lm=FALSE,
	cor=TRUE,
	stars=FALSE,
	rug=TRUE,
	hist.col="gray85",
	col="gray65")

dev.off()

##### PRINT OUTPUT #####
sink("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\programs\\_LOGS\\09_create_figure_3_log.txt")

print(startTime)
print(Sys.time())

sink()

