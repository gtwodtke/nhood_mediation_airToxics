################################################
################################################
##                                            ##
## PROGRAM NAME: 07_create_v04_eclsb_mi.R     ##
## AUTHOR: GW                                 ##
## DESCRIPTION:                               ##
##                                            ##
##  multiply impute missing data              ##
##                                            ##
################################################
################################################

##### LOAD LIBRARIES #####
rm(list=ls())
library(haven)
library(foreign)
library(dplyr)
library(tidyr)
library(Hmisc)
library(ggplot2)
library(mice)
library(RfEmpImp)
library(randomForest)
library(foreach)
library(doParallel)
library(doRNG)

startTime<-Sys.time()

set.seed(8675309)

##### DEFINE CONTROL PARAMETERS #####
n.cores<-parallel::detectCores()-1
nmi<-5

##### LOAD ECLS-B #####
eclsb<-read.dta("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\eclsb\\v03_eclsb_merged.dta")
eclsb<-eclsb[order(eclsb$caseid),]

chem.labels<-attr(eclsb,"var.labels")
chem.key<-data.frame(var.name=names(eclsb),chem.labels)

##### DEFINE VARIABLE SETS #####
vars.meta<-c(
	"caseid",
	"famid",
	"strat",
	"psu",
	"sampwt",
	"zip01")

vars.base<-c(
	"gender",
	"race",
	"twinid",
	"birthwt",
	"momage01",
	"dadage01",
	"wic01",
	"foodst01",
	"medicd01",
	"tanf01")

vars.y<-c(
	"maththeta05",
	"readtheta05",
	"mentalsc01",
	"mentalsc03",
	"motorsc01",
	"motorsc03",
	"payattn01",
	"payattn03")

vars.cov01<-c(
	"nhpovrt01",
	"faminc01",
	"pared01",
	"parocc01",
	"momemp01",
	"dademp01",
	"hhtotal01",
	"biodad01",
	"married01",
	"house01",
	"rbooks01",
	"prmlang01",
	"nhpopden01",
	"urban01",
	"region01",
	"age01",
	"age05")

vars.chem01<-c(
	"co_2001",
	"no2_2001",
	"o3_2001",
	"pm10_2001",
	"pm25_2001",
	"so2_2001",
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

eclsb<-as.data.frame(eclsb[,c(vars.meta,vars.base,vars.y,vars.cov01,vars.chem01)])
eclsb.mi.vars<-eclsb[,-which(names(eclsb) %in% c(vars.meta))]

##### IMPUTE MISSING DATA #####
eclsb.mi<-imp.rfnode.cond(eclsb.mi.vars,num.imp=nmi,seed=8675309,num.threads=n.cores)

eclsb$minum<-0
for (i in 1:nmi) {
	mi<-complete(eclsb.mi,action=i)
	mi_join<-data.frame(eclsb[which(eclsb$minum==0),vars.meta],mi)
	mi_join$minum<-i
	eclsb<-rbind(eclsb,mi_join)
	}

table(eclsb$minum)
eclsb$miwt<-1/nmi
eclsb$fnlwt<-eclsb$miwt*eclsb$sampwt

##### SAVE IMPUTED DATA #####
save(eclsb,file="C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\eclsb\\v04_eclsb_mi.RData")
write.dta(eclsb,"C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\eclsb\\v04_eclsb_mi.dta",version=10)

##### PRINT OUTPUT #####
sink("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\programs\\_LOGS\\07_create_v04_eclsb_mi_log.txt")

summary(eclsb[which(eclsb$minum==0),])
summary(eclsb[which(eclsb$minum!=0),])

print(startTime)
print(Sys.time())

sink()

sink("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\programs\\_LOGS\\_RSEI_chem_labels.txt")

chem.key<-chem.key[chem.key$var.name %in% vars.chem01,]
print(chem.key)

sink()


