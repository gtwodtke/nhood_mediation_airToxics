################################################
################################################
##                                            ##
## PROGRAM NAME: 15_create_figure_8           ##
## AUTHOR: GW                                 ##
## DESCRIPTION:                               ##
##                                            ##
##  create partial dependence plots           ##
##                                            ##
################################################
################################################

rm(list=ls())

list.of.packages<-c(
	"haven",
	"foreign",
	"dplyr",
	"tidyr",
	"ggplot2",
	"ranger",
	"caret",
	"e1071",
	"SuperLearner",
	"foreach",
	"doParallel",
	"doRNG",
	"gridExtra",
	"Hmisc")

for(package.i in list.of.packages) {
	suppressPackageStartupMessages(library(package.i, character.only = TRUE))
	}

startTime<-Sys.time()

set.seed(8675309)

##### LOAD ECLS-B #####
eclsb.mi<-read.dta("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\eclsb\\v04_eclsb_mi.dta")
eclsb.mi<-as.data.frame(eclsb.mi[which(eclsb.mi$minum!=0),])
eclsb.mi<-eclsb.mi[order(eclsb.mi$caseid,eclsb.mi$minum),]

nmi<-5
ntrees<-200
n.points<-50
n.cores<-parallel::detectCores()-1

##### DEFINE VARIABLE SETS #####
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

vars.covw1<-c(
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

vars.ntxw1<-c(
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

chem.label<-c(
	"CO",
	"NO2",
	"O3",
	"PM10",
	"PM2.5",
	"SO2",
	"Acetophenone",
	"Acrylamide",
	"Acrylonitrile",
	"Allyl.Chloride",
	"Anthracene",
	"As",
	"As.Compounds",
	"Br",
	"Methyl.Bromide",
	"Cd",
	"Cd.Compounds",
	"CS2",
	"Methyl.Chloride",
	"Cresol",
	"C12Br10O",
	"Dichlorobenzene",
	"Dichloromethane",
	"C10H12",
	"Ethylbenzene",
	"Ethylene.Oxide",
	"C2Cl6",
	"H.Cyanide",
	"Pb",
	"Pb.Compounds",
	"Malathion",
	"Mn",
	"Mn.Compounds",
	"Hg",
	"Hg.Compounds",
	"C4H5N",
	"Methanol",
	"C5H12O",
	"C2HCl5",
	"Phenol",
	"PCBs",
	"Styrene",
	"SO2F2",
	"C2Cl4",
	"Toluene",
	"C2H3Cl3.111",
	"C2H3Cl3.112",
	"C2HCl3",
	"Triethylamine",
	"Vinyl.Chloride",
	"Xylene",
	"n.Hexane")

##### STANDARDIZE VARS #####
eclsb.mi$readtheta05<-(eclsb.mi$readtheta05-weighted.mean(eclsb.mi$readtheta05,eclsb.mi$sampwt))/sqrt(wtd.var(eclsb.mi$readtheta05,eclsb.mi$sampwt))
eclsb.mi$maththeta05<-(eclsb.mi$maththeta05-weighted.mean(eclsb.mi$maththeta05,eclsb.mi$sampwt))/sqrt(wtd.var(eclsb.mi$maththeta05,eclsb.mi$sampwt))

##### ESTIMATE ALEs #####
x.pm10<-seq(from=quantile(eclsb.mi$pm10_2001,p=0.0025),to=quantile(eclsb.mi$pm10_2001,p=0.9975),length=n.points)
x.o3<-seq(from=quantile(eclsb.mi$o3_2001,p=0.0025),to=quantile(eclsb.mi$o3_2001,p=0.9975),length=n.points)
x.nhpov<-seq(from=quantile(eclsb.mi$nhpovrt01,p=0.0025),to=quantile(eclsb.mi$nhpovrt01,p=0.9975),length=n.points)

miest.rd.o3<-miest.rd.pm10<-miest.mt.o3<-miest.mt.pm10<-matrix(data=NA,nrow=n.points,ncol=nmi)
miest.o3.nh<-miest.pm10.nh<-matrix(data=NA,nrow=n.points,ncol=nmi)

for (i in 1:nmi) {

	### LOAD MI DATA ###
	print(c("nmi=",i))
	eclsb<-eclsb.mi[which(eclsb.mi$minum==i),]

	### TUNE HYPERPARAMETERS ###
	registerDoSEQ()

	cntrl<-trainControl(method="CV",number=5)

	rf.m2.grid<-expand.grid(
		min.node.size=c(5,10,15,20),
		mtry=floor(length(c(vars.base,vars.covw1,vars.ntxw1))*c(0.3,0.4,0.5,0.6,0.7)),
		splitrule="variance")

	rf.rd.m2<-train(readtheta05~.,
		data=eclsb[,c("readtheta05",vars.base,vars.covw1,vars.ntxw1)],
		method="ranger",
		tuneGrid=rf.m2.grid,
		trControl=cntrl,
		metric="RMSE",
		respect.unordered.factors=TRUE,
		num.trees=ntrees,
		weights=eclsb$sampwt,
		seed=8675309)

	rf.mt.m2<-train(maththeta05~.,
		data=eclsb[,c("maththeta05",vars.base,vars.covw1,vars.ntxw1)],
		method="ranger",
		tuneGrid=rf.m2.grid,
		trControl=cntrl,
		metric="RMSE",
		respect.unordered.factors=TRUE,
		num.trees=ntrees,
		weights=eclsb$sampwt,
		seed=8675309)

	### COMPUTE M->Y ALEs ###
	read.train<-eclsb[,"readtheta05"]
	math.train<-eclsb[,"maththeta05"]
	x.train<-eclsb[,c(vars.base,vars.covw1,vars.ntxw1)]
	wts.train<-eclsb[,"sampwt"]
	
	m1.rd<-ranger(
		y=read.train,
		x=x.train,
		num.trees=ntrees,
		min.node.size=rf.rd.m2$bestTune[,"min.node.size"],
		mtry=rf.rd.m2$bestTune[,"mtry"],
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=wts.train,
		seed=8675309)

	m1.mt<-ranger(
		y=math.train,
		x=x.train,num.trees=ntrees,
		min.node.size=rf.mt.m2$bestTune[,"min.node.size"],
		mtry=rf.mt.m2$bestTune[,"mtry"],
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=wts.train,
		seed=8675309)

	for (j in 1:n.points) {
		gcomp<-x.train
		gcomp$o3_2001<-x.o3[j]
		yhat.rd<-predict(m1.rd,gcomp)$pred
		yhat.mt<-predict(m1.mt,gcomp)$pred
		miest.rd.o3[j,i]<-weighted.mean(yhat.rd,wts.train)
		miest.mt.o3[j,i]<-weighted.mean(yhat.mt,wts.train)
		}

	for (j in 1:n.points) {
		gcomp<-x.train
		gcomp$pm10_2001<-x.pm10[j]
		yhat.rd<-predict(m1.rd,gcomp)$pred
		yhat.mt<-predict(m1.mt,gcomp)$pred
		miest.rd.pm10[j,i]<-weighted.mean(yhat.rd,wts.train)
		miest.mt.pm10[j,i]<-weighted.mean(yhat.mt,wts.train)
		}

	### COMPUTE A->M ALEs ###
	o3.train<-eclsb[,"o3_2001"]
	pm10.train<-eclsb[,"pm10_2001"]
	c.train<-eclsb[,c(vars.base,vars.covw1)]
	
	m1.o3<-ranger(
		y=o3.train,
		x=c.train,
		num.trees=ntrees,
		min.node.size=10,
		mtry=round(length(c(vars.base,vars.covw1))/3),
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=wts.train,
		seed=8675309)

	m1.pm10<-ranger(
		y=pm10.train,
		x=c.train,
		num.trees=ntrees,
		min.node.size=10,
		mtry=round(length(c(vars.base,vars.covw1))/3),
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=wts.train,
		seed=8675309)

	for (j in 1:n.points) {
		gcomp<-c.train
		gcomp$nhpovrt01<-x.nhpov[j]
		mhat.o3<-predict(m1.o3,gcomp)$pred
		mhat.pm10<-predict(m1.pm10,gcomp)$pred
		miest.o3.nh[j,i]<-weighted.mean(mhat.o3,wts.train)
		miest.pm10.nh[j,i]<-weighted.mean(mhat.pm10,wts.train)
		}
	}

### COMBINE MI ESTIMATES ###
est.o3<-est.pm10<-est.nh<-matrix(data=NA,nrow=n.points,ncol=2)

for (j in 1:n.points) { 

	est.o3[j,1]<-mean(miest.rd.o3[j,])
	est.o3[j,2]<-mean(miest.mt.o3[j,])

	est.pm10[j,1]<-mean(miest.rd.pm10[j,])
	est.pm10[j,2]<-mean(miest.mt.pm10[j,])

	est.nh[j,1]<-mean(miest.o3.nh[j,])
	est.nh[j,2]<-mean(miest.pm10.nh[j,])
	}

### PRINT RESULTS ###
sink("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\programs\\_LOGS\\15_create_figure_8_log.txt")

output.o3<-as.data.frame(cbind(x.o3,est.o3))
output.pm10<-as.data.frame(cbind(x.pm10,est.pm10))
output.nh<-as.data.frame(cbind(x.nhpov,est.nh))

colnames(output.o3)<-colnames(output.pm10)<-c('xval','rd.score','mt.score')
colnames(output.nh)<-c('xval','o3','pm10')

print(output.o3)
print(output.pm10)
print(output.nh)

##### PLOT RESULTS #####
tiff("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\figures\\figure_8.tiff",
	height=9,
	width=9,
	units='in',
	res=600)

par(mfrow=c(2,2))

plot(output.nh$xval,output.nh$pm10,type="l",
	main="A. Dose-response Relationship \n of Neighborhood Poverty to PM10",
	xlab="Neighborhood Poverty Rate",
	ylab=bquote("PM10" ~ (ug/m^3)),
	ylim=c(20,28),
	xlim=c(0.0,0.5),
	cex.axis=0.9,
	cex.lab=0.88,
	cex.main=0.97)
rug(eclsb$nhpovrt01)

plot(output.nh$xval,output.nh$o3,type="l",
	main="B. Dose-response Relationship \n of Neighborhood Poverty to O3",
	xlab="Neighborhood Poverty Rate",
	ylab=bquote("O3" ~ (ppb)),
	ylim=c(47,52),
	xlim=c(0.0,0.5),
	cex.axis=0.9,
	cex.lab=0.88,
	cex.main=0.97)
rug(eclsb$nhpovrt01)

plot(output.pm10$xval,output.pm10$rd.score,type="l",lty="dotted",lwd=1.25,
	main="C. Dose-response Relationship \n of PM10 to Test Scores",
	xlab=bquote("PM10" ~ (ug/m^3)),
	ylab="Test Scores (SD)",
	ylim=c(-0.10,0.10),
	xlim=c(5,50),
	cex.axis=0.9,
	cex.lab=0.90,
	cex.main=0.97)
rug(eclsb$pm10_2001)
lines(output.pm10$xval,output.pm10$mt.score,lty="longdash")
labels<-c("Reading Scores","Math Scores")
type<-c("dotted","longdash")
legend("topright",inset=0.03,labels,lty=c("dotted","longdash"),lwd=c(1.25,1),seg.len=2,cex=0.95)

plot(output.o3$xval,output.o3$rd.score,type="l",lty="dotted",lwd=1.25,
	main="D. Dose-response Relationship \n of O3 to Test Scores",
	xlab=bquote("O3" ~ (ppb)),
	ylab="Test Scores (SD)",
	ylim=c(-0.10,0.10),
	cex.axis=0.9,
	cex.lab=0.90,
	cex.main=0.97)
rug(eclsb$o3_2001)
lines(output.o3$xval,output.o3$mt.score,lty="longdash")
labels<-c("Reading Scores","Math Scores")
type<-c("dotted","longdash")
legend("topright",inset=0.03,labels,lty=c("dotted","longdash"),lwd=c(1.25,1),seg.len=2,cex=0.95)

dev.off()

print(startTime)
print(Sys.time())

sink()

