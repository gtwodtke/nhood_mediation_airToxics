################################################
################################################
##                                            ##
## PROGRAM NAME: 28_create_table_S8           ##
## AUTHOR: GW                                 ##
## DESCRIPTION:                               ##
##                                            ##
##  create table of nh effect estimates       ##
##  using alternative tree-based methods      ## 
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
	"xgboost",
	"e1071",
	"SuperLearner",
	"dotwhisker",
	"tidyverse",
	"kableExtra",
	"sys",
	"foreach",
	"doParallel",
	"doRNG",
	"Hmisc")

for(package.i in list.of.packages) {
	suppressPackageStartupMessages(library(package.i, character.only = TRUE))
	}

startTime<-Sys.time()

set.seed(8675309)

#########################
#####               ##### 
#####  LOAD ECLS-B  #####
#####               #####
#########################
eclsb.mi<-read.dta("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\eclsb\\v04_eclsb_mi.dta")
eclsb.mi<-as.data.frame(eclsb.mi[which(eclsb.mi$minum!=0),])
eclsb.mi<-eclsb.mi[order(eclsb.mi$caseid,eclsb.mi$minum),]

nmi<-5
nboot<-200
ntrees<-200
n.cores<-parallel::detectCores()-1
astar<-0.25
a<-0.05

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

##### STANDARDIZE VARS ######
for (v in 1:length(vars.ntxw1)) {
	eclsb.mi[,vars.ntxw1[v]]<-(eclsb.mi[,vars.ntxw1[v]]-weighted.mean(eclsb.mi[,vars.ntxw1[v]],eclsb.mi$sampwt))/sqrt(wtd.var(eclsb.mi[,vars.ntxw1[v]],eclsb.mi$sampwt))
	}

eclsb.mi$readtheta05<-(eclsb.mi$readtheta05-weighted.mean(eclsb.mi$readtheta05,eclsb.mi$sampwt))/sqrt(wtd.var(eclsb.mi$readtheta05,eclsb.mi$sampwt))
eclsb.mi$maththeta05<-(eclsb.mi$maththeta05-weighted.mean(eclsb.mi$maththeta05,eclsb.mi$sampwt))/sqrt(wtd.var(eclsb.mi$maththeta05,eclsb.mi$sampwt))

##### TUNE RF HYPERPARAMETERS #####
registerDoSEQ()

cntrl.tr<-trainControl(method="CV",number=5)

eclsb.tune<-eclsb.mi[which(eclsb.mi$minum==1),]

read.tune<-eclsb.tune[,"readtheta05"]
math.tune<-eclsb.tune[,"maththeta05"]
c.a.tune<-eclsb.tune[,c(vars.base,vars.covw1)]
c.a.m.tune<-eclsb.tune[,c(vars.base,vars.covw1,vars.ntxw1)]

rf2.m1.grid<-expand.grid(
	min.node.size=c(5,10,15,20),
	mtry=floor(length(c(vars.base,vars.covw1))*c(0.3,0.4,0.5,0.6,0.7)),
	sample.frac=c(0.6,0.7,0.8,0.9),                       
	rmse.rd=NA,
	rmse.mt=NA)

rf2.m2.grid<-expand.grid(
	min.node.size=c(5,10,15,20),
	mtry=floor(length(c(vars.base,vars.covw1,vars.ntxw1))*c(0.3,0.4,0.5,0.6,0.7)),
	sample.frac=c(0.6,0.7,0.8,0.9), 
	rmse.rd=NA,
	rmse.mt=NA)

rf2.m3.grid<-expand.grid(
	min.node.size=c(5,10,15,20),
	mtry=floor(length(c(vars.base,vars.covw1))*c(0.3,0.4,0.5,0.6,0.7)),
	sample.frac=c(0.6,0.7,0.8,0.9),                       
	rmse.rd=NA,
	rmse.mt=NA)

for(i in seq_len(nrow(rf2.m1.grid))) {

	rf2.rd.m1<-ranger(
		y=read.tune,
		x=c.a.tune,
		num.trees=ntrees,
		min.node.size=rf2.m1.grid$min.node.size[i],
		mtry=rf2.m1.grid$mtry[i],
		sample.fraction=rf2.m1.grid$sample.frac[i],
		replace=FALSE,
		splitrule="variance",
		case.weights=eclsb.tune$sampwt,
		respect.unordered.factors=TRUE,
		seed=8675309,
		verbose=FALSE)

	rf2.mt.m1<-ranger(
		y=math.tune,
		x=c.a.tune,
		num.trees=ntrees,
		min.node.size=rf2.m1.grid$min.node.size[i],
		mtry=rf2.m1.grid$mtry[i],
		sample.fraction=rf2.m1.grid$sample.frac[i],
		replace=FALSE,
		splitrule="variance",
		case.weights=eclsb.tune$sampwt,
		respect.unordered.factors=TRUE,
		seed=8675309,
		verbose=FALSE)

	rf2.m1.grid$rmse.rd[i]<-sqrt(rf2.rd.m1$prediction.error)
	rf2.m1.grid$rmse.mt[i]<-sqrt(rf2.mt.m1$prediction.error)
	}

for(i in seq_len(nrow(rf2.m2.grid))) {

	rf2.rd.m2<-ranger(
		y=read.tune,
		x=c.a.m.tune,
		num.trees=ntrees,
		min.node.size=rf2.m2.grid$min.node.size[i],
		mtry=rf2.m2.grid$mtry[i],
		sample.fraction=rf2.m2.grid$sample.frac[i],
		replace=FALSE,
		splitrule="variance",
		case.weights=eclsb.tune$sampwt,
		respect.unordered.factors=TRUE,
		seed=8675309,
		verbose=FALSE)

	rf2.mt.m2<-ranger(
		y=math.tune,
		x=c.a.m.tune,
		num.trees=ntrees,
		min.node.size=rf2.m2.grid$min.node.size[i],
		mtry=rf2.m2.grid$mtry[i],
		sample.fraction=rf2.m2.grid$sample.frac[i],
		replace=FALSE,
		splitrule="variance",
		case.weights=eclsb.tune$sampwt,
		respect.unordered.factors=TRUE,
		seed=8675309,
		verbose=FALSE)

	rf2.m2.grid$rmse.rd[i]<-sqrt(rf2.rd.m2$prediction.error)
	rf2.m2.grid$rmse.mt[i]<-sqrt(rf2.mt.m2$prediction.error)
	}

rf2.rd.m1.bestTune<-rf2.m1.grid[which(rf2.m1.grid$rmse.rd==min(rf2.m1.grid$rmse.rd)),c("min.node.size","mtry","sample.frac")]
rf2.mt.m1.bestTune<-rf2.m1.grid[which(rf2.m1.grid$rmse.mt==min(rf2.m1.grid$rmse.mt)),c("min.node.size","mtry","sample.frac")]
rf2.rd.m2.bestTune<-rf2.m2.grid[which(rf2.m2.grid$rmse.rd==min(rf2.m2.grid$rmse.rd)),c("min.node.size","mtry","sample.frac")]
rf2.mt.m2.bestTune<-rf2.m2.grid[which(rf2.m2.grid$rmse.mt==min(rf2.m2.grid$rmse.mt)),c("min.node.size","mtry","sample.frac")]

rf2.rd.m2.tune<-ranger(readtheta05~.,
	data=eclsb.tune[,c("readtheta05",vars.base,vars.covw1,vars.ntxw1)],
	num.trees=ntrees,
	min.node.size=rf2.rd.m2.bestTune[,"min.node.size"],
	mtry=rf2.rd.m2.bestTune[,"mtry"],
	sample.fraction=rf2.rd.m2.bestTune[,"sample.frac"],
	splitrule="variance",
	respect.unordered.factors=TRUE,
	case.weights=eclsb.tune$sampwt,
	seed=8675309)

rf2.mt.m2.tune<-ranger(maththeta05~.,
	data=eclsb.tune[,c("maththeta05",vars.base,vars.covw1,vars.ntxw1)],
	num.trees=ntrees,
	min.node.size=rf2.mt.m2.bestTune[,"min.node.size"],
	mtry=rf2.mt.m2.bestTune[,"mtry"],
	sample.fraction=rf2.mt.m2.bestTune[,"sample.frac"],
	splitrule="variance",
	respect.unordered.factors=TRUE,
	case.weights=eclsb.tune$sampwt,
	seed=8675309)

eclsb.tune.imp<-eclsb.tune[,c(vars.base,vars.covw1,vars.ntxw1)]
eclsb.tune.imp$nhpovrt01<-astar	
read.astar.m.tune<-predict(rf2.rd.m2.tune,eclsb.tune.imp)$pred
math.astar.m.tune<-predict(rf2.mt.m2.tune,eclsb.tune.imp)$pred

for(i in seq_len(nrow(rf2.m3.grid))) {

	rf2.rd.m3<-ranger(
		y=read.astar.m.tune,
		x=c.a.tune,
		num.trees=ntrees,
		min.node.size=rf2.m3.grid$min.node.size[i],
		mtry=rf2.m3.grid$mtry[i],
		sample.fraction=rf2.m3.grid$sample.frac[i],
		replace=FALSE,
		splitrule="variance",
		case.weights=eclsb.tune$sampwt,
		respect.unordered.factors=TRUE,
		seed=8675309,
		verbose=FALSE)

	rf2.mt.m3<-ranger(
		y=math.astar.m.tune,
		x=c.a.tune,
		num.trees=ntrees,
		min.node.size=rf2.m3.grid$min.node.size[i],
		mtry=rf2.m3.grid$mtry[i],
		sample.fraction=rf2.m3.grid$sample.frac[i],
		replace=FALSE,
		splitrule="variance",
		case.weights=eclsb.tune$sampwt,
		respect.unordered.factors=TRUE,
		seed=8675309,
		verbose=FALSE)

	rf2.m3.grid$rmse.rd[i]<-sqrt(rf2.rd.m3$prediction.error)
	rf2.m3.grid$rmse.mt[i]<-sqrt(rf2.mt.m3$prediction.error)
	}

rf2.rd.m3.bestTune<-rf2.m3.grid[which(rf2.m3.grid$rmse.rd==min(rf2.m3.grid$rmse.rd)),c("min.node.size","mtry","sample.frac")]
rf2.mt.m3.bestTune<-rf2.m3.grid[which(rf2.m3.grid$rmse.mt==min(rf2.m3.grid$rmse.mt)),c("min.node.size","mtry","sample.frac")]

rfo.m1.grid<-expand.grid(
	min.node.size=c(5,10,15,20),
	mtry=floor(length(c(vars.base,vars.covw1))*c(0.3,0.4,0.5,0.6,0.7)),
	splitrule="variance")

rfo.m2.grid<-expand.grid(
	min.node.size=c(5,10,15,20),
	mtry=floor(length(c(vars.base,vars.covw1,vars.ntxw1))*c(0.3,0.4,0.5,0.6,0.7)),
	splitrule="variance")

rfo.m3.grid<-expand.grid(
	min.node.size=c(5,10,15,20),
	mtry=floor(length(c(vars.base,vars.covw1))*c(0.3,0.4,0.5,0.6,0.7)),
	splitrule="variance")

rfo.rd.m1<-train(readtheta05~.,
	data=eclsb.tune[,c("readtheta05",vars.base,vars.covw1)],
	method="ranger",
	tuneGrid=rfo.m1.grid,
	trControl=cntrl.tr,
	metric="RMSE",
	respect.unordered.factors=TRUE,
	num.trees=ntrees,
	weights=eclsb.tune$sampwt,
	seed=8675309)

rfo.mt.m1<-train(maththeta05~.,
	data=eclsb.tune[,c("maththeta05",vars.base,vars.covw1)],
	method="ranger",
	tuneGrid=rfo.m1.grid,
	trControl=cntrl.tr,
	metric="RMSE",
	respect.unordered.factors=TRUE,
	num.trees=ntrees,
	weights=eclsb.tune$sampwt,
	seed=8675309)

rfo.rd.m2<-train(readtheta05~.,
	data=eclsb.tune[,c("readtheta05",vars.base,vars.covw1,vars.ntxw1)],
	method="ranger",
	tuneGrid=rfo.m2.grid,
	trControl=cntrl.tr,
	metric="RMSE",
	respect.unordered.factors=TRUE,
	num.trees=ntrees,
	weights=eclsb.tune$sampwt,
	seed=8675309)

rfo.mt.m2<-train(maththeta05~.,
	data=eclsb.tune[,c("maththeta05",vars.base,vars.covw1,vars.ntxw1)],
	method="ranger",
	tuneGrid=rfo.m2.grid,
	trControl=cntrl.tr,
	metric="RMSE",
	respect.unordered.factors=TRUE,
	num.trees=ntrees,
	weights=eclsb.tune$sampwt,
	seed=8675309)

rfo.rd.m2.tune<-ranger(readtheta05~.,
	data=eclsb.tune[,c("readtheta05",vars.base,vars.covw1,vars.ntxw1)],
	num.trees=ntrees,
	min.node.size=rfo.rd.m2$bestTune[,"min.node.size"],
	mtry=rfo.rd.m2$bestTune[,"mtry"],
	splitrule="variance",
	respect.unordered.factors=TRUE,
	case.weights=eclsb.tune$sampwt,
	seed=8675309)

rfo.mt.m2.tune<-ranger(maththeta05~.,
	data=eclsb.tune[,c("maththeta05",vars.base,vars.covw1,vars.ntxw1)],
	num.trees=ntrees,
	min.node.size=rfo.mt.m2$bestTune[,"min.node.size"],
	mtry=rfo.mt.m2$bestTune[,"mtry"],
	splitrule="variance",
	respect.unordered.factors=TRUE,
	case.weights=eclsb.tune$sampwt,
	seed=8675309)

eclsb.tune.imp<-eclsb.tune[,c(vars.base,vars.covw1,vars.ntxw1)]
eclsb.tune.imp$nhpovrt01<-astar	
read.astar.m.tune<-predict(rfo.rd.m2.tune,eclsb.tune.imp)$pred
math.astar.m.tune<-predict(rfo.mt.m2.tune,eclsb.tune.imp)$pred

rfo.rd.m3<-train(
	y=read.astar.m.tune,
	x=eclsb.tune[,c(vars.base,vars.covw1)],
	method="ranger",
	tuneGrid=rfo.m3.grid,
	trControl=cntrl.tr,
	metric="RMSE",
	respect.unordered.factors=TRUE,
	num.trees=ntrees,
	weights=eclsb.tune$sampwt,
	seed=8675309)

rfo.mt.m3<-train(
	y=math.astar.m.tune,
	x=eclsb.tune[,c(vars.base,vars.covw1)],
	method="ranger",
	tuneGrid=rfo.m3.grid,
	trControl=cntrl.tr,
	metric="RMSE",
	respect.unordered.factors=TRUE,
	num.trees=ntrees,
	weights=eclsb.tune$sampwt,
	seed=8675309)

##### CREATE LEARNERS #####
cntrl.sl.cv4<-SuperLearner.CV.control(V=4)
cntrl.sl.cv2<-SuperLearner.CV.control(V=2)

sl.rf1.m1<-create.Learner(
	"SL.ranger",
	name_prefix="sl.rf1.m1",
	params=list(
		num.trees=ntrees,
		min.node.size=5,
		mtry=floor(length(c(vars.base,vars.covw1))/3),
		respect.unordered.factors=TRUE,
		splitrule="variance"))

sl.rf2.rd.m1<-create.Learner(
	"SL.ranger",
	name_prefix="sl.rf2.rd.m1",
	params=list(
		num.trees=ntrees,
		min.node.size=rf2.rd.m1.bestTune[1,"min.node.size"],
		mtry=rf2.rd.m1.bestTune[1,"mtry"],
		sample.fraction=rf2.rd.m1.bestTune[1,"sample.frac"],
		replace=FALSE,
		respect.unordered.factors=TRUE,
		splitrule="variance"))

sl.rf2.mt.m1<-create.Learner(
	"SL.ranger",
	name_prefix="sl.rf2.mt.m1",
	params=list(
		num.trees=ntrees,
		min.node.size=rf2.mt.m1.bestTune[1,"min.node.size"],
		mtry=rf2.mt.m1.bestTune[1,"mtry"],
		sample.fraction=rf2.mt.m1.bestTune[1,"sample.frac"],
		replace=FALSE,
		respect.unordered.factors=TRUE,
		splitrule="variance"))

sl.rfo.rd.m1<-create.Learner(
	"SL.ranger",
	name_prefix="sl.rfo.rd.m1",
	params=list(
		num.trees=ntrees,
		min.node.size=rfo.rd.m1$bestTune[,"min.node.size"],
		mtry=rfo.rd.m1$bestTune[,"mtry"],
		respect.unordered.factors=TRUE,
		splitrule="variance"))

sl.rfo.mt.m1<-create.Learner(
	"SL.ranger",
	name_prefix="sl.rfo.mt.m1",
	params=list(
		num.trees=ntrees,
		min.node.size=rfo.mt.m1$bestTune[,"min.node.size"],
		mtry=rfo.mt.m1$bestTune[,"mtry"],
		respect.unordered.factors=TRUE,
		splitrule="variance"))

sl.rf1.m2<-create.Learner(
	"SL.ranger",
	name_prefix="sl.rf1.m2",
	params=list(
		num.trees=ntrees,
		min.node.size=5,
		mtry=floor(length(c(vars.base,vars.covw1,vars.ntxw1))/3),
		respect.unordered.factors=TRUE,
		splitrule="variance"))	

sl.rf2.rd.m2<-create.Learner(
	"SL.ranger",
	name_prefix="sl.rf2.rd.m2",
	params=list(
		num.trees=ntrees,
		min.node.size=rf2.rd.m2.bestTune[1,"min.node.size"],
		mtry=rf2.rd.m2.bestTune[1,"mtry"],
		sample.fraction=rf2.rd.m2.bestTune[1,"sample.frac"],
		replace=FALSE,
		respect.unordered.factors=TRUE,
		splitrule="variance"))

sl.rf2.mt.m2<-create.Learner(
	"SL.ranger",
	name_prefix="sl.rf2.mt.m2",
	params=list(
		num.trees=ntrees,
		min.node.size=rf2.mt.m2.bestTune[1,"min.node.size"],
		mtry=rf2.mt.m2.bestTune[1,"mtry"],
		sample.fraction=rf2.mt.m2.bestTune[1,"sample.frac"],
		replace=FALSE,
		respect.unordered.factors=TRUE,
		splitrule="variance"))

sl.rfo.rd.m2<-create.Learner(
	"SL.ranger",
	name_prefix="sl.rfo.rd.m2",
	params=list(
		num.trees=ntrees,
		min.node.size=rfo.rd.m2$bestTune[,"min.node.size"],
		mtry=rfo.rd.m2$bestTune[,"mtry"],
		respect.unordered.factors=TRUE,
		splitrule="variance"))

sl.rfo.mt.m2<-create.Learner(
	"SL.ranger",
	name_prefix="sl.rfo.mt.m2",
	params=list(
		num.trees=ntrees,
		min.node.size=rfo.mt.m2$bestTune[,"min.node.size"],
		mtry=rfo.mt.m2$bestTune[,"mtry"],
		respect.unordered.factors=TRUE,
		splitrule="variance"))

sl.rf1.m3<-create.Learner(
	"SL.ranger",
	name_prefix="sl.rf1.m3",
	params=list(
		num.trees=ntrees,
		min.node.size=5,
		mtry=floor(length(c(vars.base,vars.covw1))/3),
		respect.unordered.factors=TRUE,
		splitrule="variance"))

sl.rf2.rd.m3<-create.Learner(
	"SL.ranger",
	name_prefix="sl.rf2.rd.m3",
	params=list(
		num.trees=ntrees,
		min.node.size=rf2.rd.m3.bestTune[1,"min.node.size"],
		mtry=rf2.rd.m3.bestTune[1,"mtry"],
		sample.fraction=rf2.rd.m3.bestTune[1,"sample.frac"],
		replace=FALSE,
		respect.unordered.factors=TRUE,
		splitrule="variance"))

sl.rf2.mt.m3<-create.Learner(
	"SL.ranger",
	name_prefix="sl.rf2.mt.m3",
	params=list(
		num.trees=ntrees,
		min.node.size=rf2.mt.m3.bestTune[1,"min.node.size"],
		mtry=rf2.mt.m3.bestTune[1,"mtry"],
		sample.fraction=rf2.mt.m3.bestTune[1,"sample.frac"],
		replace=FALSE,
		respect.unordered.factors=TRUE,
		splitrule="variance"))

sl.rfo.rd.m3<-create.Learner(
	"SL.ranger",
	name_prefix="sl.rfo.rd.m3",
	params=list(
		num.trees=ntrees,
		min.node.size=rfo.rd.m3$bestTune[,"min.node.size"],
		mtry=rfo.rd.m3$bestTune[,"mtry"],
		respect.unordered.factors=TRUE,
		splitrule="variance"))

sl.rfo.mt.m3<-create.Learner(
	"SL.ranger",
	name_prefix="sl.rfo.mt.m3",
	params=list(
		num.trees=ntrees,
		min.node.size=rfo.mt.m3$bestTune[,"min.node.size"],
		mtry=rfo.mt.m3$bestTune[,"mtry"],
		respect.unordered.factors=TRUE,
		splitrule="variance"))

##### PARALLELIZATION #####
my.cluster<-parallel::makeCluster(n.cores,type="PSOCK")
doParallel::registerDoParallel(cl=my.cluster)
foreach::getDoParRegistered()
clusterEvalQ(cl=my.cluster, library(SuperLearner))
registerDoRNG(8675309)

##### ESTIMATE NHOOD EFFECTS #####
miest.rd.rf1<-miest.mt.rf1<-matrix(data=NA,nrow=nmi,ncol=3)
miest.rd.rf2<-miest.mt.rf2<-matrix(data=NA,nrow=nmi,ncol=3)
miest.rd.sl1<-miest.mt.sl1<-matrix(data=NA,nrow=nmi,ncol=3)

boot.ci.rd.rf1<-boot.ci.mt.rf1<-NULL
boot.ci.rd.rf2<-boot.ci.mt.rf2<-NULL
boot.ci.rd.sl1<-boot.ci.mt.sl1<-NULL

for (i in 1:nmi) {

	### LOAD MI DATA ###
	print(c("nmi=",i))
	eclsb<-eclsb.mi[which(eclsb.mi$minum==i),]

	### COMPUTE ESTIMATES ###
	read.train<-eclsb[,"readtheta05"]
	math.train<-eclsb[,"maththeta05"]
	c.a.train<-eclsb[,c(vars.base,vars.covw1)]
	c.a.m.train<-eclsb[,c(vars.base,vars.covw1,vars.ntxw1)]
	wts.train<-eclsb[,"sampwt"]

	# STEP 1: ESTIMATE E(Y(astar)) and E(Y(a)) #
	read.giv.c.a<-SuperLearner(
		Y=read.train,
		X=c.a.train,
		SL.library=c("sl.rf1.m1_1","sl.rf2.rd.m1_1","sl.rfo.rd.m1_1"),
		cvControl=cntrl.sl.cv4,
		obsWeights=wts.train)

	math.giv.c.a<-SuperLearner(
		Y=math.train,
		X=c.a.train,
		SL.library=c("sl.rf1.m1_1","sl.rf2.mt.m1_1","sl.rfo.mt.m1_1"),
		cvControl=cntrl.sl.cv4,
		obsWeights=wts.train)

	c.a.imp<-eclsb[,c(vars.base,vars.covw1)]

	c.a.imp$nhpovrt01<-astar	
	read.astar.pred<-predict(read.giv.c.a,c.a.imp)
	math.astar.pred<-predict(math.giv.c.a,c.a.imp)
	read.astar.rf1<-read.astar.pred$library.predict[,1]
	read.astar.rf2<-read.astar.pred$library.predict[,2]
	read.astar.sl1<-read.astar.pred$pred
	math.astar.rf1<-math.astar.pred$library.predict[,1]
	math.astar.rf2<-math.astar.pred$library.predict[,2]
	math.astar.sl1<-math.astar.pred$pred

	c.a.imp$nhpovrt01<-a
	read.a.pred<-predict(read.giv.c.a,c.a.imp)
	math.a.pred<-predict(math.giv.c.a,c.a.imp)
	read.a.rf1<-read.a.pred$library.predict[,1]
	read.a.rf2<-read.a.pred$library.predict[,2]
	read.a.sl1<-read.a.pred$pred
	math.a.rf1<-math.a.pred$library.predict[,1]
	math.a.rf2<-math.a.pred$library.predict[,2]
	math.a.sl1<-math.a.pred$pred

	# STEP 2: ESTIMATE E(Y(astar,m(a)) #
	read.giv.c.a.m<-SuperLearner(
		Y=read.train,
		X=c.a.m.train,
		SL.library=c("sl.rf1.m2_1","sl.rf2.rd.m2_1","sl.rfo.rd.m2_1"),
		cvControl=cntrl.sl.cv4,
		obsWeights=wts.train)

	math.giv.c.a.m<-SuperLearner(
		Y=math.train,
		X=c.a.m.train,
		SL.library=c("sl.rf1.m2_1","sl.rf2.mt.m2_1","sl.rfo.mt.m2_1"),
		cvControl=cntrl.sl.cv4,
		obsWeights=wts.train)

	c.a.m.train.imp<-eclsb[,c(vars.base,vars.covw1,vars.ntxw1)]
	c.a.m.train.imp$nhpovrt01<-astar	
	read.astar.m.pred<-predict(read.giv.c.a.m,c.a.m.train.imp)
	math.astar.m.pred<-predict(math.giv.c.a.m,c.a.m.train.imp)
	read.astar.m.rf1<-read.astar.m.pred$library.predict[,1]
	read.astar.m.rf2<-read.astar.m.pred$library.predict[,2]
	read.astar.m.sl1<-read.astar.m.pred$pred
	math.astar.m.rf1<-math.astar.m.pred$library.predict[,1]
	math.astar.m.rf2<-math.astar.m.pred$library.predict[,2]
	math.astar.m.sl1<-math.astar.m.pred$pred
		
	read.uhat.astar.m.rf1<-SuperLearner(
		Y=read.astar.m.rf1,
		X=c.a.train,
		SL.library="sl.rf1.m3_1",
		cvControl=cntrl.sl.cv2,
		obsWeights=wts.train)

	math.uhat.astar.m.rf1<-SuperLearner(
		Y=math.astar.m.rf1,
		X=c.a.train,
		SL.library="sl.rf1.m3_1",
		cvControl=cntrl.sl.cv2,
		obsWeights=wts.train)

	read.uhat.astar.m.rf2<-SuperLearner(
		Y=read.astar.m.rf2,
		X=c.a.train,
		SL.library="sl.rf2.rd.m3_1",
		cvControl=cntrl.sl.cv2,
		obsWeights=wts.train)

	math.uhat.astar.m.rf2<-SuperLearner(
		Y=math.astar.m.rf2,
		X=c.a.train,
		SL.library="sl.rf2.mt.m3_1",
		cvControl=cntrl.sl.cv2,
		obsWeights=wts.train)

	read.uhat.astar.m.sl1<-SuperLearner(
		Y=read.astar.m.sl1,
		X=c.a.train,
		SL.library=c("sl.rf1.m3_1","sl.rf2.rd.m3_1","sl.rfo.rd.m3_1"),
		cvControl=cntrl.sl.cv4,
		obsWeights=wts.train)

	math.uhat.astar.m.sl1<-SuperLearner(
		Y=math.astar.m.sl1,
		X=c.a.train,
		SL.library=c("sl.rf1.m3_1","sl.rf2.mt.m3_1","sl.rfo.mt.m3_1"),
		cvControl=cntrl.sl.cv4,
		obsWeights=wts.train)

	c.a.imp<-eclsb[,c(vars.base,vars.covw1)]
	c.a.imp$nhpovrt01<-a	
	read.astar.mofa.rf1<-predict(read.uhat.astar.m.rf1,c.a.imp)$pred
	read.astar.mofa.rf2<-predict(read.uhat.astar.m.rf2,c.a.imp)$pred
	read.astar.mofa.sl1<-predict(read.uhat.astar.m.sl1,c.a.imp)$pred
	math.astar.mofa.rf1<-predict(math.uhat.astar.m.rf1,c.a.imp)$pred
	math.astar.mofa.rf2<-predict(math.uhat.astar.m.rf2,c.a.imp)$pred
	math.astar.mofa.sl1<-predict(math.uhat.astar.m.sl1,c.a.imp)$pred

	# STEP 3: COMPUTE ATE, NDE, NIE #
	miest.rd.rf1[i,1]<-weighted.mean(read.astar.rf1,wts.train)-weighted.mean(read.a.rf1,wts.train)
	miest.rd.rf1[i,2]<-weighted.mean(read.astar.mofa.rf1,wts.train)-weighted.mean(read.a.rf1,wts.train)
	miest.rd.rf1[i,3]<-weighted.mean(read.astar.rf1,wts.train)-weighted.mean(read.astar.mofa.rf1,wts.train)

	miest.mt.rf1[i,1]<-weighted.mean(math.astar.rf1,wts.train)-weighted.mean(math.a.rf1,wts.train)
	miest.mt.rf1[i,2]<-weighted.mean(math.astar.mofa.rf1,wts.train)-weighted.mean(math.a.rf1,wts.train)
	miest.mt.rf1[i,3]<-weighted.mean(math.astar.rf1,wts.train)-weighted.mean(math.astar.mofa.rf1,wts.train)

	miest.rd.rf2[i,1]<-weighted.mean(read.astar.rf2,wts.train)-weighted.mean(read.a.rf2,wts.train)
	miest.rd.rf2[i,2]<-weighted.mean(read.astar.mofa.rf2,wts.train)-weighted.mean(read.a.rf2,wts.train)
	miest.rd.rf2[i,3]<-weighted.mean(read.astar.rf2,wts.train)-weighted.mean(read.astar.mofa.rf2,wts.train)

	miest.mt.rf2[i,1]<-weighted.mean(math.astar.rf2,wts.train)-weighted.mean(math.a.rf2,wts.train)
	miest.mt.rf2[i,2]<-weighted.mean(math.astar.mofa.rf2,wts.train)-weighted.mean(math.a.rf2,wts.train)
	miest.mt.rf2[i,3]<-weighted.mean(math.astar.rf2,wts.train)-weighted.mean(math.astar.mofa.rf2,wts.train)

	miest.rd.sl1[i,1]<-weighted.mean(read.astar.sl1,wts.train)-weighted.mean(read.a.sl1,wts.train)
	miest.rd.sl1[i,2]<-weighted.mean(read.astar.mofa.sl1,wts.train)-weighted.mean(read.a.sl1,wts.train)
	miest.rd.sl1[i,3]<-weighted.mean(read.astar.sl1,wts.train)-weighted.mean(read.astar.mofa.sl1,wts.train)

	miest.mt.sl1[i,1]<-weighted.mean(math.astar.sl1,wts.train)-weighted.mean(math.a.sl1,wts.train)
	miest.mt.sl1[i,2]<-weighted.mean(math.astar.mofa.sl1,wts.train)-weighted.mean(math.a.sl1,wts.train)
	miest.mt.sl1[i,3]<-weighted.mean(math.astar.sl1,wts.train)-weighted.mean(math.astar.mofa.sl1,wts.train)

	### COMPUTE BOOTSTRAP CIs ###
	clusterExport(cl=my.cluster,
		list(
			"vars.ntxw1",
			"vars.covw1",
			"vars.base",
			"nmi",
			"nboot",
			"astar",
			"a",
			"eclsb",
			"sl.rf1.m1",
			"sl.rf2.rd.m1",
			"sl.rf2.mt.m1",
			"sl.rfo.rd.m1",
			"sl.rfo.mt.m1",
			"sl.rf1.m2",
			"sl.rf2.rd.m2",
			"sl.rf2.mt.m2",
			"sl.rfo.rd.m2",
			"sl.rfo.mt.m2",
			"sl.rf1.m3",
			"sl.rf2.rd.m3",
			"sl.rf2.mt.m3",
			"sl.rfo.rd.m3",
			"sl.rfo.mt.m3",
			"sl.rf1.m1_1",
			"sl.rf2.rd.m1_1",
			"sl.rf2.mt.m1_1",
			"sl.rfo.rd.m1_1",
			"sl.rfo.mt.m1_1",
			"sl.rf1.m2_1",
			"sl.rf2.rd.m2_1",
			"sl.rf2.mt.m2_1",
			"sl.rfo.rd.m2_1",
			"sl.rfo.mt.m2_1",
			"sl.rf1.m3_1",
			"sl.rf2.rd.m3_1",
			"sl.rf2.mt.m3_1",
			"sl.rfo.rd.m3_1",
			"sl.rfo.mt.m3_1",
			"cntrl.sl.cv4",
			"cntrl.sl.cv2"),
		envir=environment())

	boot.est.w1<-foreach(h=1:nboot, .combine=cbind) %dopar% {

		boot.eclsb<-NULL

		for (s in 1:length(unique(eclsb$strat))) {
			eclsb.strata<-eclsb[which(eclsb$strat==s),]
			idboot.1<-sample(unique(eclsb.strata$psu),length(unique(eclsb.strata$psu))-1,replace=T)
			idboot.2<-table(idboot.1)
			boot.eclsb.strata<-NULL
			for (g in 1:max(idboot.2)) {
				boot.data<-eclsb.strata[eclsb.strata$psu %in% names(idboot.2[idboot.2 %in% g]),]
				for (l in 1:g) { 
					boot.eclsb.strata<-rbind(boot.eclsb.strata,boot.data) 
					}
				}
		boot.eclsb<-rbind(boot.eclsb,boot.eclsb.strata,boot.eclsb.strata)
		}
			
		boot.read.train<-boot.eclsb[,"readtheta05"]
		boot.math.train<-boot.eclsb[,"maththeta05"]
		boot.c.a.train<-boot.eclsb[,c(vars.base,vars.covw1)]
		boot.c.a.m.train<-boot.eclsb[,c(vars.base,vars.covw1,vars.ntxw1)]
		boot.wts.train<-boot.eclsb[,"sampwt"]

		boot.read.giv.c.a<-SuperLearner(
			Y=boot.read.train,
			X=boot.c.a.train,
			SL.library=c("sl.rf1.m1_1","sl.rf2.rd.m1_1","sl.rfo.rd.m1_1"),
			cvControl=cntrl.sl.cv4,
			obsWeights=boot.wts.train)

		boot.math.giv.c.a<-SuperLearner(
			Y=boot.math.train,
			X=boot.c.a.train,
			SL.library=c("sl.rf1.m1_1","sl.rf2.mt.m1_1","sl.rfo.mt.m1_1"),
			cvControl=cntrl.sl.cv4,
			obsWeights=boot.wts.train)

		boot.c.a.imp<-boot.eclsb[,c(vars.base,vars.covw1)]

		boot.c.a.imp$nhpovrt01<-astar	
		boot.read.astar.pred<-predict(boot.read.giv.c.a,boot.c.a.imp)
		boot.math.astar.pred<-predict(boot.math.giv.c.a,boot.c.a.imp)
		boot.read.astar.rf1<-boot.read.astar.pred$library.predict[,1]
		boot.read.astar.rf2<-boot.read.astar.pred$library.predict[,2]
		boot.read.astar.sl1<-boot.read.astar.pred$pred
		boot.math.astar.rf1<-boot.math.astar.pred$library.predict[,1]
		boot.math.astar.rf2<-boot.math.astar.pred$library.predict[,2]
		boot.math.astar.sl1<-boot.math.astar.pred$pred

		boot.c.a.imp$nhpovrt01<-a
		boot.read.a.pred<-predict(boot.read.giv.c.a,boot.c.a.imp)
		boot.math.a.pred<-predict(boot.math.giv.c.a,boot.c.a.imp)
		boot.read.a.rf1<-boot.read.a.pred$library.predict[,1]
		boot.read.a.rf2<-boot.read.a.pred$library.predict[,2]
		boot.read.a.sl1<-boot.read.a.pred$pred
		boot.math.a.rf1<-boot.math.a.pred$library.predict[,1]
		boot.math.a.rf2<-boot.math.a.pred$library.predict[,2]
		boot.math.a.sl1<-boot.math.a.pred$pred
	
		boot.read.giv.c.a.m<-SuperLearner(
			Y=boot.read.train,
			X=boot.c.a.m.train,
			SL.library=c("sl.rf1.m2_1","sl.rf2.rd.m2_1","sl.rfo.rd.m2_1"),
			cvControl=cntrl.sl.cv4,
			obsWeights=boot.wts.train)

		boot.math.giv.c.a.m<-SuperLearner(
			Y=boot.math.train,
			X=boot.c.a.m.train,
			SL.library=c("sl.rf1.m2_1","sl.rf2.mt.m2_1","sl.rfo.mt.m2_1"),
			cvControl=cntrl.sl.cv4,
			obsWeights=boot.wts.train)

		boot.c.a.m.train.imp<-boot.eclsb[,c(vars.base,vars.covw1,vars.ntxw1)]
		boot.c.a.m.train.imp$nhpovrt01<-astar	
		boot.read.astar.m.pred<-predict(boot.read.giv.c.a.m,boot.c.a.m.train.imp)
		boot.math.astar.m.pred<-predict(boot.math.giv.c.a.m,boot.c.a.m.train.imp)
		boot.read.astar.m.rf1<-boot.read.astar.m.pred$library.predict[,1]
		boot.read.astar.m.rf2<-boot.read.astar.m.pred$library.predict[,2]
		boot.read.astar.m.sl1<-boot.read.astar.m.pred$pred
		boot.math.astar.m.rf1<-boot.math.astar.m.pred$library.predict[,1]
		boot.math.astar.m.rf2<-boot.math.astar.m.pred$library.predict[,2]
		boot.math.astar.m.sl1<-boot.math.astar.m.pred$pred
		
		boot.read.uhat.astar.m.rf1<-SuperLearner(
			Y=boot.read.astar.m.rf1,
			X=boot.c.a.train,
			SL.library="sl.rf1.m3_1",
			cvControl=cntrl.sl.cv2,
			obsWeights=boot.wts.train)

		boot.math.uhat.astar.m.rf1<-SuperLearner(
			Y=boot.math.astar.m.rf1,
			X=boot.c.a.train,
			SL.library="sl.rf1.m3_1",
			cvControl=cntrl.sl.cv2,
			obsWeights=boot.wts.train)

		boot.read.uhat.astar.m.rf2<-SuperLearner(
			Y=boot.read.astar.m.rf2,
			X=boot.c.a.train,
			SL.library="sl.rf2.rd.m3_1",
			cvControl=cntrl.sl.cv2,
			obsWeights=boot.wts.train)

		boot.math.uhat.astar.m.rf2<-SuperLearner(
			Y=boot.math.astar.m.rf2,
			X=boot.c.a.train,
			SL.library="sl.rf2.mt.m3_1",
			cvControl=cntrl.sl.cv2,
			obsWeights=boot.wts.train)

		boot.read.uhat.astar.m.sl1<-SuperLearner(
			Y=boot.read.astar.m.sl1,
			X=boot.c.a.train,
			SL.library=c("sl.rf1.m3_1","sl.rf2.rd.m3_1","sl.rfo.rd.m3_1"),
			cvControl=cntrl.sl.cv4,
			obsWeights=boot.wts.train)

		boot.math.uhat.astar.m.sl1<-SuperLearner(
			Y=boot.math.astar.m.sl1,
			X=boot.c.a.train,
			SL.library=c("sl.rf1.m3_1","sl.rf2.mt.m3_1","sl.rfo.mt.m3_1"),
			cvControl=cntrl.sl.cv4,
			obsWeights=boot.wts.train)

		boot.c.a.imp<-boot.eclsb[,c(vars.base,vars.covw1)]
		boot.c.a.imp$nhpovrt01<-a	
		boot.read.astar.mofa.rf1<-predict(boot.read.uhat.astar.m.rf1,boot.c.a.imp)$pred
		boot.read.astar.mofa.rf2<-predict(boot.read.uhat.astar.m.rf2,boot.c.a.imp)$pred
		boot.read.astar.mofa.sl1<-predict(boot.read.uhat.astar.m.sl1,boot.c.a.imp)$pred
		boot.math.astar.mofa.rf1<-predict(boot.math.uhat.astar.m.rf1,boot.c.a.imp)$pred
		boot.math.astar.mofa.rf2<-predict(boot.math.uhat.astar.m.rf2,boot.c.a.imp)$pred
		boot.math.astar.mofa.sl1<-predict(boot.math.uhat.astar.m.sl1,boot.c.a.imp)$pred

		boot.ate.rd.rf1<-weighted.mean(boot.read.astar.rf1,boot.wts.train)-weighted.mean(boot.read.a.rf1,boot.wts.train)
		boot.nde.rd.rf1<-weighted.mean(boot.read.astar.mofa.rf1,boot.wts.train)-weighted.mean(boot.read.a.rf1,boot.wts.train)
		boot.nie.rd.rf1<-weighted.mean(boot.read.astar.rf1,boot.wts.train)-weighted.mean(boot.read.astar.mofa.rf1,boot.wts.train)
		boot.ate.mt.rf1<-weighted.mean(boot.math.astar.rf1,boot.wts.train)-weighted.mean(boot.math.a.rf1,boot.wts.train)
		boot.nde.mt.rf1<-weighted.mean(boot.math.astar.mofa.rf1,boot.wts.train)-weighted.mean(boot.math.a.rf1,boot.wts.train)
		boot.nie.mt.rf1<-weighted.mean(boot.math.astar.rf1,boot.wts.train)-weighted.mean(boot.math.astar.mofa.rf1,boot.wts.train)

		boot.ate.rd.rf2<-weighted.mean(boot.read.astar.rf2,boot.wts.train)-weighted.mean(boot.read.a.rf2,boot.wts.train)
		boot.nde.rd.rf2<-weighted.mean(boot.read.astar.mofa.rf2,boot.wts.train)-weighted.mean(boot.read.a.rf2,boot.wts.train)
		boot.nie.rd.rf2<-weighted.mean(boot.read.astar.rf2,boot.wts.train)-weighted.mean(boot.read.astar.mofa.rf2,boot.wts.train)
		boot.ate.mt.rf2<-weighted.mean(boot.math.astar.rf2,boot.wts.train)-weighted.mean(boot.math.a.rf2,boot.wts.train)
		boot.nde.mt.rf2<-weighted.mean(boot.math.astar.mofa.rf2,boot.wts.train)-weighted.mean(boot.math.a.rf2,boot.wts.train)
		boot.nie.mt.rf2<-weighted.mean(boot.math.astar.rf2,boot.wts.train)-weighted.mean(boot.math.astar.mofa.rf2,boot.wts.train)

		boot.ate.rd.sl1<-weighted.mean(boot.read.astar.sl1,boot.wts.train)-weighted.mean(boot.read.a.sl1,boot.wts.train)
		boot.nde.rd.sl1<-weighted.mean(boot.read.astar.mofa.sl1,boot.wts.train)-weighted.mean(boot.read.a.sl1,boot.wts.train)
		boot.nie.rd.sl1<-weighted.mean(boot.read.astar.sl1,boot.wts.train)-weighted.mean(boot.read.astar.mofa.sl1,boot.wts.train)
		boot.ate.mt.sl1<-weighted.mean(boot.math.astar.sl1,boot.wts.train)-weighted.mean(boot.math.a.sl1,boot.wts.train)
		boot.nde.mt.sl1<-weighted.mean(boot.math.astar.mofa.sl1,boot.wts.train)-weighted.mean(boot.math.a.sl1,boot.wts.train)
		boot.nie.mt.sl1<-weighted.mean(boot.math.astar.sl1,boot.wts.train)-weighted.mean(boot.math.astar.mofa.sl1,boot.wts.train)

		return(
			list(
				boot.ate.rd.rf1,boot.nde.rd.rf1,boot.nie.rd.rf1,
				boot.ate.mt.rf1,boot.nde.mt.rf1,boot.nie.mt.rf1,
				boot.ate.rd.rf2,boot.nde.rd.rf2,boot.nie.rd.rf2,
				boot.ate.mt.rf2,boot.nde.mt.rf2,boot.nie.mt.rf2,
				boot.ate.rd.sl1,boot.nde.rd.sl1,boot.nie.rd.sl1,
				boot.ate.mt.sl1,boot.nde.mt.sl1,boot.nie.mt.sl1))
		}

	boot.est.w1<-matrix(unlist(boot.est.w1),ncol=18,byrow=TRUE)

	boot.ci.rd.rf1<-rbind(boot.ci.rd.rf1,boot.est.w1[,1:3])
	boot.ci.mt.rf1<-rbind(boot.ci.mt.rf1,boot.est.w1[,4:6])

	boot.ci.rd.rf2<-rbind(boot.ci.rd.rf2,boot.est.w1[,7:9])
	boot.ci.mt.rf2<-rbind(boot.ci.mt.rf2,boot.est.w1[,10:12])

	boot.ci.rd.sl1<-rbind(boot.ci.rd.sl1,boot.est.w1[,13:15])
	boot.ci.mt.sl1<-rbind(boot.ci.mt.sl1,boot.est.w1[,16:18])
	}

stopCluster(my.cluster)
rm(my.cluster)

### COMBINE MI ESTIMATES ###
est.rd.rf1<-est.mt.rf1<-matrix(data=NA,nrow=3,ncol=3)
est.rd.rf2<-est.mt.rf2<-matrix(data=NA,nrow=3,ncol=3)
est.rd.sl1<-est.mt.sl1<-matrix(data=NA,nrow=3,ncol=3)

for (i in 1:3) { 

	est.rd.rf1[i,1]<-round(mean(miest.rd.rf1[,i]),digits=3)
	est.rd.rf1[i,2]<-round(quantile(boot.ci.rd.rf1[,i],prob=0.025),digits=3)
	est.rd.rf1[i,3]<-round(quantile(boot.ci.rd.rf1[,i],prob=0.975),digits=3)

	est.mt.rf1[i,1]<-round(mean(miest.mt.rf1[,i]),digits=3)
	est.mt.rf1[i,2]<-round(quantile(boot.ci.mt.rf1[,i],prob=0.025),digits=3)
	est.mt.rf1[i,3]<-round(quantile(boot.ci.mt.rf1[,i],prob=0.975),digits=3)

	est.rd.rf2[i,1]<-round(mean(miest.rd.rf2[,i]),digits=3)
	est.rd.rf2[i,2]<-round(quantile(boot.ci.rd.rf2[,i],prob=0.025),digits=3)
	est.rd.rf2[i,3]<-round(quantile(boot.ci.rd.rf2[,i],prob=0.975),digits=3)

	est.mt.rf2[i,1]<-round(mean(miest.mt.rf2[,i]),digits=3)
	est.mt.rf2[i,2]<-round(quantile(boot.ci.mt.rf2[,i],prob=0.025),digits=3)
	est.mt.rf2[i,3]<-round(quantile(boot.ci.mt.rf2[,i],prob=0.975),digits=3)
	
	est.rd.sl1[i,1]<-round(mean(miest.rd.sl1[,i]),digits=3)
	est.rd.sl1[i,2]<-round(quantile(boot.ci.rd.sl1[,i],prob=0.025),digits=3)
	est.rd.sl1[i,3]<-round(quantile(boot.ci.rd.sl1[,i],prob=0.975),digits=3)

	est.mt.sl1[i,1]<-round(mean(miest.mt.sl1[,i]),digits=3)
	est.mt.sl1[i,2]<-round(quantile(boot.ci.mt.sl1[,i],prob=0.025),digits=3)
	est.mt.sl1[i,3]<-round(quantile(boot.ci.mt.sl1[,i],prob=0.975),digits=3)
	}

rlabel<-c('ATE','NDE','NIE')

output.rd.rf1<-data.frame(est.rd.rf1,row.names=rlabel)
output.mt.rf1<-data.frame(est.mt.rf1,row.names=rlabel)
output.rd.rf2<-data.frame(est.rd.rf2,row.names=rlabel)
output.mt.rf2<-data.frame(est.mt.rf2,row.names=rlabel)
output.rd.sl1<-data.frame(est.rd.sl1,row.names=rlabel)
output.mt.sl1<-data.frame(est.mt.sl1,row.names=rlabel)

colnames(output.rd.rf1)<-colnames(output.mt.rf1)<-c('estimate','ll.pct.95ci','ul.pct.95ci')
colnames(output.rd.rf2)<-colnames(output.mt.rf2)<-c('estimate','ll.pct.95ci','ul.pct.95ci')
colnames(output.rd.sl1)<-colnames(output.mt.sl1)<-c('estimate','ll.pct.95ci','ul.pct.95ci')

###########################
#####                 #####
#####  PRINT RESULTS  #####
#####                 #####
###########################
sink("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\programs\\_LOGS\\28_create_table_S8_log.txt")

cat("===========================================\n")
cat("RF W/ DEFAULT HYPERPARAMETERS\n")
cat("===========================================\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Reading Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.rd.rf1)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Math Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.mt.rf1)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("===========================================\n")
cat(" \n")
cat(" \n")
cat(" \n")
cat("===========================================\n")
cat("RF W/ REPLACEMENT SAMPLING FRACTION S\n")
cat("===========================================\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Reading Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.rd.rf2)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Math Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.mt.rf2)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("===========================================\n")
cat(" \n")
cat(" \n")
cat(" \n")
cat("===========================================\n")
cat("SUPER LEARNER\n")
cat("===========================================\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Reading Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.rd.sl1)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Math Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.mt.sl1)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("===========================================\n")

print(startTime)
print(Sys.time())

sink()
