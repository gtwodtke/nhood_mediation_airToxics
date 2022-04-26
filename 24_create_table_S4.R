################################################
################################################
##                                            ##
## PROGRAM NAME: 24_create_table_S4           ##
## AUTHOR: GW                                 ##
## DESCRIPTION:                               ##
##                                            ##
##  create table of nh effect estimates       ##
##  based on alternative Tx contrasts         ## 
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

##### LOAD ECLS-B #####
eclsb.mi<-read.dta("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\eclsb\\v04_eclsb_mi.dta")
eclsb.mi<-as.data.frame(eclsb.mi[which(eclsb.mi$minum!=0),])
eclsb.mi<-eclsb.mi[order(eclsb.mi$caseid,eclsb.mi$minum),]

nmi<-5
nboot<-1000
ntrees<-200
n.cores<-parallel::detectCores()-1

a5<-0.05
a10<-0.10
a15<-0.15
a20<-0.20
a30<-0.30

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

##### ESTIMATE MARGINAL EFFECTS #####
miest.rd.10v5<-miest.mt.10v5<-matrix(data=NA,nrow=nmi,ncol=3)
miest.rd.15v5<-miest.mt.15v5<-matrix(data=NA,nrow=nmi,ncol=3)
miest.rd.20v5<-miest.mt.20v5<-matrix(data=NA,nrow=nmi,ncol=3)
miest.rd.30v5<-miest.mt.30v5<-matrix(data=NA,nrow=nmi,ncol=3)

boot.ci.rd.10v5<-boot.ci.mt.10v5<-NULL
boot.ci.rd.15v5<-boot.ci.mt.15v5<-NULL
boot.ci.rd.20v5<-boot.ci.mt.20v5<-NULL
boot.ci.rd.30v5<-boot.ci.mt.30v5<-NULL

for (i in 1:nmi) {

	### LOAD MI DATA ###
	print(i)
	eclsb<-eclsb.mi[which(eclsb.mi$minum==i),]

	### TUNE HYPERPARAMETERS ###
	registerDoSEQ()

	cntrl<-trainControl(method="CV",number=5)

	rf.m1.grid<-expand.grid(
		min.node.size=c(5,10,15,20),
		mtry=floor(length(c(vars.base,vars.covw1))*c(0.3,0.4,0.5,0.6,0.7)),
		splitrule="variance")

	rf.m2.grid<-expand.grid(
		min.node.size=c(5,10,15,20),
		mtry=floor(length(c(vars.base,vars.covw1,vars.ntxw1))*c(0.3,0.4,0.5,0.6,0.7)),
		splitrule="variance")

	rf.m3.grid<-expand.grid(
		min.node.size=c(5,10,15,20),
		mtry=floor(length(c(vars.base,vars.covw1))*c(0.3,0.4,0.5,0.6,0.7)),
		splitrule="variance")

	rf.rd.m1<-train(readtheta05~.,
		data=eclsb[,c("readtheta05",vars.base,vars.covw1)],
		method="ranger",
		tuneGrid=rf.m1.grid,
		trControl=cntrl,
		metric="RMSE",
		respect.unordered.factors=TRUE,
		num.trees=ntrees,
		weights=eclsb$sampwt,
		seed=8675309)

	rf.mt.m1<-train(maththeta05~.,
		data=eclsb[,c("maththeta05",vars.base,vars.covw1)],
		method="ranger",
		tuneGrid=rf.m1.grid,
		trControl=cntrl,
		metric="RMSE",
		respect.unordered.factors=TRUE,
		num.trees=ntrees,
		weights=eclsb$sampwt,
		seed=8675309)

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

	rd.m2.tune<-ranger(readtheta05~.,
		data=eclsb[,c("readtheta05",vars.base,vars.covw1,vars.ntxw1)],
		num.trees=ntrees,
		min.node.size=rf.rd.m2$bestTune[,"min.node.size"],
		mtry=rf.rd.m2$bestTune[,"mtry"],
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=eclsb$sampwt,
		seed=8675309)

	mt.m2.tune<-ranger(maththeta05~.,
		data=eclsb[,c("maththeta05",vars.base,vars.covw1,vars.ntxw1)],
		num.trees=ntrees,
		min.node.size=rf.mt.m2$bestTune[,"min.node.size"],
		mtry=rf.mt.m2$bestTune[,"mtry"],
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=eclsb$sampwt,
		seed=8675309)

	eclsb.tune.imp<-eclsb[,c(vars.base,vars.covw1,vars.ntxw1)]
	eclsb.tune.imp$nhpovrt01<-a20	
	read.astar.m.tune<-predict(rd.m2.tune,eclsb.tune.imp)$pred
	math.astar.m.tune<-predict(mt.m2.tune,eclsb.tune.imp)$pred

	rf.rd.m3<-train(
		y=read.astar.m.tune,
		x=eclsb[,c(vars.base,vars.covw1)],
		method="ranger",
		tuneGrid=rf.m3.grid,
		trControl=cntrl,
		metric="RMSE",
		respect.unordered.factors=TRUE,
		num.trees=ntrees,
		weights=eclsb$sampwt,
		seed=8675309)

	rf.mt.m3<-train(
		y=math.astar.m.tune,
		x=eclsb[,c(vars.base,vars.covw1)],
		method="ranger",
		tuneGrid=rf.m3.grid,
		trControl=cntrl,
		metric="RMSE",
		respect.unordered.factors=TRUE,
		num.trees=ntrees,
		weights=eclsb$sampwt,
		seed=8675309)

	### COMPUTE ESTIMATES ###
	read.train<-eclsb[,"readtheta05"]
	math.train<-eclsb[,"maththeta05"]
	c.a.train<-eclsb[,c(vars.base,vars.covw1)]
	c.a.m.train<-eclsb[,c(vars.base,vars.covw1,vars.ntxw1)]
	wts.train<-eclsb[,"sampwt"]

	# STEP 1: ESTIMATE E(Y(astar)) and E(Y(a)) #
	read.giv.c.a<-ranger(
		y=read.train,
		x=c.a.train,
		num.trees=ntrees,
		min.node.size=rf.rd.m1$bestTune[,"min.node.size"],
		mtry=rf.rd.m1$bestTune[,"mtry"],
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=wts.train,
		seed=8675309)

	math.giv.c.a<-ranger(
		y=math.train,
		x=c.a.train,
		num.trees=ntrees,
		min.node.size=rf.mt.m1$bestTune[,"min.node.size"],
		mtry=rf.mt.m1$bestTune[,"mtry"],
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=wts.train,
		seed=8675309)

	c.a.imp<-eclsb[,c(vars.base,vars.covw1)]
	c.a.imp$nhpovrt01<-a5	
	read.a5<-predict(read.giv.c.a,c.a.imp)$pred
	math.a5<-predict(math.giv.c.a,c.a.imp)$pred
	c.a.imp$nhpovrt01<-a10	
	read.a10<-predict(read.giv.c.a,c.a.imp)$pred
	math.a10<-predict(math.giv.c.a,c.a.imp)$pred
	c.a.imp$nhpovrt01<-a15	
	read.a15<-predict(read.giv.c.a,c.a.imp)$pred
	math.a15<-predict(math.giv.c.a,c.a.imp)$pred
	c.a.imp$nhpovrt01<-a20	
	read.a20<-predict(read.giv.c.a,c.a.imp)$pred
	math.a20<-predict(math.giv.c.a,c.a.imp)$pred
	c.a.imp$nhpovrt01<-a30	
	read.a30<-predict(read.giv.c.a,c.a.imp)$pred
	math.a30<-predict(math.giv.c.a,c.a.imp)$pred

	# STEP 2: ESTIMATE E(Y(astar,m(a)) #
	read.giv.c.a.m<-ranger(
		y=read.train,
		x=c.a.m.train,
		num.trees=ntrees,
		min.node.size=rf.rd.m2$bestTune[,"min.node.size"],
		mtry=rf.rd.m2$bestTune[,"mtry"],
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=wts.train,
		seed=8675309)

	math.giv.c.a.m<-ranger(
		y=math.train,
		x=c.a.m.train,
		num.trees=ntrees,
		min.node.size=rf.mt.m2$bestTune[,"min.node.size"],
		mtry=rf.mt.m2$bestTune[,"mtry"],
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=wts.train,
		seed=8675309)

	c.a.m.imp<-eclsb[,c(vars.base,vars.covw1,vars.ntxw1)]
	c.a.m.imp$nhpovrt01<-a10	
	read.a10.m<-predict(read.giv.c.a.m,c.a.m.imp)$pred
	math.a10.m<-predict(math.giv.c.a.m,c.a.m.imp)$pred
	c.a.m.imp$nhpovrt01<-a15	
	read.a15.m<-predict(read.giv.c.a.m,c.a.m.imp)$pred
	math.a15.m<-predict(math.giv.c.a.m,c.a.m.imp)$pred
	c.a.m.imp$nhpovrt01<-a20	
	read.a20.m<-predict(read.giv.c.a.m,c.a.m.imp)$pred
	math.a20.m<-predict(math.giv.c.a.m,c.a.m.imp)$pred
	c.a.m.imp$nhpovrt01<-a30	
	read.a30.m<-predict(read.giv.c.a.m,c.a.m.imp)$pred
	math.a30.m<-predict(math.giv.c.a.m,c.a.m.imp)$pred

	read.uhat.a10.m<-ranger(
		y=read.a10.m,
		x=c.a.train,
		num.trees=ntrees,
		min.node.size=rf.rd.m3$bestTune[,"min.node.size"],
		mtry=rf.rd.m3$bestTune[,"mtry"],
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=wts.train,
		seed=8675309)

	math.uhat.a10.m<-ranger(
		y=math.a10.m,
		x=c.a.train,
		num.trees=ntrees,
		min.node.size=rf.mt.m3$bestTune[,"min.node.size"],
		mtry=rf.mt.m3$bestTune[,"mtry"],
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=wts.train,
		seed=8675309)

	read.uhat.a15.m<-ranger(
		y=read.a15.m,
		x=c.a.train,
		num.trees=ntrees,
		min.node.size=rf.rd.m3$bestTune[,"min.node.size"],
		mtry=rf.rd.m3$bestTune[,"mtry"],
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=wts.train,
		seed=8675309)

	math.uhat.a15.m<-ranger(
		y=math.a15.m,
		x=c.a.train,
		num.trees=ntrees,
		min.node.size=rf.mt.m3$bestTune[,"min.node.size"],
		mtry=rf.mt.m3$bestTune[,"mtry"],
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=wts.train,
		seed=8675309)

	read.uhat.a20.m<-ranger(
		y=read.a20.m,
		x=c.a.train,
		num.trees=ntrees,
		min.node.size=rf.rd.m3$bestTune[,"min.node.size"],
		mtry=rf.rd.m3$bestTune[,"mtry"],
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=wts.train,
		seed=8675309)

	math.uhat.a20.m<-ranger(
		y=math.a20.m,
		x=c.a.train,
		num.trees=ntrees,
		min.node.size=rf.mt.m3$bestTune[,"min.node.size"],
		mtry=rf.mt.m3$bestTune[,"mtry"],
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=wts.train,
		seed=8675309)

	read.uhat.a30.m<-ranger(
		y=read.a30.m,
		x=c.a.train,
		num.trees=ntrees,
		min.node.size=rf.rd.m3$bestTune[,"min.node.size"],
		mtry=rf.rd.m3$bestTune[,"mtry"],
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=wts.train,
		seed=8675309)

	math.uhat.a30.m<-ranger(
		y=math.a30.m,
		x=c.a.train,
		num.trees=ntrees,
		min.node.size=rf.mt.m3$bestTune[,"min.node.size"],
		mtry=rf.mt.m3$bestTune[,"mtry"],
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=wts.train,
		seed=8675309)

	c.a.imp<-eclsb[,c(vars.base,vars.covw1)]
	c.a.imp$nhpovrt01<-a5	
	read.a10.mofa<-predict(read.uhat.a10.m,c.a.imp)$pred
	math.a10.mofa<-predict(math.uhat.a10.m,c.a.imp)$pred
	read.a15.mofa<-predict(read.uhat.a15.m,c.a.imp)$pred
	math.a15.mofa<-predict(math.uhat.a15.m,c.a.imp)$pred
	read.a20.mofa<-predict(read.uhat.a20.m,c.a.imp)$pred
	math.a20.mofa<-predict(math.uhat.a20.m,c.a.imp)$pred
	read.a30.mofa<-predict(read.uhat.a30.m,c.a.imp)$pred
	math.a30.mofa<-predict(math.uhat.a30.m,c.a.imp)$pred

	# STEP 3: COMPUTE ATE, NDE, NIE #
	miest.rd.10v5[i,1]<-weighted.mean(read.a10,wts.train)-weighted.mean(read.a5,wts.train)
	miest.rd.10v5[i,2]<-weighted.mean(read.a10.mofa,wts.train)-weighted.mean(read.a5,wts.train)
	miest.rd.10v5[i,3]<-weighted.mean(read.a10,wts.train)-weighted.mean(read.a10.mofa,wts.train)

	miest.mt.10v5[i,1]<-weighted.mean(math.a10,wts.train)-weighted.mean(math.a5,wts.train)
	miest.mt.10v5[i,2]<-weighted.mean(math.a10.mofa,wts.train)-weighted.mean(math.a5,wts.train)
	miest.mt.10v5[i,3]<-weighted.mean(math.a10,wts.train)-weighted.mean(math.a10.mofa,wts.train)

	miest.rd.15v5[i,1]<-weighted.mean(read.a15,wts.train)-weighted.mean(read.a5,wts.train)
	miest.rd.15v5[i,2]<-weighted.mean(read.a15.mofa,wts.train)-weighted.mean(read.a5,wts.train)
	miest.rd.15v5[i,3]<-weighted.mean(read.a15,wts.train)-weighted.mean(read.a15.mofa,wts.train)

	miest.mt.15v5[i,1]<-weighted.mean(math.a15,wts.train)-weighted.mean(math.a5,wts.train)
	miest.mt.15v5[i,2]<-weighted.mean(math.a15.mofa,wts.train)-weighted.mean(math.a5,wts.train)
	miest.mt.15v5[i,3]<-weighted.mean(math.a15,wts.train)-weighted.mean(math.a15.mofa,wts.train)

	miest.rd.20v5[i,1]<-weighted.mean(read.a20,wts.train)-weighted.mean(read.a5,wts.train)
	miest.rd.20v5[i,2]<-weighted.mean(read.a20.mofa,wts.train)-weighted.mean(read.a5,wts.train)
	miest.rd.20v5[i,3]<-weighted.mean(read.a20,wts.train)-weighted.mean(read.a20.mofa,wts.train)

	miest.mt.20v5[i,1]<-weighted.mean(math.a20,wts.train)-weighted.mean(math.a5,wts.train)
	miest.mt.20v5[i,2]<-weighted.mean(math.a20.mofa,wts.train)-weighted.mean(math.a5,wts.train)
	miest.mt.20v5[i,3]<-weighted.mean(math.a20,wts.train)-weighted.mean(math.a20.mofa,wts.train)

	miest.rd.30v5[i,1]<-weighted.mean(read.a30,wts.train)-weighted.mean(read.a5,wts.train)
	miest.rd.30v5[i,2]<-weighted.mean(read.a30.mofa,wts.train)-weighted.mean(read.a5,wts.train)
	miest.rd.30v5[i,3]<-weighted.mean(read.a30,wts.train)-weighted.mean(read.a30.mofa,wts.train)

	miest.mt.30v5[i,1]<-weighted.mean(math.a30,wts.train)-weighted.mean(math.a5,wts.train)
	miest.mt.30v5[i,2]<-weighted.mean(math.a30.mofa,wts.train)-weighted.mean(math.a5,wts.train)
	miest.mt.30v5[i,3]<-weighted.mean(math.a30,wts.train)-weighted.mean(math.a30.mofa,wts.train)

	### COMPUTE BOOTSTRAP CIs ###
	my.cluster<-parallel::makeCluster(n.cores,type="PSOCK")
	doParallel::registerDoParallel(cl=my.cluster)
	foreach::getDoParRegistered()
	clusterEvalQ(cl=my.cluster, library(ranger))
	registerDoRNG(8675309)

	clusterExport(cl=my.cluster,
		list(
			"vars.ntxw1",
			"vars.covw1",
			"vars.base",
			"nmi",
			"nboot",
			"a5",
			"a10",
			"a15",
			"a20",
			"a30",
			"eclsb",
			"ntrees",
			"rf.rd.m1",
			"rf.mt.m1",
			"rf.rd.m2",
			"rf.mt.m2",
			"rf.rd.m3",
			"rf.mt.m3"),
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

		boot.read.giv.c.a<-ranger(
			y=boot.read.train,
			x=boot.c.a.train,
			num.trees=ntrees,
			min.node.size=rf.rd.m1$bestTune[,"min.node.size"],
			mtry=rf.rd.m1$bestTune[,"mtry"],
			splitrule="variance",
			respect.unordered.factors=TRUE,
			case.weights=boot.wts.train,
			seed=8675309)

		boot.math.giv.c.a<-ranger(
			y=boot.math.train,
			x=boot.c.a.train,
			num.trees=ntrees,
			min.node.size=rf.mt.m1$bestTune[,"min.node.size"],
			mtry=rf.mt.m1$bestTune[,"mtry"],
			splitrule="variance",
			respect.unordered.factors=TRUE,
			case.weights=boot.wts.train,
			seed=8675309)

		boot.c.a.imp<-boot.eclsb[,c(vars.base,vars.covw1)]
		boot.c.a.imp$nhpovrt01<-a5
		boot.read.a5<-predict(boot.read.giv.c.a,boot.c.a.imp)$pred
		boot.math.a5<-predict(boot.math.giv.c.a,boot.c.a.imp)$pred
		boot.c.a.imp$nhpovrt01<-a10	
		boot.read.a10<-predict(boot.read.giv.c.a,boot.c.a.imp)$pred
		boot.math.a10<-predict(boot.math.giv.c.a,boot.c.a.imp)$pred
		boot.c.a.imp$nhpovrt01<-a15	
		boot.read.a15<-predict(boot.read.giv.c.a,boot.c.a.imp)$pred
		boot.math.a15<-predict(boot.math.giv.c.a,boot.c.a.imp)$pred
		boot.c.a.imp$nhpovrt01<-a20	
		boot.read.a20<-predict(boot.read.giv.c.a,boot.c.a.imp)$pred
		boot.math.a20<-predict(boot.math.giv.c.a,boot.c.a.imp)$pred
		boot.c.a.imp$nhpovrt01<-a30	
		boot.read.a30<-predict(boot.read.giv.c.a,boot.c.a.imp)$pred
		boot.math.a30<-predict(boot.math.giv.c.a,boot.c.a.imp)$pred

		boot.read.giv.c.a.m<-ranger(
			y=boot.read.train,
			x=boot.c.a.m.train,
			num.trees=ntrees,
			min.node.size=rf.rd.m2$bestTune[,"min.node.size"],
			mtry=rf.rd.m2$bestTune[,"mtry"],
			splitrule="variance",
			respect.unordered.factors=TRUE,
			case.weights=boot.wts.train,
			seed=8675309)

		boot.math.giv.c.a.m<-ranger(
			y=boot.math.train,
			x=boot.c.a.m.train,
			num.trees=ntrees,
			min.node.size=rf.mt.m2$bestTune[,"min.node.size"],
			mtry=rf.mt.m2$bestTune[,"mtry"],
			splitrule="variance",
			respect.unordered.factors=TRUE,
			case.weights=boot.wts.train,
			seed=8675309)

		boot.c.a.m.imp<-boot.eclsb[,c(vars.base,vars.covw1,vars.ntxw1)]
		boot.c.a.m.imp$nhpovrt01<-a10	
		boot.read.a10.m<-predict(boot.read.giv.c.a.m,boot.c.a.m.imp)$pred
		boot.math.a10.m<-predict(boot.math.giv.c.a.m,boot.c.a.m.imp)$pred
		boot.c.a.m.imp$nhpovrt01<-a15	
		boot.read.a15.m<-predict(boot.read.giv.c.a.m,boot.c.a.m.imp)$pred
		boot.math.a15.m<-predict(boot.math.giv.c.a.m,boot.c.a.m.imp)$pred
		boot.c.a.m.imp$nhpovrt01<-a20	
		boot.read.a20.m<-predict(boot.read.giv.c.a.m,boot.c.a.m.imp)$pred
		boot.math.a20.m<-predict(boot.math.giv.c.a.m,boot.c.a.m.imp)$pred
		boot.c.a.m.imp$nhpovrt01<-a30	
		boot.read.a30.m<-predict(boot.read.giv.c.a.m,boot.c.a.m.imp)$pred
		boot.math.a30.m<-predict(boot.math.giv.c.a.m,boot.c.a.m.imp)$pred

		boot.read.uhat.a10.m<-ranger(
			y=boot.read.a10.m,
			x=boot.c.a.train,
			num.trees=ntrees,
			min.node.size=rf.rd.m3$bestTune[,"min.node.size"],
			mtry=rf.rd.m3$bestTune[,"mtry"],
			splitrule="variance",
			respect.unordered.factors=TRUE,
			case.weights=boot.wts.train,
			seed=8675309)

		boot.math.uhat.a10.m<-ranger(
			y=boot.math.a10.m,
			x=boot.c.a.train,
			num.trees=ntrees,
			min.node.size=rf.mt.m3$bestTune[,"min.node.size"],
			mtry=rf.mt.m3$bestTune[,"mtry"],
			splitrule="variance",
			respect.unordered.factors=TRUE,
			case.weights=boot.wts.train,
			seed=8675309)

		boot.read.uhat.a15.m<-ranger(
			y=boot.read.a15.m,
			x=boot.c.a.train,
			num.trees=ntrees,
			min.node.size=rf.rd.m3$bestTune[,"min.node.size"],
			mtry=rf.rd.m3$bestTune[,"mtry"],
			splitrule="variance",
			respect.unordered.factors=TRUE,
			case.weights=boot.wts.train,
			seed=8675309)

		boot.math.uhat.a15.m<-ranger(
			y=boot.math.a15.m,
			x=boot.c.a.train,
			num.trees=ntrees,
			min.node.size=rf.mt.m3$bestTune[,"min.node.size"],
			mtry=rf.mt.m3$bestTune[,"mtry"],
			splitrule="variance",
			respect.unordered.factors=TRUE,
			case.weights=boot.wts.train,
			seed=8675309)

		boot.read.uhat.a20.m<-ranger(
			y=boot.read.a20.m,
			x=boot.c.a.train,
			num.trees=ntrees,
			min.node.size=rf.rd.m3$bestTune[,"min.node.size"],
			mtry=rf.rd.m3$bestTune[,"mtry"],
			splitrule="variance",
			respect.unordered.factors=TRUE,
			case.weights=boot.wts.train,
			seed=8675309)

		boot.math.uhat.a20.m<-ranger(
			y=boot.math.a20.m,
			x=boot.c.a.train,
			num.trees=ntrees,
			min.node.size=rf.mt.m3$bestTune[,"min.node.size"],
			mtry=rf.mt.m3$bestTune[,"mtry"],
			splitrule="variance",
			respect.unordered.factors=TRUE,
			case.weights=boot.wts.train,
			seed=8675309)

		boot.read.uhat.a30.m<-ranger(
			y=boot.read.a30.m,
			x=boot.c.a.train,
			num.trees=ntrees,
			min.node.size=rf.rd.m3$bestTune[,"min.node.size"],
			mtry=rf.rd.m3$bestTune[,"mtry"],
			splitrule="variance",
			respect.unordered.factors=TRUE,
			case.weights=boot.wts.train,
			seed=8675309)

		boot.math.uhat.a30.m<-ranger(
			y=boot.math.a30.m,
			x=boot.c.a.train,
			num.trees=ntrees,
			min.node.size=rf.mt.m3$bestTune[,"min.node.size"],
			mtry=rf.mt.m3$bestTune[,"mtry"],
			splitrule="variance",
			respect.unordered.factors=TRUE,
			case.weights=boot.wts.train,
			seed=8675309)

		boot.c.a.imp<-boot.eclsb[,c(vars.base,vars.covw1)]
		boot.c.a.imp$nhpovrt01<-a5	
		boot.read.a10.mofa<-predict(boot.read.uhat.a10.m,boot.c.a.imp)$pred
		boot.math.a10.mofa<-predict(boot.math.uhat.a10.m,boot.c.a.imp)$pred
		boot.read.a15.mofa<-predict(boot.read.uhat.a15.m,boot.c.a.imp)$pred
		boot.math.a15.mofa<-predict(boot.math.uhat.a15.m,boot.c.a.imp)$pred
		boot.read.a20.mofa<-predict(boot.read.uhat.a20.m,boot.c.a.imp)$pred
		boot.math.a20.mofa<-predict(boot.math.uhat.a20.m,boot.c.a.imp)$pred
		boot.read.a30.mofa<-predict(boot.read.uhat.a30.m,boot.c.a.imp)$pred
		boot.math.a30.mofa<-predict(boot.math.uhat.a30.m,boot.c.a.imp)$pred

		boot.ate.rd.10v5<-weighted.mean(boot.read.a10,boot.wts.train)-weighted.mean(boot.read.a5,boot.wts.train)
		boot.nde.rd.10v5<-weighted.mean(boot.read.a10.mofa,boot.wts.train)-weighted.mean(boot.read.a5,boot.wts.train)
		boot.nie.rd.10v5<-weighted.mean(boot.read.a10,boot.wts.train)-weighted.mean(boot.read.a10.mofa,boot.wts.train)

		boot.ate.mt.10v5<-weighted.mean(boot.math.a10,boot.wts.train)-weighted.mean(boot.math.a5,boot.wts.train)
		boot.nde.mt.10v5<-weighted.mean(boot.math.a10.mofa,boot.wts.train)-weighted.mean(boot.math.a5,boot.wts.train)
		boot.nie.mt.10v5<-weighted.mean(boot.math.a10,boot.wts.train)-weighted.mean(boot.math.a10.mofa,boot.wts.train)

		boot.ate.rd.15v5<-weighted.mean(boot.read.a15,boot.wts.train)-weighted.mean(boot.read.a5,boot.wts.train)
		boot.nde.rd.15v5<-weighted.mean(boot.read.a15.mofa,boot.wts.train)-weighted.mean(boot.read.a5,boot.wts.train)
		boot.nie.rd.15v5<-weighted.mean(boot.read.a15,boot.wts.train)-weighted.mean(boot.read.a15.mofa,boot.wts.train)

		boot.ate.mt.15v5<-weighted.mean(boot.math.a15,boot.wts.train)-weighted.mean(boot.math.a5,boot.wts.train)
		boot.nde.mt.15v5<-weighted.mean(boot.math.a15.mofa,boot.wts.train)-weighted.mean(boot.math.a5,boot.wts.train)
		boot.nie.mt.15v5<-weighted.mean(boot.math.a15,boot.wts.train)-weighted.mean(boot.math.a15.mofa,boot.wts.train)

		boot.ate.rd.20v5<-weighted.mean(boot.read.a20,boot.wts.train)-weighted.mean(boot.read.a5,boot.wts.train)
		boot.nde.rd.20v5<-weighted.mean(boot.read.a20.mofa,boot.wts.train)-weighted.mean(boot.read.a5,boot.wts.train)
		boot.nie.rd.20v5<-weighted.mean(boot.read.a20,boot.wts.train)-weighted.mean(boot.read.a20.mofa,boot.wts.train)

		boot.ate.mt.20v5<-weighted.mean(boot.math.a20,boot.wts.train)-weighted.mean(boot.math.a5,boot.wts.train)
		boot.nde.mt.20v5<-weighted.mean(boot.math.a20.mofa,boot.wts.train)-weighted.mean(boot.math.a5,boot.wts.train)
		boot.nie.mt.20v5<-weighted.mean(boot.math.a20,boot.wts.train)-weighted.mean(boot.math.a20.mofa,boot.wts.train)

		boot.ate.rd.30v5<-weighted.mean(boot.read.a30,boot.wts.train)-weighted.mean(boot.read.a5,boot.wts.train)
		boot.nde.rd.30v5<-weighted.mean(boot.read.a30.mofa,boot.wts.train)-weighted.mean(boot.read.a5,boot.wts.train)
		boot.nie.rd.30v5<-weighted.mean(boot.read.a30,boot.wts.train)-weighted.mean(boot.read.a30.mofa,boot.wts.train)

		boot.ate.mt.30v5<-weighted.mean(boot.math.a30,boot.wts.train)-weighted.mean(boot.math.a5,boot.wts.train)
		boot.nde.mt.30v5<-weighted.mean(boot.math.a30.mofa,boot.wts.train)-weighted.mean(boot.math.a5,boot.wts.train)
		boot.nie.mt.30v5<-weighted.mean(boot.math.a30,boot.wts.train)-weighted.mean(boot.math.a30.mofa,boot.wts.train)
		
		return(
			list(
				boot.ate.rd.10v5,boot.nde.rd.10v5,boot.nie.rd.10v5,
				boot.ate.mt.10v5,boot.nde.mt.10v5,boot.nie.mt.10v5,
				boot.ate.rd.15v5,boot.nde.rd.15v5,boot.nie.rd.15v5,
				boot.ate.mt.15v5,boot.nde.mt.15v5,boot.nie.mt.15v5,
				boot.ate.rd.20v5,boot.nde.rd.20v5,boot.nie.rd.20v5,
				boot.ate.mt.20v5,boot.nde.mt.20v5,boot.nie.mt.20v5,
				boot.ate.rd.30v5,boot.nde.rd.30v5,boot.nie.rd.30v5,
				boot.ate.mt.30v5,boot.nde.mt.30v5,boot.nie.mt.30v5))
		}

	boot.est.w1<-matrix(unlist(boot.est.w1),ncol=24,byrow=TRUE)

	boot.ci.rd.10v5<-rbind(boot.ci.rd.10v5,boot.est.w1[,1:3])
	boot.ci.mt.10v5<-rbind(boot.ci.mt.10v5,boot.est.w1[,4:6])

	boot.ci.rd.15v5<-rbind(boot.ci.rd.15v5,boot.est.w1[,7:9])
	boot.ci.mt.15v5<-rbind(boot.ci.mt.15v5,boot.est.w1[,10:12])

	boot.ci.rd.20v5<-rbind(boot.ci.rd.20v5,boot.est.w1[,13:15])
	boot.ci.mt.20v5<-rbind(boot.ci.mt.20v5,boot.est.w1[,16:18])

	boot.ci.rd.30v5<-rbind(boot.ci.rd.30v5,boot.est.w1[,19:21])
	boot.ci.mt.30v5<-rbind(boot.ci.mt.30v5,boot.est.w1[,22:24])

	stopCluster(my.cluster)
	rm(my.cluster)
	}

### COMBINE MI ESTIMATES ###
est.rd.10v5<-est.mt.10v5<-matrix(data=NA,nrow=3,ncol=3)
est.rd.15v5<-est.mt.15v5<-matrix(data=NA,nrow=3,ncol=3)
est.rd.20v5<-est.mt.20v5<-matrix(data=NA,nrow=3,ncol=3)
est.rd.30v5<-est.mt.30v5<-matrix(data=NA,nrow=3,ncol=3)

for (i in 1:3) { 
	
	est.rd.10v5[i,1]<-round(mean(miest.rd.10v5[,i]),digits=3)
	est.rd.10v5[i,2]<-round(quantile(boot.ci.rd.10v5[,i],prob=0.025),digits=3)
	est.rd.10v5[i,3]<-round(quantile(boot.ci.rd.10v5[,i],prob=0.975),digits=3)

	est.mt.10v5[i,1]<-round(mean(miest.mt.10v5[,i]),digits=3)
	est.mt.10v5[i,2]<-round(quantile(boot.ci.mt.10v5[,i],prob=0.025),digits=3)
	est.mt.10v5[i,3]<-round(quantile(boot.ci.mt.10v5[,i],prob=0.975),digits=3)

	est.rd.15v5[i,1]<-round(mean(miest.rd.15v5[,i]),digits=3)
	est.rd.15v5[i,2]<-round(quantile(boot.ci.rd.15v5[,i],prob=0.025),digits=3)
	est.rd.15v5[i,3]<-round(quantile(boot.ci.rd.15v5[,i],prob=0.975),digits=3)

	est.mt.15v5[i,1]<-round(mean(miest.mt.15v5[,i]),digits=3)
	est.mt.15v5[i,2]<-round(quantile(boot.ci.mt.15v5[,i],prob=0.025),digits=3)
	est.mt.15v5[i,3]<-round(quantile(boot.ci.mt.15v5[,i],prob=0.975),digits=3)

	est.rd.20v5[i,1]<-round(mean(miest.rd.20v5[,i]),digits=3)
	est.rd.20v5[i,2]<-round(quantile(boot.ci.rd.20v5[,i],prob=0.025),digits=3)
	est.rd.20v5[i,3]<-round(quantile(boot.ci.rd.20v5[,i],prob=0.975),digits=3)

	est.mt.20v5[i,1]<-round(mean(miest.mt.20v5[,i]),digits=3)
	est.mt.20v5[i,2]<-round(quantile(boot.ci.mt.20v5[,i],prob=0.025),digits=3)
	est.mt.20v5[i,3]<-round(quantile(boot.ci.mt.20v5[,i],prob=0.975),digits=3)

	est.rd.30v5[i,1]<-round(mean(miest.rd.30v5[,i]),digits=3)
	est.rd.30v5[i,2]<-round(quantile(boot.ci.rd.30v5[,i],prob=0.025),digits=3)
	est.rd.30v5[i,3]<-round(quantile(boot.ci.rd.30v5[,i],prob=0.975),digits=3)

	est.mt.30v5[i,1]<-round(mean(miest.mt.30v5[,i]),digits=3)
	est.mt.30v5[i,2]<-round(quantile(boot.ci.mt.30v5[,i],prob=0.025),digits=3)
	est.mt.30v5[i,3]<-round(quantile(boot.ci.mt.30v5[,i],prob=0.975),digits=3)
	}

rlabel<-c('ATE','NDE','NIE')
output.rd.10v5<-data.frame(est.rd.10v5,row.names=rlabel)
output.mt.10v5<-data.frame(est.mt.10v5,row.names=rlabel)
output.rd.15v5<-data.frame(est.rd.15v5,row.names=rlabel)
output.mt.15v5<-data.frame(est.mt.15v5,row.names=rlabel)
output.rd.20v5<-data.frame(est.rd.20v5,row.names=rlabel)
output.mt.20v5<-data.frame(est.mt.20v5,row.names=rlabel)
output.rd.30v5<-data.frame(est.rd.30v5,row.names=rlabel)
output.mt.30v5<-data.frame(est.mt.30v5,row.names=rlabel)

colnames(output.rd.10v5)<-colnames(output.mt.10v5)<-c('estimate','ll.pct.95ci','ul.pct.95ci')
colnames(output.rd.15v5)<-colnames(output.mt.15v5)<-c('estimate','ll.pct.95ci','ul.pct.95ci')
colnames(output.rd.20v5)<-colnames(output.mt.20v5)<-c('estimate','ll.pct.95ci','ul.pct.95ci')
colnames(output.rd.30v5)<-colnames(output.mt.30v5)<-c('estimate','ll.pct.95ci','ul.pct.95ci')

### PRINT RESULTS ###
sink("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\programs\\_LOGS\\24_create_table_S4_log.txt")

cat("===========================================\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Reading Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("*** 10% vs. 5% Poverty ***\n")
print(output.rd.10v5)
cat(" \n")
cat("*** 15% vs. 5% Poverty ***\n")
print(output.rd.15v5)
cat(" \n")
cat("*** 20% vs. 5% Poverty ***\n")
print(output.rd.20v5)
cat(" \n")
cat("*** 30% vs. 5% Poverty ***\n")
print(output.rd.30v5)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Math Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("*** 10% vs. 5% Poverty ***\n")
print(output.mt.10v5)
cat(" \n")
cat("*** 15% vs. 5% Poverty ***\n")
print(output.mt.15v5)
cat(" \n")
cat("*** 20% vs. 5% Poverty ***\n")
print(output.mt.20v5)
cat(" \n")
cat("*** 30% vs. 5% Poverty ***\n")
print(output.mt.30v5)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("===========================================\n")

print(startTime)
print(Sys.time())

sink()

