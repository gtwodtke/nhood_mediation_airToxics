################################################
################################################
##                                            ##
## PROGRAM NAME: 27_create_table_S7           ##
## AUTHOR: GW                                 ##
## DESCRIPTION:                               ##
##                                            ##
##  create table of nh effect estimates by    ##
##  region using random forests and reg       ##
##  imputation                                ## 
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
astar<-0.25
a<-0.05

##### DEFINE COVARIATE SETS #####
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

##### STANDARDIZE ######
for (v in 1:length(vars.ntxw1)) {
	eclsb.mi[,vars.ntxw1[v]]<-(eclsb.mi[,vars.ntxw1[v]]-weighted.mean(eclsb.mi[,vars.ntxw1[v]],eclsb.mi$sampwt))/sqrt(wtd.var(eclsb.mi[,vars.ntxw1[v]],eclsb.mi$sampwt))
	}

eclsb.mi$readtheta05<-(eclsb.mi$readtheta05-weighted.mean(eclsb.mi$readtheta05,eclsb.mi$sampwt))/sqrt(wtd.var(eclsb.mi$readtheta05,eclsb.mi$sampwt))
eclsb.mi$maththeta05<-(eclsb.mi$maththeta05-weighted.mean(eclsb.mi$maththeta05,eclsb.mi$sampwt))/sqrt(wtd.var(eclsb.mi$maththeta05,eclsb.mi$sampwt))

##### ESTIMATE NHOOD EFFECTS #####
miest.rd.ne<-miest.mt.ne<-matrix(data=NA,nrow=nmi,ncol=3)
miest.rd.mw<-miest.mt.mw<-matrix(data=NA,nrow=nmi,ncol=3)
miest.rd.so<-miest.mt.so<-matrix(data=NA,nrow=nmi,ncol=3)
miest.rd.we<-miest.mt.we<-matrix(data=NA,nrow=nmi,ncol=3)

boot.ci.rd.ne<-boot.ci.mt.ne<-NULL
boot.ci.rd.mw<-boot.ci.mt.mw<-NULL
boot.ci.rd.so<-boot.ci.mt.so<-NULL
boot.ci.rd.we<-boot.ci.mt.we<-NULL

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
	eclsb.tune.imp$nhpovrt01<-astar	
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
	c.a.imp$nhpovrt01<-astar	
	eclsb$read.astar<-predict(read.giv.c.a,c.a.imp)$pred
	eclsb$math.astar<-predict(math.giv.c.a,c.a.imp)$pred
	c.a.imp$nhpovrt01<-a	
	eclsb$read.a<-predict(read.giv.c.a,c.a.imp)$pred
	eclsb$math.a<-predict(math.giv.c.a,c.a.imp)$pred

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
	c.a.m.imp$nhpovrt01<-astar	
	read.astar.m<-predict(read.giv.c.a.m,c.a.m.imp)$pred
	math.astar.m<-predict(math.giv.c.a.m,c.a.m.imp)$pred
		
	read.uhat.astar.m<-ranger(
		y=read.astar.m,
		x=c.a.train,
		num.trees=ntrees,
		min.node.size=rf.rd.m3$bestTune[,"min.node.size"],
		mtry=rf.rd.m3$bestTune[,"mtry"],
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=wts.train,
		seed=8675309)

	math.uhat.astar.m<-ranger(
		y=math.astar.m,
		x=c.a.train,
		num.trees=ntrees,
		min.node.size=rf.mt.m3$bestTune[,"min.node.size"],
		mtry=rf.mt.m3$bestTune[,"mtry"],
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=wts.train,
		seed=8675309)

	c.a.imp<-eclsb[,c(vars.base,vars.covw1)]
	c.a.imp$nhpovrt01<-a	
	eclsb$read.astar.mofa<-predict(read.uhat.astar.m,c.a.imp)$pred
	eclsb$math.astar.mofa<-predict(math.uhat.astar.m,c.a.imp)$pred

	# STEP 3: COMPUTE ATE, NDE, NIE #
	miest.rd.ne[i,1]<-weighted.mean(eclsb$read.astar[eclsb$region01=="Northeast"],eclsb$sampwt[eclsb$region01=="Northeast"])-weighted.mean(eclsb$read.a[eclsb$region01=="Northeast"],eclsb$sampwt[eclsb$region01=="Northeast"])
	miest.rd.ne[i,2]<-weighted.mean(eclsb$read.astar.mofa[eclsb$region01=="Northeast"],eclsb$sampwt[eclsb$region01=="Northeast"])-weighted.mean(eclsb$read.a[eclsb$region01=="Northeast"],eclsb$sampwt[eclsb$region01=="Northeast"])
	miest.rd.ne[i,3]<-weighted.mean(eclsb$read.astar[eclsb$region01=="Northeast"],eclsb$sampwt[eclsb$region01=="Northeast"])-weighted.mean(eclsb$read.astar.mofa[eclsb$region01=="Northeast"],eclsb$sampwt[eclsb$region01=="Northeast"])

	miest.mt.ne[i,1]<-weighted.mean(eclsb$math.astar[eclsb$region01=="Northeast"],eclsb$sampwt[eclsb$region01=="Northeast"])-weighted.mean(eclsb$math.a[eclsb$region01=="Northeast"],eclsb$sampwt[eclsb$region01=="Northeast"])
	miest.mt.ne[i,2]<-weighted.mean(eclsb$math.astar.mofa[eclsb$region01=="Northeast"],eclsb$sampwt[eclsb$region01=="Northeast"])-weighted.mean(eclsb$math.a[eclsb$region01=="Northeast"],eclsb$sampwt[eclsb$region01=="Northeast"])
	miest.mt.ne[i,3]<-weighted.mean(eclsb$math.astar[eclsb$region01=="Northeast"],eclsb$sampwt[eclsb$region01=="Northeast"])-weighted.mean(eclsb$math.astar.mofa[eclsb$region01=="Northeast"],eclsb$sampwt[eclsb$region01=="Northeast"])

	miest.rd.mw[i,1]<-weighted.mean(eclsb$read.astar[eclsb$region01=="Midwest"],eclsb$sampwt[eclsb$region01=="Midwest"])-weighted.mean(eclsb$read.a[eclsb$region01=="Midwest"],eclsb$sampwt[eclsb$region01=="Midwest"])
	miest.rd.mw[i,2]<-weighted.mean(eclsb$read.astar.mofa[eclsb$region01=="Midwest"],eclsb$sampwt[eclsb$region01=="Midwest"])-weighted.mean(eclsb$read.a[eclsb$region01=="Midwest"],eclsb$sampwt[eclsb$region01=="Midwest"])
	miest.rd.mw[i,3]<-weighted.mean(eclsb$read.astar[eclsb$region01=="Midwest"],eclsb$sampwt[eclsb$region01=="Midwest"])-weighted.mean(eclsb$read.astar.mofa[eclsb$region01=="Midwest"],eclsb$sampwt[eclsb$region01=="Midwest"])

	miest.mt.mw[i,1]<-weighted.mean(eclsb$math.astar[eclsb$region01=="Midwest"],eclsb$sampwt[eclsb$region01=="Midwest"])-weighted.mean(eclsb$math.a[eclsb$region01=="Midwest"],eclsb$sampwt[eclsb$region01=="Midwest"])
	miest.mt.mw[i,2]<-weighted.mean(eclsb$math.astar.mofa[eclsb$region01=="Midwest"],eclsb$sampwt[eclsb$region01=="Midwest"])-weighted.mean(eclsb$math.a[eclsb$region01=="Midwest"],eclsb$sampwt[eclsb$region01=="Midwest"])
	miest.mt.mw[i,3]<-weighted.mean(eclsb$math.astar[eclsb$region01=="Midwest"],eclsb$sampwt[eclsb$region01=="Midwest"])-weighted.mean(eclsb$math.astar.mofa[eclsb$region01=="Midwest"],eclsb$sampwt[eclsb$region01=="Midwest"])

	miest.rd.so[i,1]<-weighted.mean(eclsb$read.astar[eclsb$region01=="South"],eclsb$sampwt[eclsb$region01=="South"])-weighted.mean(eclsb$read.a[eclsb$region01=="South"],eclsb$sampwt[eclsb$region01=="South"])
	miest.rd.so[i,2]<-weighted.mean(eclsb$read.astar.mofa[eclsb$region01=="South"],eclsb$sampwt[eclsb$region01=="South"])-weighted.mean(eclsb$read.a[eclsb$region01=="South"],eclsb$sampwt[eclsb$region01=="South"])
	miest.rd.so[i,3]<-weighted.mean(eclsb$read.astar[eclsb$region01=="South"],eclsb$sampwt[eclsb$region01=="South"])-weighted.mean(eclsb$read.astar.mofa[eclsb$region01=="South"],eclsb$sampwt[eclsb$region01=="South"])

	miest.mt.so[i,1]<-weighted.mean(eclsb$math.astar[eclsb$region01=="South"],eclsb$sampwt[eclsb$region01=="South"])-weighted.mean(eclsb$math.a[eclsb$region01=="South"],eclsb$sampwt[eclsb$region01=="South"])
	miest.mt.so[i,2]<-weighted.mean(eclsb$math.astar.mofa[eclsb$region01=="South"],eclsb$sampwt[eclsb$region01=="South"])-weighted.mean(eclsb$math.a[eclsb$region01=="South"],eclsb$sampwt[eclsb$region01=="South"])
	miest.mt.so[i,3]<-weighted.mean(eclsb$math.astar[eclsb$region01=="South"],eclsb$sampwt[eclsb$region01=="South"])-weighted.mean(eclsb$math.astar.mofa[eclsb$region01=="South"],eclsb$sampwt[eclsb$region01=="South"])

	miest.rd.we[i,1]<-weighted.mean(eclsb$read.astar[eclsb$region01=="West"],eclsb$sampwt[eclsb$region01=="West"])-weighted.mean(eclsb$read.a[eclsb$region01=="West"],eclsb$sampwt[eclsb$region01=="West"])
	miest.rd.we[i,2]<-weighted.mean(eclsb$read.astar.mofa[eclsb$region01=="West"],eclsb$sampwt[eclsb$region01=="West"])-weighted.mean(eclsb$read.a[eclsb$region01=="West"],eclsb$sampwt[eclsb$region01=="West"])
	miest.rd.we[i,3]<-weighted.mean(eclsb$read.astar[eclsb$region01=="West"],eclsb$sampwt[eclsb$region01=="West"])-weighted.mean(eclsb$read.astar.mofa[eclsb$region01=="West"],eclsb$sampwt[eclsb$region01=="West"])

	miest.mt.we[i,1]<-weighted.mean(eclsb$math.astar[eclsb$region01=="West"],eclsb$sampwt[eclsb$region01=="West"])-weighted.mean(eclsb$math.a[eclsb$region01=="West"],eclsb$sampwt[eclsb$region01=="West"])
	miest.mt.we[i,2]<-weighted.mean(eclsb$math.astar.mofa[eclsb$region01=="West"],eclsb$sampwt[eclsb$region01=="West"])-weighted.mean(eclsb$math.a[eclsb$region01=="West"],eclsb$sampwt[eclsb$region01=="West"])
	miest.mt.we[i,3]<-weighted.mean(eclsb$math.astar[eclsb$region01=="West"],eclsb$sampwt[eclsb$region01=="West"])-weighted.mean(eclsb$math.astar.mofa[eclsb$region01=="West"],eclsb$sampwt[eclsb$region01=="West"])

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
			"astar",
			"a",
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
		boot.c.a.imp$nhpovrt01<-astar	
		boot.eclsb$read.astar<-predict(boot.read.giv.c.a,boot.c.a.imp)$pred
		boot.eclsb$math.astar<-predict(boot.math.giv.c.a,boot.c.a.imp)$pred
		boot.c.a.imp$nhpovrt01<-a	
		boot.eclsb$read.a<-predict(boot.read.giv.c.a,boot.c.a.imp)$pred
		boot.eclsb$math.a<-predict(boot.math.giv.c.a,boot.c.a.imp)$pred

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
		boot.c.a.m.imp$nhpovrt01<-astar	
		boot.read.astar.m<-predict(boot.read.giv.c.a.m,boot.c.a.m.imp)$pred
		boot.math.astar.m<-predict(boot.math.giv.c.a.m,boot.c.a.m.imp)$pred

		boot.read.uhat.astar.m<-ranger(
			y=boot.read.astar.m,
			x=boot.c.a.train,
			num.trees=ntrees,
			min.node.size=rf.rd.m3$bestTune[,"min.node.size"],
			mtry=rf.rd.m3$bestTune[,"mtry"],
			splitrule="variance",
			respect.unordered.factors=TRUE,
			case.weights=boot.wts.train,
			seed=8675309)

		boot.math.uhat.astar.m<-ranger(
			y=boot.math.astar.m,
			x=boot.c.a.train,
			num.trees=ntrees,
			min.node.size=rf.mt.m3$bestTune[,"min.node.size"],
			mtry=rf.mt.m3$bestTune[,"mtry"],
			splitrule="variance",
			respect.unordered.factors=TRUE,
			case.weights=boot.wts.train,
			seed=8675309)

		boot.c.a.imp<-boot.eclsb[,c(vars.base,vars.covw1)]
		boot.c.a.imp$nhpovrt01<-a	
		boot.eclsb$read.astar.mofa<-predict(boot.read.uhat.astar.m,boot.c.a.imp)$pred
		boot.eclsb$math.astar.mofa<-predict(boot.math.uhat.astar.m,boot.c.a.imp)$pred

		boot.ate.rd.ne<-weighted.mean(boot.eclsb$read.astar[boot.eclsb$region01=="Northeast"],boot.eclsb$sampwt[boot.eclsb$region01=="Northeast"])-weighted.mean(boot.eclsb$read.a[boot.eclsb$region01=="Northeast"],boot.eclsb$sampwt[boot.eclsb$region01=="Northeast"])
		boot.nde.rd.ne<-weighted.mean(boot.eclsb$read.astar.mofa[boot.eclsb$region01=="Northeast"],boot.eclsb$sampwt[boot.eclsb$region01=="Northeast"])-weighted.mean(boot.eclsb$read.a[boot.eclsb$region01=="Northeast"],boot.eclsb$sampwt[boot.eclsb$region01=="Northeast"])
		boot.nie.rd.ne<-weighted.mean(boot.eclsb$read.astar[boot.eclsb$region01=="Northeast"],boot.eclsb$sampwt[boot.eclsb$region01=="Northeast"])-weighted.mean(boot.eclsb$read.astar.mofa[boot.eclsb$region01=="Northeast"],boot.eclsb$sampwt[boot.eclsb$region01=="Northeast"])

		boot.ate.mt.ne<-weighted.mean(boot.eclsb$math.astar[boot.eclsb$region01=="Northeast"],boot.eclsb$sampwt[boot.eclsb$region01=="Northeast"])-weighted.mean(boot.eclsb$math.a[boot.eclsb$region01=="Northeast"],boot.eclsb$sampwt[boot.eclsb$region01=="Northeast"])
		boot.nde.mt.ne<-weighted.mean(boot.eclsb$math.astar.mofa[boot.eclsb$region01=="Northeast"],boot.eclsb$sampwt[boot.eclsb$region01=="Northeast"])-weighted.mean(boot.eclsb$math.a[boot.eclsb$region01=="Northeast"],boot.eclsb$sampwt[boot.eclsb$region01=="Northeast"])
		boot.nie.mt.ne<-weighted.mean(boot.eclsb$math.astar[boot.eclsb$region01=="Northeast"],boot.eclsb$sampwt[boot.eclsb$region01=="Northeast"])-weighted.mean(boot.eclsb$math.astar.mofa[boot.eclsb$region01=="Northeast"],boot.eclsb$sampwt[boot.eclsb$region01=="Northeast"])

		boot.ate.rd.mw<-weighted.mean(boot.eclsb$read.astar[boot.eclsb$region01=="Midwest"],boot.eclsb$sampwt[boot.eclsb$region01=="Midwest"])-weighted.mean(boot.eclsb$read.a[boot.eclsb$region01=="Midwest"],boot.eclsb$sampwt[boot.eclsb$region01=="Midwest"])
		boot.nde.rd.mw<-weighted.mean(boot.eclsb$read.astar.mofa[boot.eclsb$region01=="Midwest"],boot.eclsb$sampwt[boot.eclsb$region01=="Midwest"])-weighted.mean(boot.eclsb$read.a[boot.eclsb$region01=="Midwest"],boot.eclsb$sampwt[boot.eclsb$region01=="Midwest"])
		boot.nie.rd.mw<-weighted.mean(boot.eclsb$read.astar[boot.eclsb$region01=="Midwest"],boot.eclsb$sampwt[boot.eclsb$region01=="Midwest"])-weighted.mean(boot.eclsb$read.astar.mofa[boot.eclsb$region01=="Midwest"],boot.eclsb$sampwt[boot.eclsb$region01=="Midwest"])

		boot.ate.mt.mw<-weighted.mean(boot.eclsb$math.astar[boot.eclsb$region01=="Midwest"],boot.eclsb$sampwt[boot.eclsb$region01=="Midwest"])-weighted.mean(boot.eclsb$math.a[boot.eclsb$region01=="Midwest"],boot.eclsb$sampwt[boot.eclsb$region01=="Midwest"])
		boot.nde.mt.mw<-weighted.mean(boot.eclsb$math.astar.mofa[boot.eclsb$region01=="Midwest"],boot.eclsb$sampwt[boot.eclsb$region01=="Midwest"])-weighted.mean(boot.eclsb$math.a[boot.eclsb$region01=="Midwest"],boot.eclsb$sampwt[boot.eclsb$region01=="Midwest"])
		boot.nie.mt.mw<-weighted.mean(boot.eclsb$math.astar[boot.eclsb$region01=="Midwest"],boot.eclsb$sampwt[boot.eclsb$region01=="Midwest"])-weighted.mean(boot.eclsb$math.astar.mofa[boot.eclsb$region01=="Midwest"],boot.eclsb$sampwt[boot.eclsb$region01=="Midwest"])

		boot.ate.rd.so<-weighted.mean(boot.eclsb$read.astar[boot.eclsb$region01=="South"],boot.eclsb$sampwt[boot.eclsb$region01=="South"])-weighted.mean(boot.eclsb$read.a[boot.eclsb$region01=="South"],boot.eclsb$sampwt[boot.eclsb$region01=="South"])
		boot.nde.rd.so<-weighted.mean(boot.eclsb$read.astar.mofa[boot.eclsb$region01=="South"],boot.eclsb$sampwt[boot.eclsb$region01=="South"])-weighted.mean(boot.eclsb$read.a[boot.eclsb$region01=="South"],boot.eclsb$sampwt[boot.eclsb$region01=="South"])
		boot.nie.rd.so<-weighted.mean(boot.eclsb$read.astar[boot.eclsb$region01=="South"],boot.eclsb$sampwt[boot.eclsb$region01=="South"])-weighted.mean(boot.eclsb$read.astar.mofa[boot.eclsb$region01=="South"],boot.eclsb$sampwt[boot.eclsb$region01=="South"])

		boot.ate.mt.so<-weighted.mean(boot.eclsb$math.astar[boot.eclsb$region01=="South"],boot.eclsb$sampwt[boot.eclsb$region01=="South"])-weighted.mean(boot.eclsb$math.a[boot.eclsb$region01=="South"],boot.eclsb$sampwt[boot.eclsb$region01=="South"])
		boot.nde.mt.so<-weighted.mean(boot.eclsb$math.astar.mofa[boot.eclsb$region01=="South"],boot.eclsb$sampwt[boot.eclsb$region01=="South"])-weighted.mean(boot.eclsb$math.a[boot.eclsb$region01=="South"],boot.eclsb$sampwt[boot.eclsb$region01=="South"])
		boot.nie.mt.so<-weighted.mean(boot.eclsb$math.astar[boot.eclsb$region01=="South"],boot.eclsb$sampwt[boot.eclsb$region01=="South"])-weighted.mean(boot.eclsb$math.astar.mofa[boot.eclsb$region01=="South"],boot.eclsb$sampwt[boot.eclsb$region01=="South"])

		boot.ate.rd.we<-weighted.mean(boot.eclsb$read.astar[boot.eclsb$region01=="West"],boot.eclsb$sampwt[boot.eclsb$region01=="West"])-weighted.mean(boot.eclsb$read.a[boot.eclsb$region01=="West"],boot.eclsb$sampwt[boot.eclsb$region01=="West"])
		boot.nde.rd.we<-weighted.mean(boot.eclsb$read.astar.mofa[boot.eclsb$region01=="West"],boot.eclsb$sampwt[boot.eclsb$region01=="West"])-weighted.mean(boot.eclsb$read.a[boot.eclsb$region01=="West"],boot.eclsb$sampwt[boot.eclsb$region01=="West"])
		boot.nie.rd.we<-weighted.mean(boot.eclsb$read.astar[boot.eclsb$region01=="West"],boot.eclsb$sampwt[boot.eclsb$region01=="West"])-weighted.mean(boot.eclsb$read.astar.mofa[boot.eclsb$region01=="West"],boot.eclsb$sampwt[boot.eclsb$region01=="West"])

		boot.ate.mt.we<-weighted.mean(boot.eclsb$math.astar[boot.eclsb$region01=="West"],boot.eclsb$sampwt[boot.eclsb$region01=="West"])-weighted.mean(boot.eclsb$math.a[boot.eclsb$region01=="West"],boot.eclsb$sampwt[boot.eclsb$region01=="West"])
		boot.nde.mt.we<-weighted.mean(boot.eclsb$math.astar.mofa[boot.eclsb$region01=="West"],boot.eclsb$sampwt[boot.eclsb$region01=="West"])-weighted.mean(boot.eclsb$math.a[boot.eclsb$region01=="West"],boot.eclsb$sampwt[boot.eclsb$region01=="West"])
		boot.nie.mt.we<-weighted.mean(boot.eclsb$math.astar[boot.eclsb$region01=="West"],boot.eclsb$sampwt[boot.eclsb$region01=="West"])-weighted.mean(boot.eclsb$math.astar.mofa[boot.eclsb$region01=="West"],boot.eclsb$sampwt[boot.eclsb$region01=="West"])
			
		return(list(
			boot.ate.rd.ne,boot.nde.rd.ne,boot.nie.rd.ne,
			boot.ate.mt.ne,boot.nde.mt.ne,boot.nie.mt.ne,
			boot.ate.rd.mw,boot.nde.rd.mw,boot.nie.rd.mw,
			boot.ate.mt.mw,boot.nde.mt.mw,boot.nie.mt.mw,
			boot.ate.rd.so,boot.nde.rd.so,boot.nie.rd.so,
			boot.ate.mt.so,boot.nde.mt.so,boot.nie.mt.so,
			boot.ate.rd.we,boot.nde.rd.we,boot.nie.rd.we,
			boot.ate.mt.we,boot.nde.mt.we,boot.nie.mt.we))
		}

	boot.est.w1<-matrix(unlist(boot.est.w1),ncol=24,byrow=TRUE)

	boot.ci.rd.ne<-rbind(boot.ci.rd.ne,boot.est.w1[,1:3])
	boot.ci.mt.ne<-rbind(boot.ci.mt.ne,boot.est.w1[,4:6])

	boot.ci.rd.mw<-rbind(boot.ci.rd.mw,boot.est.w1[,7:9])
	boot.ci.mt.mw<-rbind(boot.ci.mt.mw,boot.est.w1[,10:12])

	boot.ci.rd.so<-rbind(boot.ci.rd.so,boot.est.w1[,13:15])
	boot.ci.mt.so<-rbind(boot.ci.mt.so,boot.est.w1[,16:18])

	boot.ci.rd.we<-rbind(boot.ci.rd.so,boot.est.w1[,19:21])
	boot.ci.mt.we<-rbind(boot.ci.mt.so,boot.est.w1[,22:24])

	stopCluster(my.cluster)
	rm(my.cluster)
	}

## COMBINE MI ESTIMATES ###
est.rd.ne<-est.mt.ne<-matrix(data=NA,nrow=3,ncol=3)
est.rd.mw<-est.mt.mw<-matrix(data=NA,nrow=3,ncol=3)
est.rd.so<-est.mt.so<-matrix(data=NA,nrow=3,ncol=3)
est.rd.we<-est.mt.we<-matrix(data=NA,nrow=3,ncol=3)

for (i in 1:3) { 
	
	est.rd.ne[i,1]<-round(mean(miest.rd.ne[,i]),digits=3)
	est.rd.ne[i,2]<-round(quantile(boot.ci.rd.ne[,i],prob=0.025),digits=3)
	est.rd.ne[i,3]<-round(quantile(boot.ci.rd.ne[,i],prob=0.975),digits=3)

	est.mt.ne[i,1]<-round(mean(miest.mt.ne[,i]),digits=3)
	est.mt.ne[i,2]<-round(quantile(boot.ci.mt.ne[,i],prob=0.025),digits=3)
	est.mt.ne[i,3]<-round(quantile(boot.ci.mt.ne[,i],prob=0.975),digits=3)

	est.rd.mw[i,1]<-round(mean(miest.rd.mw[,i]),digits=3)
	est.rd.mw[i,2]<-round(quantile(boot.ci.rd.mw[,i],prob=0.025),digits=3)
	est.rd.mw[i,3]<-round(quantile(boot.ci.rd.mw[,i],prob=0.975),digits=3)

	est.mt.mw[i,1]<-round(mean(miest.mt.mw[,i]),digits=3)
	est.mt.mw[i,2]<-round(quantile(boot.ci.mt.mw[,i],prob=0.025),digits=3)
	est.mt.mw[i,3]<-round(quantile(boot.ci.mt.mw[,i],prob=0.975),digits=3)

	est.rd.so[i,1]<-round(mean(miest.rd.so[,i]),digits=3)
	est.rd.so[i,2]<-round(quantile(boot.ci.rd.so[,i],prob=0.025),digits=3)
	est.rd.so[i,3]<-round(quantile(boot.ci.rd.so[,i],prob=0.975),digits=3)

	est.mt.so[i,1]<-round(mean(miest.mt.so[,i]),digits=3)
	est.mt.so[i,2]<-round(quantile(boot.ci.mt.so[,i],prob=0.025),digits=3)
	est.mt.so[i,3]<-round(quantile(boot.ci.mt.so[,i],prob=0.975),digits=3)

	est.rd.we[i,1]<-round(mean(miest.rd.we[,i]),digits=3)
	est.rd.we[i,2]<-round(quantile(boot.ci.rd.we[,i],prob=0.025),digits=3)
	est.rd.we[i,3]<-round(quantile(boot.ci.rd.we[,i],prob=0.975),digits=3)

	est.mt.we[i,1]<-round(mean(miest.mt.we[,i]),digits=3)
	est.mt.we[i,2]<-round(quantile(boot.ci.mt.we[,i],prob=0.025),digits=3)
	est.mt.we[i,3]<-round(quantile(boot.ci.mt.we[,i],prob=0.975),digits=3)
	}

rlabel<-c('ATE','NDE','NIE')
output.rd.ne<-data.frame(est.rd.ne,row.names=rlabel)
output.mt.ne<-data.frame(est.mt.ne,row.names=rlabel)
output.rd.mw<-data.frame(est.rd.mw,row.names=rlabel)
output.mt.mw<-data.frame(est.mt.mw,row.names=rlabel)
output.rd.so<-data.frame(est.rd.so,row.names=rlabel)
output.mt.so<-data.frame(est.mt.so,row.names=rlabel)
output.rd.we<-data.frame(est.rd.we,row.names=rlabel)
output.mt.we<-data.frame(est.mt.we,row.names=rlabel)

colnames(output.rd.ne)<-colnames(output.mt.ne)<-c('estimate','ll.pct.95ci','ul.pct.95ci')
colnames(output.rd.mw)<-colnames(output.mt.mw)<-c('estimate','ll.pct.95ci','ul.pct.95ci')
colnames(output.rd.so)<-colnames(output.mt.so)<-c('estimate','ll.pct.95ci','ul.pct.95ci')
colnames(output.rd.we)<-colnames(output.mt.we)<-c('estimate','ll.pct.95ci','ul.pct.95ci')

### PRINT RESULTS ###
sink("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\programs\\_LOGS\\27_create_table_S7_log.txt")

table(eclsb.mi[which(eclsb.mi$minum==1),"region01"])

cat("===========================================\n")
cat("Northeast\n")
cat("===========================================\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Reading Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.rd.ne)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Math Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.mt.ne)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("===========================================\n")
cat(" \n")
cat(" \n")
cat(" \n")
cat("===========================================\n")
cat("Midwest\n")
cat("===========================================\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Reading Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.rd.mw)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Math Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.mt.mw)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("===========================================\n")
cat(" \n")
cat(" \n")
cat(" \n")
cat("===========================================\n")
cat("South\n")
cat("===========================================\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Reading Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.rd.so)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Math Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.mt.so)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("===========================================\n")
cat(" \n")
cat(" \n")
cat(" \n")
cat("===========================================\n")
cat("West\n")
cat("===========================================\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Reading Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.rd.we)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Math Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.mt.we)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("===========================================\n")

print(startTime)
print(Sys.time())

sink()

