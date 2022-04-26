################################################
################################################
##                                            ##
## PROGRAM NAME: 25_create_table_S5           ##
## AUTHOR: GW                                 ##
## DESCRIPTION:                               ##
##                                            ##
##  create table of nh effect estimates by    ##
##  race using random forests and reg         ##
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
miest.rd.wht<-miest.mt.wht<-matrix(data=NA,nrow=nmi,ncol=3)
miest.rd.blk<-miest.mt.blk<-matrix(data=NA,nrow=nmi,ncol=3)
miest.rd.hsp<-miest.mt.hsp<-matrix(data=NA,nrow=nmi,ncol=3)
miest.rd.oth<-miest.mt.oth<-matrix(data=NA,nrow=nmi,ncol=3)

boot.ci.rd.wht<-boot.ci.mt.wht<-NULL
boot.ci.rd.blk<-boot.ci.mt.blk<-NULL
boot.ci.rd.hsp<-boot.ci.mt.hsp<-NULL
boot.ci.rd.oth<-boot.ci.mt.oth<-NULL

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

	miest.rd.wht[i,1]<-weighted.mean(eclsb$read.astar[eclsb$race=="White, non-Hispanic"],eclsb$sampwt[eclsb$race=="White, non-Hispanic"])-weighted.mean(eclsb$read.a[eclsb$race=="White, non-Hispanic"],eclsb$sampwt[eclsb$race=="White, non-Hispanic"])
	miest.rd.wht[i,2]<-weighted.mean(eclsb$read.astar.mofa[eclsb$race=="White, non-Hispanic"],eclsb$sampwt[eclsb$race=="White, non-Hispanic"])-weighted.mean(eclsb$read.a[eclsb$race=="White, non-Hispanic"],eclsb$sampwt[eclsb$race=="White, non-Hispanic"])
	miest.rd.wht[i,3]<-weighted.mean(eclsb$read.astar[eclsb$race=="White, non-Hispanic"],eclsb$sampwt[eclsb$race=="White, non-Hispanic"])-weighted.mean(eclsb$read.astar.mofa[eclsb$race=="White, non-Hispanic"],eclsb$sampwt[eclsb$race=="White, non-Hispanic"])

	miest.mt.wht[i,1]<-weighted.mean(eclsb$math.astar[eclsb$race=="White, non-Hispanic"],eclsb$sampwt[eclsb$race=="White, non-Hispanic"])-weighted.mean(eclsb$math.a[eclsb$race=="White, non-Hispanic"],eclsb$sampwt[eclsb$race=="White, non-Hispanic"])
	miest.mt.wht[i,2]<-weighted.mean(eclsb$math.astar.mofa[eclsb$race=="White, non-Hispanic"],eclsb$sampwt[eclsb$race=="White, non-Hispanic"])-weighted.mean(eclsb$math.a[eclsb$race=="White, non-Hispanic"],eclsb$sampwt[eclsb$race=="White, non-Hispanic"])
	miest.mt.wht[i,3]<-weighted.mean(eclsb$math.astar[eclsb$race=="White, non-Hispanic"],eclsb$sampwt[eclsb$race=="White, non-Hispanic"])-weighted.mean(eclsb$math.astar.mofa[eclsb$race=="White, non-Hispanic"],eclsb$sampwt[eclsb$race=="White, non-Hispanic"])

	miest.rd.blk[i,1]<-weighted.mean(eclsb$read.astar[eclsb$race=="Black, non-Hispanic"],eclsb$sampwt[eclsb$race=="Black, non-Hispanic"])-weighted.mean(eclsb$read.a[eclsb$race=="Black, non-Hispanic"],eclsb$sampwt[eclsb$race=="Black, non-Hispanic"])
	miest.rd.blk[i,2]<-weighted.mean(eclsb$read.astar.mofa[eclsb$race=="Black, non-Hispanic"],eclsb$sampwt[eclsb$race=="Black, non-Hispanic"])-weighted.mean(eclsb$read.a[eclsb$race=="Black, non-Hispanic"],eclsb$sampwt[eclsb$race=="Black, non-Hispanic"])
	miest.rd.blk[i,3]<-weighted.mean(eclsb$read.astar[eclsb$race=="Black, non-Hispanic"],eclsb$sampwt[eclsb$race=="Black, non-Hispanic"])-weighted.mean(eclsb$read.astar.mofa[eclsb$race=="Black, non-Hispanic"],eclsb$sampwt[eclsb$race=="Black, non-Hispanic"])

	miest.mt.blk[i,1]<-weighted.mean(eclsb$math.astar[eclsb$race=="Black, non-Hispanic"],eclsb$sampwt[eclsb$race=="Black, non-Hispanic"])-weighted.mean(eclsb$math.a[eclsb$race=="Black, non-Hispanic"],eclsb$sampwt[eclsb$race=="Black, non-Hispanic"])
	miest.mt.blk[i,2]<-weighted.mean(eclsb$math.astar.mofa[eclsb$race=="Black, non-Hispanic"],eclsb$sampwt[eclsb$race=="Black, non-Hispanic"])-weighted.mean(eclsb$math.a[eclsb$race=="Black, non-Hispanic"],eclsb$sampwt[eclsb$race=="Black, non-Hispanic"])
	miest.mt.blk[i,3]<-weighted.mean(eclsb$math.astar[eclsb$race=="Black, non-Hispanic"],eclsb$sampwt[eclsb$race=="Black, non-Hispanic"])-weighted.mean(eclsb$math.astar.mofa[eclsb$race=="Black, non-Hispanic"],eclsb$sampwt[eclsb$race=="Black, non-Hispanic"])

	miest.rd.hsp[i,1]<-weighted.mean(eclsb$read.astar[eclsb$race=="Hispanic"],eclsb$sampwt[eclsb$race=="Hispanic"])-weighted.mean(eclsb$read.a[eclsb$race=="Hispanic"],eclsb$sampwt[eclsb$race=="Hispanic"])
	miest.rd.hsp[i,2]<-weighted.mean(eclsb$read.astar.mofa[eclsb$race=="Hispanic"],eclsb$sampwt[eclsb$race=="Hispanic"])-weighted.mean(eclsb$read.a[eclsb$race=="Hispanic"],eclsb$sampwt[eclsb$race=="Hispanic"])
	miest.rd.hsp[i,3]<-weighted.mean(eclsb$read.astar[eclsb$race=="Hispanic"],eclsb$sampwt[eclsb$race=="Hispanic"])-weighted.mean(eclsb$read.astar.mofa[eclsb$race=="Hispanic"],eclsb$sampwt[eclsb$race=="Hispanic"])

	miest.mt.hsp[i,1]<-weighted.mean(eclsb$math.astar[eclsb$race=="Hispanic"],eclsb$sampwt[eclsb$race=="Hispanic"])-weighted.mean(eclsb$math.a[eclsb$race=="Hispanic"],eclsb$sampwt[eclsb$race=="Hispanic"])
	miest.mt.hsp[i,2]<-weighted.mean(eclsb$math.astar.mofa[eclsb$race=="Hispanic"],eclsb$sampwt[eclsb$race=="Hispanic"])-weighted.mean(eclsb$math.a[eclsb$race=="Hispanic"],eclsb$sampwt[eclsb$race=="Hispanic"])
	miest.mt.hsp[i,3]<-weighted.mean(eclsb$math.astar[eclsb$race=="Hispanic"],eclsb$sampwt[eclsb$race=="Hispanic"])-weighted.mean(eclsb$math.astar.mofa[eclsb$race=="Hispanic"],eclsb$sampwt[eclsb$race=="Hispanic"])

	miest.rd.oth[i,1]<-weighted.mean(eclsb$read.astar[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"],eclsb$sampwt[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"])-weighted.mean(eclsb$read.a[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"],eclsb$sampwt[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"])
	miest.rd.oth[i,2]<-weighted.mean(eclsb$read.astar.mofa[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"],eclsb$sampwt[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"])-weighted.mean(eclsb$read.a[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"],eclsb$sampwt[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"])
	miest.rd.oth[i,3]<-weighted.mean(eclsb$read.astar[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"],eclsb$sampwt[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"])-weighted.mean(eclsb$read.astar.mofa[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"],eclsb$sampwt[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"])

	miest.mt.oth[i,1]<-weighted.mean(eclsb$math.astar[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"],eclsb$sampwt[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"])-weighted.mean(eclsb$math.a[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"],eclsb$sampwt[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"])
	miest.mt.oth[i,2]<-weighted.mean(eclsb$math.astar.mofa[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"],eclsb$sampwt[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"])-weighted.mean(eclsb$math.a[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"],eclsb$sampwt[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"])
	miest.mt.oth[i,3]<-weighted.mean(eclsb$math.astar[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"],eclsb$sampwt[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"])-weighted.mean(eclsb$math.astar.mofa[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"],eclsb$sampwt[eclsb$race=="Asian, non-Hispanic" | eclsb$race=="Other"])

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

		boot.ate.rd.wht<-weighted.mean(boot.eclsb$read.astar[boot.eclsb$race=="White, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="White, non-Hispanic"])-weighted.mean(boot.eclsb$read.a[boot.eclsb$race=="White, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="White, non-Hispanic"])
		boot.nde.rd.wht<-weighted.mean(boot.eclsb$read.astar.mofa[boot.eclsb$race=="White, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="White, non-Hispanic"])-weighted.mean(boot.eclsb$read.a[boot.eclsb$race=="White, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="White, non-Hispanic"])
		boot.nie.rd.wht<-weighted.mean(boot.eclsb$read.astar[boot.eclsb$race=="White, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="White, non-Hispanic"])-weighted.mean(boot.eclsb$read.astar.mofa[boot.eclsb$race=="White, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="White, non-Hispanic"])

		boot.ate.mt.wht<-weighted.mean(boot.eclsb$math.astar[boot.eclsb$race=="White, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="White, non-Hispanic"])-weighted.mean(boot.eclsb$math.a[boot.eclsb$race=="White, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="White, non-Hispanic"])
		boot.nde.mt.wht<-weighted.mean(boot.eclsb$math.astar.mofa[boot.eclsb$race=="White, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="White, non-Hispanic"])-weighted.mean(boot.eclsb$math.a[boot.eclsb$race=="White, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="White, non-Hispanic"])
		boot.nie.mt.wht<-weighted.mean(boot.eclsb$math.astar[boot.eclsb$race=="White, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="White, non-Hispanic"])-weighted.mean(boot.eclsb$math.astar.mofa[boot.eclsb$race=="White, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="White, non-Hispanic"])

		boot.ate.rd.blk<-weighted.mean(boot.eclsb$read.astar[boot.eclsb$race=="Black, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Black, non-Hispanic"])-weighted.mean(boot.eclsb$read.a[boot.eclsb$race=="Black, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Black, non-Hispanic"])
		boot.nde.rd.blk<-weighted.mean(boot.eclsb$read.astar.mofa[boot.eclsb$race=="Black, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Black, non-Hispanic"])-weighted.mean(boot.eclsb$read.a[boot.eclsb$race=="Black, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Black, non-Hispanic"])
		boot.nie.rd.blk<-weighted.mean(boot.eclsb$read.astar[boot.eclsb$race=="Black, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Black, non-Hispanic"])-weighted.mean(boot.eclsb$read.astar.mofa[boot.eclsb$race=="Black, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Black, non-Hispanic"])

		boot.ate.mt.blk<-weighted.mean(boot.eclsb$math.astar[boot.eclsb$race=="Black, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Black, non-Hispanic"])-weighted.mean(boot.eclsb$math.a[boot.eclsb$race=="Black, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Black, non-Hispanic"])
		boot.nde.mt.blk<-weighted.mean(boot.eclsb$math.astar.mofa[boot.eclsb$race=="Black, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Black, non-Hispanic"])-weighted.mean(boot.eclsb$math.a[boot.eclsb$race=="Black, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Black, non-Hispanic"])
		boot.nie.mt.blk<-weighted.mean(boot.eclsb$math.astar[boot.eclsb$race=="Black, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Black, non-Hispanic"])-weighted.mean(boot.eclsb$math.astar.mofa[boot.eclsb$race=="Black, non-Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Black, non-Hispanic"])

		boot.ate.rd.hsp<-weighted.mean(boot.eclsb$read.astar[boot.eclsb$race=="Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Hispanic"])-weighted.mean(boot.eclsb$read.a[boot.eclsb$race=="Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Hispanic"])
		boot.nde.rd.hsp<-weighted.mean(boot.eclsb$read.astar.mofa[boot.eclsb$race=="Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Hispanic"])-weighted.mean(boot.eclsb$read.a[boot.eclsb$race=="Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Hispanic"])
		boot.nie.rd.hsp<-weighted.mean(boot.eclsb$read.astar[boot.eclsb$race=="Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Hispanic"])-weighted.mean(boot.eclsb$read.astar.mofa[boot.eclsb$race=="Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Hispanic"])

		boot.ate.mt.hsp<-weighted.mean(boot.eclsb$math.astar[boot.eclsb$race=="Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Hispanic"])-weighted.mean(boot.eclsb$math.a[boot.eclsb$race=="Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Hispanic"])
		boot.nde.mt.hsp<-weighted.mean(boot.eclsb$math.astar.mofa[boot.eclsb$race=="Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Hispanic"])-weighted.mean(boot.eclsb$math.a[boot.eclsb$race=="Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Hispanic"])
		boot.nie.mt.hsp<-weighted.mean(boot.eclsb$math.astar[boot.eclsb$race=="Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Hispanic"])-weighted.mean(boot.eclsb$math.astar.mofa[boot.eclsb$race=="Hispanic"],boot.eclsb$sampwt[boot.eclsb$race=="Hispanic"])

		boot.ate.rd.oth<-weighted.mean(boot.eclsb$read.astar[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"],boot.eclsb$sampwt[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"])-weighted.mean(boot.eclsb$read.a[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"],boot.eclsb$sampwt[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"])
		boot.nde.rd.oth<-weighted.mean(boot.eclsb$read.astar.mofa[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"],boot.eclsb$sampwt[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"])-weighted.mean(boot.eclsb$read.a[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"],boot.eclsb$sampwt[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"])
		boot.nie.rd.oth<-weighted.mean(boot.eclsb$read.astar[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"],boot.eclsb$sampwt[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"])-weighted.mean(boot.eclsb$read.astar.mofa[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"],boot.eclsb$sampwt[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"])

		boot.ate.mt.oth<-weighted.mean(boot.eclsb$math.astar[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"],boot.eclsb$sampwt[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"])-weighted.mean(boot.eclsb$math.a[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"],boot.eclsb$sampwt[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"])
		boot.nde.mt.oth<-weighted.mean(boot.eclsb$math.astar.mofa[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"],boot.eclsb$sampwt[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"])-weighted.mean(boot.eclsb$math.a[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"],boot.eclsb$sampwt[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"])
		boot.nie.mt.oth<-weighted.mean(boot.eclsb$math.astar[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"],boot.eclsb$sampwt[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"])-weighted.mean(boot.eclsb$math.astar.mofa[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"],boot.eclsb$sampwt[boot.eclsb$race=="Asian, non-Hispanic" | boot.eclsb$race=="Other"])
			
		return(list(
			boot.ate.rd.wht,boot.nde.rd.wht,boot.nie.rd.wht,
			boot.ate.mt.wht,boot.nde.mt.wht,boot.nie.mt.wht,
			boot.ate.rd.blk,boot.nde.rd.blk,boot.nie.rd.blk,
			boot.ate.mt.blk,boot.nde.mt.blk,boot.nie.mt.blk,
			boot.ate.rd.hsp,boot.nde.rd.hsp,boot.nie.rd.hsp,
			boot.ate.mt.hsp,boot.nde.mt.hsp,boot.nie.mt.hsp,
			boot.ate.rd.oth,boot.nde.rd.oth,boot.nie.rd.oth,
			boot.ate.mt.oth,boot.nde.mt.oth,boot.nie.mt.oth))
		}

	boot.est.w1<-matrix(unlist(boot.est.w1),ncol=24,byrow=TRUE)

	boot.ci.rd.wht<-rbind(boot.ci.rd.wht,boot.est.w1[,1:3])
	boot.ci.mt.wht<-rbind(boot.ci.mt.wht,boot.est.w1[,4:6])

	boot.ci.rd.blk<-rbind(boot.ci.rd.blk,boot.est.w1[,7:9])
	boot.ci.mt.blk<-rbind(boot.ci.mt.blk,boot.est.w1[,10:12])

	boot.ci.rd.hsp<-rbind(boot.ci.rd.hsp,boot.est.w1[,13:15])
	boot.ci.mt.hsp<-rbind(boot.ci.mt.hsp,boot.est.w1[,16:18])

	boot.ci.rd.oth<-rbind(boot.ci.rd.hsp,boot.est.w1[,19:21])
	boot.ci.mt.oth<-rbind(boot.ci.mt.hsp,boot.est.w1[,22:24])

	stopCluster(my.cluster)
	rm(my.cluster)
	}

## COMBINE MI ESTIMATES ###
est.rd.wht<-est.mt.wht<-matrix(data=NA,nrow=3,ncol=3)
est.rd.blk<-est.mt.blk<-matrix(data=NA,nrow=3,ncol=3)
est.rd.hsp<-est.mt.hsp<-matrix(data=NA,nrow=3,ncol=3)
est.rd.oth<-est.mt.oth<-matrix(data=NA,nrow=3,ncol=3)

for (i in 1:3) { 
	
	est.rd.wht[i,1]<-round(mean(miest.rd.wht[,i]),digits=3)
	est.rd.wht[i,2]<-round(quantile(boot.ci.rd.wht[,i],prob=0.025),digits=3)
	est.rd.wht[i,3]<-round(quantile(boot.ci.rd.wht[,i],prob=0.975),digits=3)

	est.mt.wht[i,1]<-round(mean(miest.mt.wht[,i]),digits=3)
	est.mt.wht[i,2]<-round(quantile(boot.ci.mt.wht[,i],prob=0.025),digits=3)
	est.mt.wht[i,3]<-round(quantile(boot.ci.mt.wht[,i],prob=0.975),digits=3)

	est.rd.blk[i,1]<-round(mean(miest.rd.blk[,i]),digits=3)
	est.rd.blk[i,2]<-round(quantile(boot.ci.rd.blk[,i],prob=0.025),digits=3)
	est.rd.blk[i,3]<-round(quantile(boot.ci.rd.blk[,i],prob=0.975),digits=3)

	est.mt.blk[i,1]<-round(mean(miest.mt.blk[,i]),digits=3)
	est.mt.blk[i,2]<-round(quantile(boot.ci.mt.blk[,i],prob=0.025),digits=3)
	est.mt.blk[i,3]<-round(quantile(boot.ci.mt.blk[,i],prob=0.975),digits=3)

	est.rd.hsp[i,1]<-round(mean(miest.rd.hsp[,i]),digits=3)
	est.rd.hsp[i,2]<-round(quantile(boot.ci.rd.hsp[,i],prob=0.025),digits=3)
	est.rd.hsp[i,3]<-round(quantile(boot.ci.rd.hsp[,i],prob=0.975),digits=3)

	est.mt.hsp[i,1]<-round(mean(miest.mt.hsp[,i]),digits=3)
	est.mt.hsp[i,2]<-round(quantile(boot.ci.mt.hsp[,i],prob=0.025),digits=3)
	est.mt.hsp[i,3]<-round(quantile(boot.ci.mt.hsp[,i],prob=0.975),digits=3)

	est.rd.oth[i,1]<-round(mean(miest.rd.oth[,i]),digits=3)
	est.rd.oth[i,2]<-round(quantile(boot.ci.rd.oth[,i],prob=0.025),digits=3)
	est.rd.oth[i,3]<-round(quantile(boot.ci.rd.oth[,i],prob=0.975),digits=3)

	est.mt.oth[i,1]<-round(mean(miest.mt.oth[,i]),digits=3)
	est.mt.oth[i,2]<-round(quantile(boot.ci.mt.oth[,i],prob=0.025),digits=3)
	est.mt.oth[i,3]<-round(quantile(boot.ci.mt.oth[,i],prob=0.975),digits=3)
	}

rlabel<-c('ATE','NDE','NIE')
output.rd.wht<-data.frame(est.rd.wht,row.names=rlabel)
output.mt.wht<-data.frame(est.mt.wht,row.names=rlabel)
output.rd.blk<-data.frame(est.rd.blk,row.names=rlabel)
output.mt.blk<-data.frame(est.mt.blk,row.names=rlabel)
output.rd.hsp<-data.frame(est.rd.hsp,row.names=rlabel)
output.mt.hsp<-data.frame(est.mt.hsp,row.names=rlabel)
output.rd.oth<-data.frame(est.rd.oth,row.names=rlabel)
output.mt.oth<-data.frame(est.mt.oth,row.names=rlabel)

colnames(output.rd.wht)<-colnames(output.mt.wht)<-c('estimate','ll.pct.95ci','ul.pct.95ci')
colnames(output.rd.blk)<-colnames(output.mt.blk)<-c('estimate','ll.pct.95ci','ul.pct.95ci')
colnames(output.rd.hsp)<-colnames(output.mt.hsp)<-c('estimate','ll.pct.95ci','ul.pct.95ci')
colnames(output.rd.oth)<-colnames(output.mt.oth)<-c('estimate','ll.pct.95ci','ul.pct.95ci')

### PRINT RESULTS ###
sink("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\programs\\_LOGS\\25_create_table_S5_log.txt")

table(eclsb.mi[which(eclsb.mi$minum==1),"race"])

cat("===========================================\n")
cat("WHITES\n")
cat("===========================================\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Reading Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.rd.wht)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Math Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.mt.wht)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("===========================================\n")
cat(" \n")
cat(" \n")
cat(" \n")
cat("===========================================\n")
cat("BLACKS\n")
cat("===========================================\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Reading Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.rd.blk)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Math Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.mt.blk)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("===========================================\n")
cat(" \n")
cat(" \n")
cat(" \n")
cat("===========================================\n")
cat("HISPANICS\n")
cat("===========================================\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Reading Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.rd.hsp)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Math Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.mt.hsp)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("===========================================\n")
cat(" \n")
cat(" \n")
cat(" \n")
cat("===========================================\n")
cat("ASIANS + OTHER RACE\n")
cat("===========================================\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Reading Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.rd.oth)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Math Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.mt.oth)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("===========================================\n")

print(startTime)
print(Sys.time())

sink()

