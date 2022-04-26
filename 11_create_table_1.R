################################################
################################################
##                                            ##
## PROGRAM NAME: 11_create_table_1            ##
## AUTHOR: GW                                 ##
## DESCRIPTION:                               ##
##                                            ##
##  create table of nh effect estimates       ##
##  using random forests and reg imputation   ## 
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
miest.rd.w1<-miest.mt.w1<-matrix(data=NA,nrow=nmi,ncol=3)
boot.ci.rd<-boot.ci.mt<-NULL

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
	read.astar<-predict(read.giv.c.a,c.a.imp)$pred
	math.astar<-predict(math.giv.c.a,c.a.imp)$pred
	c.a.imp$nhpovrt01<-a	
	read.a<-predict(read.giv.c.a,c.a.imp)$pred
	math.a<-predict(math.giv.c.a,c.a.imp)$pred

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
	read.astar.mofa<-predict(read.uhat.astar.m,c.a.imp)$pred
	math.astar.mofa<-predict(math.uhat.astar.m,c.a.imp)$pred

	# STEP 3: COMPUTE ATE, NDE, NIE #
	miest.rd.w1[i,1]<-weighted.mean(read.astar,wts.train)-weighted.mean(read.a,wts.train)
	miest.rd.w1[i,2]<-weighted.mean(read.astar.mofa,wts.train)-weighted.mean(read.a,wts.train)
	miest.rd.w1[i,3]<-weighted.mean(read.astar,wts.train)-weighted.mean(read.astar.mofa,wts.train)

	miest.mt.w1[i,1]<-weighted.mean(math.astar,wts.train)-weighted.mean(math.a,wts.train)
	miest.mt.w1[i,2]<-weighted.mean(math.astar.mofa,wts.train)-weighted.mean(math.a,wts.train)
	miest.mt.w1[i,3]<-weighted.mean(math.astar,wts.train)-weighted.mean(math.astar.mofa,wts.train)

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
		boot.read.astar<-predict(boot.read.giv.c.a,boot.c.a.imp)$pred
		boot.math.astar<-predict(boot.math.giv.c.a,boot.c.a.imp)$pred
		boot.c.a.imp$nhpovrt01<-a	
		boot.read.a<-predict(boot.read.giv.c.a,boot.c.a.imp)$pred
		boot.math.a<-predict(boot.math.giv.c.a,boot.c.a.imp)$pred

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
		boot.read.astar.mofa<-predict(boot.read.uhat.astar.m,boot.c.a.imp)$pred
		boot.math.astar.mofa<-predict(boot.math.uhat.astar.m,boot.c.a.imp)$pred

		boot.ate.rd.w1<-weighted.mean(boot.read.astar,boot.wts.train)-weighted.mean(boot.read.a,boot.wts.train)
		boot.nde.rd.w1<-weighted.mean(boot.read.astar.mofa,boot.wts.train)-weighted.mean(boot.read.a,boot.wts.train)
		boot.nie.rd.w1<-weighted.mean(boot.read.astar,boot.wts.train)-weighted.mean(boot.read.astar.mofa,boot.wts.train)

		boot.ate.mt.w1<-weighted.mean(boot.math.astar,boot.wts.train)-weighted.mean(boot.math.a,boot.wts.train)
		boot.nde.mt.w1<-weighted.mean(boot.math.astar.mofa,boot.wts.train)-weighted.mean(boot.math.a,boot.wts.train)
		boot.nie.mt.w1<-weighted.mean(boot.math.astar,boot.wts.train)-weighted.mean(boot.math.astar.mofa,boot.wts.train)
			
		return(list(
			boot.ate.rd.w1,
			boot.nde.rd.w1,
			boot.nie.rd.w1,
			boot.ate.mt.w1,
			boot.nde.mt.w1,
			boot.nie.mt.w1))
		}

	boot.est.w1<-matrix(unlist(boot.est.w1),ncol=6,byrow=TRUE)

	boot.ci.rd<-rbind(boot.ci.rd,boot.est.w1[,1:3])
	boot.ci.mt<-rbind(boot.ci.mt,boot.est.w1[,4:6])

	stopCluster(my.cluster)
	rm(my.cluster)
	}

## COMBINE MI ESTIMATES ###
est.rd.w1<-est.mt.w1<-matrix(data=NA,nrow=3,ncol=3)

for (i in 1:3) { 
	
	est.rd.w1[i,1]<-round(mean(miest.rd.w1[,i]),digits=3)
	est.rd.w1[i,2]<-round(quantile(boot.ci.rd[,i],prob=0.025),digits=3)
	est.rd.w1[i,3]<-round(quantile(boot.ci.rd[,i],prob=0.975),digits=3)

	est.mt.w1[i,1]<-round(mean(miest.mt.w1[,i]),digits=3)
	est.mt.w1[i,2]<-round(quantile(boot.ci.mt[,i],prob=0.025),digits=3)
	est.mt.w1[i,3]<-round(quantile(boot.ci.mt[,i],prob=0.975),digits=3)
	}

rlabel<-c('ATE','NDE','NIE')
output.rd.w1<-data.frame(est.rd.w1,row.names=rlabel)
output.mt.w1<-data.frame(est.mt.w1,row.names=rlabel)
colnames(output.rd.w1)<-colnames(output.mt.w1)<-c('estimate','ll.pct.95ci','ul.pct.95ci')

### PRINT RESULTS ###
sink("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\programs\\_LOGS\\11_create_table_1_log.txt")

cat("===========================================\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Reading Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.rd.w1)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Math Test Scores\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
print(output.mt.w1)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("===========================================\n")

print(startTime)
print(Sys.time())

sink()

save(output.rd.w1,output.mt.w1,file="C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\_TEMP\\efx_est.RData")

