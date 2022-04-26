################################################
################################################
##                                            ##
## PROGRAM NAME: 19_create_figure_S1          ##
## AUTHOR: GW, BP                             ##
## DESCRIPTION:                               ##
##                                            ##
##  create plots of all predictors' SHAP      ##
##  importance for test scores                ##
##                                            ##
################################################
################################################

rm(list=ls())

list.of.packages <- c(
	"foreach",
	"doParallel",
	"tidyverse",
	"doRNG",
	"sys",
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
	"gridExtra",
	"fastshap",
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
shap.sim<-50
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

var.label<-c(
	"Gender",
	"Race",
	"Twin",
	"Birth.weight",
	"Maternal.age",
	"Paternal.age",
	"WIC",
	"SNAP",
	"Medicaid",
	"TANF",
	"Nh.poverty",
	"Family.income",
	"Parental.education",
	"Parental.occupation",
	"Maternal.employment",
	"Paternal.employment",
	"Total.HH.members",
	"Biologicial.father",
	"Mother.married",
	"Owns.home",
	"Parental.involvement",
	"English.primary",
	"Population.density",
	"Urbanicity",
	"Region",
	"Age.at.baseline",
	"Age.at.assessment",
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
for (v in 1:length(vars.ntxw1)) {
	eclsb.mi[,vars.ntxw1[v]]<-(eclsb.mi[,vars.ntxw1[v]]-weighted.mean(eclsb.mi[,vars.ntxw1[v]],eclsb.mi$sampwt))/sqrt(wtd.var(eclsb.mi[,vars.ntxw1[v]],eclsb.mi$sampwt))
	}

eclsb.mi$readtheta05<-(eclsb.mi$readtheta05-weighted.mean(eclsb.mi$readtheta05,eclsb.mi$sampwt))/sqrt(wtd.var(eclsb.mi$readtheta05,eclsb.mi$sampwt))
eclsb.mi$maththeta05<-(eclsb.mi$maththeta05-weighted.mean(eclsb.mi$maththeta05,eclsb.mi$sampwt))/sqrt(wtd.var(eclsb.mi$maththeta05,eclsb.mi$sampwt))

##### COMPUTE SHAP IMPORTANCE #####
miest.rd.w1<-miest.mt.w1<-matrix(data=NA,nrow=length(c(vars.ntxw1, vars.covw1, vars.base)),ncol=nmi)

pfun <- function(object, newdata) { predict(object,data=newdata)$pred }

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

	### FIT MODELS ###
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
		x=x.train,
		num.trees=ntrees,
		min.node.size=rf.mt.m2$bestTune[,"min.node.size"],
		mtry=rf.mt.m2$bestTune[,"mtry"],
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=wts.train,
		seed=8675309)

	### SETUP PARALLEL PROCESSING ###
	cl<-makeCluster(n.cores,type="PSOCK")
	registerDoParallel(cl)
	clusterEvalQ(cl,library(ranger))
	clusterExport(cl,c("eclsb","shap.sim","pfun","m1.rd","m1.mt","x.train"),envir=environment())
	registerDoRNG(8675309)

	### COMPUTE SHAP VALUES ###
	shap.rd<-explain(m1.rd, X=x.train, pred_wrapper=pfun, nsim=shap.sim, .parallel=TRUE)
	shap.mt<-explain(m1.mt, X=x.train, pred_wrapper=pfun, nsim=shap.sim, .parallel=TRUE)

	stopCluster(cl)
	rm(cl)

	### COMPUTE SHAP IMPORTANCE ###
	shap.chem.rd<-shap.rd[c(vars.base, vars.covw1, vars.ntxw1)]
	shap.chem.mt<-shap.mt[c(vars.base, vars.covw1, vars.ntxw1)]
	shap.chem.imp.rd<-data.frame(label=var.label, importance=apply(shap.chem.rd,MARGIN=2,FUN=function(x) mean(abs(x))))
	shap.chem.imp.mt<-data.frame(label=var.label, importance=apply(shap.chem.mt,MARGIN=2,FUN=function(x) mean(abs(x))))
	miest.rd.w1[,i]<-shap.chem.imp.rd[,"importance"]
	miest.mt.w1[,i]<-shap.chem.imp.mt[,"importance"]
	}

### COMBINE MI ESTIMATES ###
est.rd.w1<-est.mt.w1<-matrix(data=NA,nrow=length(c(vars.ntxw1, vars.covw1, vars.base)),ncol=2)

for (i in 1:length(c(vars.ntxw1, vars.covw1, vars.base))) { 

	est.rd.w1[i,1]<-est.mt.w1[i,1]<-var.label[i]
	est.rd.w1[i,2]<-round(mean(miest.rd.w1[i,]),digits=4)
	est.mt.w1[i,2]<-round(mean(miest.mt.w1[i,]),digits=4)
	}

### PRINT RESULTS ###
sink("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\programs\\_LOGS\\19_create_figure_S1_log.txt")

output.rd.w1<-as.data.frame(est.rd.w1)
output.mt.w1<-as.data.frame(est.mt.w1)

colnames(output.rd.w1)<-colnames(output.mt.w1)<-c('label','importance')
rownames(output.rd.w1)<-rownames(output.mt.w1)<-var.label

output.rd.w1$importance<-as.numeric(output.rd.w1$importance)
output.mt.w1$importance<-as.numeric(output.mt.w1$importance)

output.rd.w1<-output.rd.w1[order(output.rd.w1$importance,decreasing=T),]
output.mt.w1<-output.mt.w1[order(output.mt.w1$importance,decreasing=T),]

print(output.rd.w1)
print(output.mt.w1)

##### PLOT RESULTS #####
#load("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\_TEMP\\shap_est_allvars.RData")

theme_set(theme_bw())

plot1<-ggplot(output.rd.w1[1:39,], aes(reorder(label, importance), importance)) +
		geom_col() +
		coord_flip() +
		xlab(" ") +
		ylab("mean(|SHAP value|)") +
		ylim(0.0, 0.20)+
		ggtitle("A. Reading Scores")

plot2<-ggplot(output.mt.w1[1:39,], aes(reorder(label, importance), importance)) +
		geom_col() +
		coord_flip() +
		xlab(" ") +
		ylab("mean(|SHAP value|)") +
		ylim(0.0, 0.20)+
		ggtitle("B. Math Scores")

tiff("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\figures\\figure_S1.1.tiff",
	height=9,
	width=9,
	units='in',
	res=600)

grid.arrange(plot1,plot2,ncol=2)

dev.off()

plot3<-ggplot(output.rd.w1[40:79,], aes(reorder(label, importance), importance)) +
  geom_col() +
  coord_flip() +
  xlab(" ") +
  ylab("mean(|SHAP value|)") +
  ylim(0.0, 0.20)+
  ggtitle("A. Reading Scores")

plot4<-ggplot(output.mt.w1[40:79,], aes(reorder(label, importance), importance)) +
  geom_col() +
  coord_flip() +
  xlab(" ") +
  ylab("mean(|SHAP value|)") +
  ylim(0.0, 0.20)+
  ggtitle("B. Math Scores")

tiff("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\figures\\figure_S1.2.tiff",
     height=9,
     width=9,
     units='in',
     res=600)

grid.arrange(plot3,plot4,ncol=2)

dev.off()

print(startTime)
print(Sys.time())

sink()

save(output.rd.w1,output.mt.w1,file="C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\_TEMP\\shap_est_allvars.RData")

