################################################
################################################
##                                            ##
## PROGRAM NAME: 20_create_figure_S2          ##
## AUTHOR: BP, GW                             ##
## DESCRIPTION:                               ##
##                                            ##
##  create plots of covariates, no mediators  ##
##  SHAP importance for test scores           ##
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
  "Age.at.assessment")

##### STANDARDIZE VARS #####
eclsb.mi$readtheta05<-(eclsb.mi$readtheta05-weighted.mean(eclsb.mi$readtheta05,eclsb.mi$sampwt))/sqrt(wtd.var(eclsb.mi$readtheta05,eclsb.mi$sampwt))
eclsb.mi$maththeta05<-(eclsb.mi$maththeta05-weighted.mean(eclsb.mi$maththeta05,eclsb.mi$sampwt))/sqrt(wtd.var(eclsb.mi$maththeta05,eclsb.mi$sampwt))

##### COMPUTE SHAP IMPORTANCE #####
miest.rd.w1<-miest.mt.w1<-matrix(data=NA,nrow=length(c(vars.covw1, vars.base)),ncol=nmi)

pfun <- function(object, newdata) { predict(object,data=newdata)$pred }

for (i in 1:nmi) {

	### LOAD MI DATA ###
	print(c("nmi=",i))
	eclsb<-eclsb.mi[which(eclsb.mi$minum==i),]

	### TUNE HYPERPARAMETERS ###
	registerDoSEQ()

	cntrl<-trainControl(method="CV",number=5)

	rf.m1.grid<-expand.grid(
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

	### FIT MODELS ###
	read.train<-eclsb[,"readtheta05"]
	math.train<-eclsb[,"maththeta05"]
	x.train<-eclsb[,c(vars.base,vars.covw1)]
	wts.train<-eclsb[,"sampwt"]

	m1.rd<-ranger(
		y=read.train,
		x=x.train,
		num.trees=ntrees,
		min.node.size=rf.rd.m1$bestTune[,"min.node.size"],
		mtry=rf.rd.m1$bestTune[,"mtry"],
		splitrule="variance",
		respect.unordered.factors=TRUE,
		case.weights=wts.train,
		seed=8675309)

	m1.mt<-ranger(
		y=math.train,
		x=x.train,
		num.trees=ntrees,
		min.node.size=rf.mt.m1$bestTune[,"min.node.size"],
		mtry=rf.mt.m1$bestTune[,"mtry"],
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
	shap.cov.rd<-shap.rd[c(vars.base, vars.covw1)]
	shap.cov.mt<-shap.mt[c(vars.base, vars.covw1)]
	shap.cov.imp.rd<-data.frame(label=var.label, importance=apply(shap.cov.rd,MARGIN=2,FUN=function(x) mean(abs(x))))
	shap.cov.imp.mt<-data.frame(label=var.label, importance=apply(shap.cov.mt,MARGIN=2,FUN=function(x) mean(abs(x))))
	miest.rd.w1[,i]<-shap.cov.imp.rd[,"importance"]
	miest.mt.w1[,i]<-shap.cov.imp.mt[,"importance"]
	}

### COMBINE MI ESTIMATES ###
est.rd.w1<-est.mt.w1<-matrix(data=NA,nrow=length(c(vars.covw1, vars.base)),ncol=2)

for (i in 1:length(c(vars.covw1, vars.base))) { 

	est.rd.w1[i,1]<-est.mt.w1[i,1]<-var.label[i]
	est.rd.w1[i,2]<-round(mean(miest.rd.w1[i,]),digits=4)
	est.mt.w1[i,2]<-round(mean(miest.mt.w1[i,]),digits=4)
	}

### PRINT RESULTS ###
sink("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\programs\\_LOGS\\20_create_figure_S2_log.txt")

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
#load("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\_TEMP\\shap_est_covars.RData")

theme_set(theme_bw())

plot1<-ggplot(output.rd.w1, aes(reorder(label, importance), importance)) +
		geom_col() +
		coord_flip() +
		xlab(" ") +
		ylab("mean(|SHAP value|)") +
		ylim(0.0, 0.25)+
		ggtitle("A. Reading Scores")

plot2<-ggplot(output.mt.w1, aes(reorder(label, importance), importance)) +
		geom_col() +
		coord_flip() +
		xlab(" ") +
		ylab("mean(|SHAP value|)") +
		ylim(0.0, 0.25)+
		ggtitle("B. Math Scores")

tiff("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\figures\\figure_S2.tiff",
	height=9,
	width=9,
	units='in',
	res=600)

grid.arrange(plot1,plot2,ncol=2)

dev.off()

print(startTime)
print(Sys.time())

sink()

save(output.rd.w1,output.mt.w1,file="C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\_TEMP\\shap_est_covars.RData")

registerDoSEQ()
