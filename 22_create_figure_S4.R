################################################
################################################
##                                            ##
## PROGRAM NAME: 22_create_figure_S4          ##
## AUTHOR: BP, GW                             ##
##                                            ##
## DESCRIPTION:                               ##
##                                            ##
##  create SHAP summary plots for baseline    ##
##  covariates only                           ##
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
	"SHAPforxgboost",
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

##### COMPUTE SHAP VALUES #####
pfun <- function(object, newdata) { predict(object,data=newdata)$pred }

### LOAD 1ST IMPUTATION ###
eclsb<-eclsb.mi[which(eclsb.mi$minum==1),]

### TUNE RF HYPERPARAMETERS ###
registerDoSEQ()

cntrl<-trainControl(method="CV",number=5)

rf.grid<-expand.grid(
	min.node.size=c(5,10,15,20),
	mtry=floor(length(c(vars.base,vars.covw1))*c(0.3,0.4,0.5,0.6,0.7)),
	splitrule="variance")

rf.rd.tune<-train(readtheta05~.,
	data=eclsb[,c("readtheta05",vars.base,vars.covw1)],
	method="ranger",
	tuneGrid=rf.grid,
	trControl=cntrl,
	metric="RMSE",
	respect.unordered.factors=TRUE,
	num.trees=ntrees,
	weights=eclsb$sampwt,
	seed=8675309)

rf.mt.tune<-train(maththeta05~.,
	data=eclsb[,c("maththeta05",vars.base,vars.covw1)],
	method="ranger",
	tuneGrid=rf.grid,
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
	min.node.size=rf.rd.tune$bestTune[,"min.node.size"],
	mtry=rf.rd.tune$bestTune[,"mtry"],
	splitrule="variance",
	respect.unordered.factors=TRUE,
	case.weights=wts.train,
	seed=8675309)

m1.mt<-ranger(
	y=math.train,
	x=x.train,
	num.trees=ntrees,
	min.node.size=rf.mt.tune$bestTune[,"min.node.size"],
	mtry=rf.mt.tune$bestTune[,"mtry"],
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

##### PLOT RESULTS #####
#load("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\_TEMP\\shap_valueplot_covars.RData")

theme_set(theme_bw())

numvar <- c("gender", "race", "birthwt", "momage01", "dadage01", "wic01", 
            "foodst01", "medicd01", "tanf01", "faminc01", "pared01", "momemp01", 
            "dademp01", "hhtotal01", "biodad01", "married01", "house01", 
            "rbooks01", "prmlang01", "region01", "urban01")
x.train[numvar] <- lapply(x.train[numvar], as.numeric)
colnames(shap.rd) <- colnames(shap.mt) <- colnames(x.train) <- var.label

shap_long1 <- shap.prep(shap_contrib = shap.rd, X_train=x.train)
plot1<-shap.plot.summary(shap_long1, x_bound = 1.5)+
  ggtitle("A. Reading Scores")
shap_long2 <- shap.prep(shap_contrib = shap.mt, X_train=x.train)
plot2<-shap.plot.summary(shap_long2, x_bound = 1.5)+
  ggtitle("B. Math Scores")

tiff("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\figures\\figure_S4.tiff",
	height=9,
	width=11,
	units='in',
	res=600)

grid.arrange(plot1, plot2, ncol=2)

dev.off()

sink("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\programs\\_LOGS\\22_create_figure_S4_log.txt")

print(startTime)
print(Sys.time())

sink()

save(shap.rd,shap.mt,file="C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\_TEMP\\shap_values_covars.RData")


