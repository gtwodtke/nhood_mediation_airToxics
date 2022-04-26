################################################
################################################
##                                            ##
## PROGRAM NAME: 21_create_figure_S3          ##
## AUTHOR: BP, GW                             ##
## DESCRIPTION:                               ##
##                                            ##
##  create SHAP summary plots for all vars    ##
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

##### COMPUTE SHAP VALUES #####
pfun <- function(object, newdata) { predict(object,data=newdata)$pred }

### LOAD 1ST IMPUTATION ###
eclsb<-eclsb.mi[which(eclsb.mi$minum==1),]

### TUNE RF HYPERPARAMETERS ###
registerDoSEQ()

cntrl<-trainControl(method="CV",number=5)

rf.grid<-expand.grid(
	min.node.size=c(5,10,15,20),
	mtry=floor(length(c(vars.base,vars.covw1,vars.ntxw1))*c(0.3,0.4,0.5,0.6,0.7)),
	splitrule="variance")

rf.rd.tune<-train(readtheta05~.,
	data=eclsb[,c("readtheta05",vars.base,vars.covw1,vars.ntxw1)],
	method="ranger",
	tuneGrid=rf.grid,
	trControl=cntrl,
	metric="RMSE",
	respect.unordered.factors=TRUE,
	num.trees=ntrees,
	weights=eclsb$sampwt,
	seed=8675309)

rf.mt.tune<-train(maththeta05~.,
	data=eclsb[,c("maththeta05",vars.base,vars.covw1,vars.ntxw1)],
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
x.train<-eclsb[,c(vars.base,vars.covw1,vars.ntxw1)]
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
#load("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\_TEMP\\shap_values_allvars.RData")

theme_set(theme_bw())

numvar <- c("gender", "race", "birthwt", "momage01", "dadage01", "wic01", 
            "foodst01", "medicd01", "tanf01", "faminc01", "pared01", "momemp01", 
            "dademp01", "hhtotal01", "biodad01", "married01", "house01", 
            "rbooks01", "prmlang01", "region01", "urban01")
x.train[numvar] <- lapply(x.train[numvar], as.numeric)
colnames(shap.rd) <- colnames(shap.mt) <- colnames(x.train) <- var.label

shap_long_rd <- shap.prep(shap_contrib = shap.rd, X_train=x.train)
shap_long_mt <- shap.prep(shap_contrib = shap.mt, X_train=x.train)

shap_long_rdH <- shap_long_rd[shap_long_rd$mean_value >= median(shap_long_rd$mean_value)]
shap_long_rdL <- shap_long_rd[shap_long_rd$mean_value < median(shap_long_rd$mean_value)]
shap_long_mtH <- shap_long_mt[shap_long_mt$mean_value >= median(shap_long_mt$mean_value)]
shap_long_mtL <- shap_long_mt[shap_long_mt$mean_value < median(shap_long_mt$mean_value)]

read_listH <- as.character(unique(shap_long_rdH$variable))
shap.rd.high <- subset(shap.rd, select = read_listH)
read_listL <- as.character(unique(shap_long_rdL$variable))
shap.rd.low <- subset(shap.rd, select = read_listL)
x.train_rdhigh <- x.train %>% select(one_of(dput(as.character(shap_long_rdH$variable))))
x.train_rdlow <- x.train %>% select(one_of(dput(as.character(shap_long_rdL$variable))))

math_listH <- as.character(unique(shap_long_mtH$variable))
shap.mt.high <- subset(shap.mt, select = math_listH)
math_listL <- as.character(unique(shap_long_mtL$variable))
shap.mt.low <- subset(shap.mt, select = math_listL)
x.train_mthigh <- x.train %>% select(one_of(dput(as.character(shap_long_mtH$variable))))
x.train_mtlow <- x.train %>% select(one_of(dput(as.character(shap_long_mtL$variable))))

SHAP_RDhigh <- shap.prep(shap_contrib = shap.rd.high, X_train=x.train_rdhigh)
SHAP_RDlow <- shap.prep(shap_contrib = shap.rd.low, X_train=x.train_rdlow)
SHAP_MThigh <- shap.prep(shap_contrib = shap.mt.high, X_train=x.train_mthigh)
SHAP_MTlow <- shap.prep(shap_contrib = shap.mt.low, X_train=x.train_mtlow)

plot1<-shap.plot.summary(SHAP_RDhigh, x_bound = 1)+
  ggtitle("A. Reading Scores")
plot2<-shap.plot.summary(SHAP_MThigh, x_bound = 1)+
  ggtitle("B. Math Scores")

tiff("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\figures\\figure_S3.1.tiff",
     height=10,
     width=11,
     units='in',
     res=600)

grid.arrange(plot1,plot2,ncol=2)

dev.off()

plot3<-shap.plot.summary(SHAP_RDlow, x_bound = 1)+
  ggtitle("A. Reading Scores")
plot4<-shap.plot.summary(SHAP_MTlow, x_bound = 1)+
  ggtitle("B. Math Scores")

tiff("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\figures\\figure_S3.2.tiff",
     height=10,
     width=11,
     units='in',
     res=600)

grid.arrange(plot3,plot4,ncol=2)

dev.off()

sink("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\programs\\_LOGS\\21_create_figure_S3_log.txt")

print(startTime)
print(Sys.time())

sink()

save(shap.rd,shap.mt,file="C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\_TEMP\\shap_values_allvars.RData")

