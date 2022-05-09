################################################
################################################
##                                            ##
## PROGRAM NAME: 14_create_figure_7           ##
## AUTHOR: GW                                 ##
## DESCRIPTION:                               ##
##                                            ##
##  create scatter plot of exposure-mediator  ##
##  and mediator-outcome SHAP importance      ##
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
	"ggrepel",
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

##### STANDARDIZE TOXINS ######
for (v in 1:length(vars.ntxw1)) {
	eclsb.mi[,vars.ntxw1[v]]<-(eclsb.mi[,vars.ntxw1[v]]-weighted.mean(eclsb.mi[,vars.ntxw1[v]],eclsb.mi$sampwt))/sqrt(wtd.var(eclsb.mi[,vars.ntxw1[v]],eclsb.mi$sampwt))
	}

##### SET RF HYPERPARAMETERS #####
ntrees<-200
min.node<-5
num.vars<-floor(length(c(vars.base,vars.covw1))/3)

##### ESTIMATE A->M SHAP IMPORTANCE #####
miest<-matrix(data=NA,nrow=length(vars.ntxw1),ncol=nmi)

pfun <- function(object, newdata) {	predict(object,data=newdata)$pred }

for (i in 1:nmi) {

	### LOAD MI DATA ###
	eclsb<-eclsb.mi[which(eclsb.mi$minum==i),]

	for (v in 1:length(vars.ntxw1)) {

		print(c("nmi=",i,"/","nvar=",v))

		### FIT MODELS ###
		y.train<-eclsb[,vars.ntxw1[v]]
		x.train<-eclsb[,c(vars.base,vars.covw1)]
		wts.train<-eclsb[,"sampwt"]
		
		m1<-ranger(
			y=y.train,
			x=x.train,
			num.trees=ntrees,
			min.node.size=min.node,
			mtry=num.vars,
			splitrule="variance",
			respect.unordered.factors=TRUE,
			case.weights=wts.train,
			seed=8675309)

		### SETUP PARALLEL PROCESSING ###
		cl<-makeCluster(n.cores,type="PSOCK")
		registerDoParallel(cl)
		clusterEvalQ(cl,library(ranger))
		clusterExport(cl,c("eclsb","shap.sim","pfun","m1","x.train"),envir=environment())
		registerDoRNG(8675309)

		### COMPUTE SHAPLEY VALUES ###
		shap<-explain(m1,X=x.train, pred_wrapper=pfun, nsim=shap.sim, .parallel=TRUE)

		stopCluster(cl)
		rm(cl)

		### COMPUTE SHAP IMPORTANCE ###
		shap.nh<-as.data.frame(shap["nhpovrt01"])
		miest[v,i]<-mean(abs(shap.nh[,"nhpovrt01"]))
		}
	}

### COMBINE MI ESTIMATES ###
est.nh.w1<-matrix(data=NA,nrow=length(vars.ntxw1),ncol=2)

for (i in 1:length(vars.ntxw1)) { 

	est.nh.w1[i,1]<-chem.label[i]
	est.nh.w1[i,2]<-round(mean(miest[i,]),digits=4)
	}

### PRINT RESULTS ###
sink("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\programs\\_LOGS\\14_create_figure_7_log.txt")

load("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\_TEMP\\shap_est.RData")

output.nh.w1<-as.data.frame(est.nh.w1)

colnames(output.nh.w1)<-c('label','importance')
rownames(output.nh.w1)<-vars.ntxw1

output.nh.w1$importance<-as.numeric(output.nh.w1$importance)

output.nh.w1<-output.nh.w1[order(output.nh.w1$label,decreasing=F),]
output.rd.w1<-output.rd.w1[order(output.rd.w1$label,decreasing=F),]
output.mt.w1<-output.mt.w1[order(output.mt.w1$label,decreasing=F),]

output<-cbind(output.nh.w1,output.rd.w1[,2],output.mt.w1[,2])
colnames(output)<-c('label','importance.nh','importance.rd','importance.mt')

print(output)

##### CREATE SCATTER PLOT #####
#load("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\_TEMP\\shap_est_combined.RData")

p1.base<-ggplot(output, aes(importance.nh,importance.rd,label=label))+
	xlab(bquote(poverty%->%toxin ~ "mean(|SHAP value|)")) +
	ylab(bquote(toxin%->%ability ~ "mean(|SHAP value|)")) +
	ggtitle("A. Reading Scores") +
	scale_y_continuous(breaks=seq(0,0.018,0.003),limits=c(0,0.018)) +
	scale_x_continuous(breaks=seq(0,0.18,0.03),limits=c(0,0.18)) +
	geom_point(color="black",size=0.6)+
	theme_grey(base_size=8.5)

p1.adj<-p1.base + geom_text_repel(size=1.85,max.overlaps=10)

p2.base<-ggplot(output, aes(importance.nh,importance.mt,label=label))+
	xlab(bquote(poverty%->%toxin ~ "mean(|SHAP value|)")) +
	ylab(bquote(toxin%->%ability ~ "mean(|SHAP value|)")) +
	ggtitle("B. Math Scores") +
	scale_y_continuous(breaks=seq(0,0.018,0.003),limits=c(0,0.018)) +
	scale_x_continuous(breaks=seq(0,0.18,0.03),limits=c(0,0.18)) +
	geom_point(color="black",size=0.6)+
	theme_grey(base_size=8.5)

p2.adj<-p2.base + geom_text_repel(size=1.85,max.overlaps=10)

tiff("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\figures\\figure_7.tiff",
	height=4.5,
	width=9,
	units='in',
	res=600)

grid.arrange(p1.adj,p2.adj,ncol=2)

dev.off()

print(startTime)
print(Sys.time())

sink()

save(output,file="C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\_TEMP\\shap_est_combined.RData")

