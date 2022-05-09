################################################
################################################
##                                            ##
## PROGRAM NAME: 10_create_figure_4           ##
## AUTHOR: GW                                 ##
## DESCRIPTION:                               ##
##                                            ##
##  create plot for marginal effects of       ##
##  nh poverty on toxic exposures using       ##
##                                            ##
################################################
################################################

rm(list=ls())

list.of.packages <- c(
	"foreach",
	"doParallel",
	"tidyverse",
	"forcats",
	"doRNG",
	"sys",
	"haven",
	"foreign",
	"dplyr",
	"tidyr",
	"ggplot2",
	"ranger",
	"caret",
	"SuperLearner",
	"dotwhisker",
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
	"region01",
	"urban01",
	"nhpopden01",
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

##### STANDARDIZE TOXICS ######
for (v in 1:length(vars.ntxw1)) {
	eclsb.mi[,vars.ntxw1[v]]<-(eclsb.mi[,vars.ntxw1[v]]-weighted.mean(eclsb.mi[,vars.ntxw1[v]],eclsb.mi$sampwt))/sqrt(wtd.var(eclsb.mi[,vars.ntxw1[v]],eclsb.mi$sampwt))
	}

##### SET HYPERPARAMETERS #####
ntrees<-200
min.node<-5
num.vars<-floor(length(c(vars.base,vars.covw1))/3)

##### PARALLELIZATION #####
my.cluster<-parallel::makeCluster(n.cores,type="PSOCK")
doParallel::registerDoParallel(cl=my.cluster)
foreach::getDoParRegistered()
clusterEvalQ(cl=my.cluster, library(ranger))
registerDoRNG(8675309)

##### ESTIMATE MARGINAL EFFECTS #####
miest.w1<-matrix(data=NA,nrow=length(vars.ntxw1),ncol=nmi)
boot.ci.mi<-NULL

for (i in 1:nmi) {

	### LOAD MI DATA ###
	eclsb<-eclsb.mi[which(eclsb.mi$minum==i),]
	print(c("nmi=", i))

	### COMPUTE ESTIMATES ###
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
			"num.vars",
			"min.node",
			"ntrees"),
		envir=environment())

	point.est.w1<-foreach(v=1:length(vars.ntxw1), .combine=cbind) %dopar% {

		y.train<-eclsb[,vars.ntxw1[v]]
		x.train<-eclsb[,c(vars.base,vars.covw1)]
		wts.train<-eclsb[,"sampwt"]

		m1.w1<-ranger(
			y=y.train,
			x=x.train,
			num.trees=ntrees,
			min.node.size=min.node,
			mtry=num.vars,
			splitrule="variance",
			respect.unordered.factors=TRUE,
			case.weights=wts.train,
			seed=8675309)

		x.imp<-eclsb[,c(vars.base,vars.covw1)]
		x.imp$nhpovrt01<-astar	
		uhat.astar<-predict(m1.w1,x.imp)$pred
		x.imp$nhpovrt01<-a	
		uhat.a<-predict(m1.w1,x.imp)$pred

		return(weighted.mean((uhat.astar-uhat.a),wts.train))
		}

	point.est.w1<-matrix(unlist(point.est.w1),ncol=1,byrow=TRUE)
	miest.w1[,i]<-point.est.w1[,1]

	### COMPUTE BOOTSTRAP CIs ###
	boot.ci<-NULL

	for (v in 1:length(vars.ntxw1)) {

		print(c("  var=", v))
		
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
				"num.vars",
				"min.node",
				"ntrees",
				"v"),
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
			
			boot.y.train<-boot.eclsb[,vars.ntxw1[v]]
			boot.x.train<-boot.eclsb[,c(vars.base,vars.covw1)]
			boot.wts.train<-boot.eclsb[,"sampwt"]
			
			boot.m1.w1<-ranger(
				y=boot.y.train,
				x=boot.x.train,
				num.trees=ntrees,
				min.node.size=min.node,
				splitrule="variance",
				respect.unordered.factors=TRUE,
				case.weights=boot.wts.train,
				seed=8675309)

			boot.x.imp<-boot.eclsb[,c(vars.base,vars.covw1)]
			boot.x.imp$nhpovrt01<-astar	
			boot.uhat.astar<-predict(boot.m1.w1,boot.x.imp)$pred
			boot.x.imp$nhpovrt01<-a	
			boot.uhat.a<-predict(boot.m1.w1,boot.x.imp)$pred
			boot.mrgnfx.w1<-weighted.mean((boot.uhat.astar-boot.uhat.a),boot.wts.train)

			return(list(boot.mrgnfx.w1))
			}

		boot.est.w1<-matrix(unlist(boot.est.w1),ncol=1,byrow=TRUE)
		
		boot.ci<-cbind(boot.ci,boot.est.w1[,1])
		}

	boot.ci.mi<-rbind(boot.ci.mi,boot.ci)
	}

stopCluster(my.cluster)
rm(my.cluster)

### COMBINE MI ESTIMATES ###
est.w1<-matrix(data=NA,nrow=length(vars.ntxw1),ncol=3)

for (i in 1:length(vars.ntxw1)) { 
	
	est.w1[i,1]<-round(mean(miest.w1[i,]),digits=3)
	est.w1[i,2]<-round(quantile(boot.ci.mi[,i],prob=0.025),digits=3)
	est.w1[i,3]<-round(quantile(boot.ci.mi[,i],prob=0.975),digits=3)
	}

### PRINT RESULTS ###
sink("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\programs\\_LOGS\\10_create_figure_4_log.txt")

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

output.w1<-as.data.frame(cbind(chem.label,est.w1))
colnames(output.w1)<-c('term','estimate','llci','ulci')
output.w1$estimate<-as.numeric(output.w1$estimate)
output.w1$llci<-as.numeric(output.w1$llci)
output.w1$ulci<-as.numeric(output.w1$ulci)
output.w1<-output.w1[order(output.w1$estimate,decreasing=T),]

print(output.w1)

##### CREATE DOT-WHISKER PLOT #####
#load("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\_TEMP\\tox_est.RData")

output.w1 %>%
	mutate(term=fct_reorder(term, estimate)) %>%
		ggplot(aes(x=term, y=estimate)) +
			geom_errorbar(aes(ymin=llci, ymax=ulci), width=0, color="gray55") +
			geom_point(size=1.05) +
			coord_flip() +
		 	theme_bw() +
		    	ggtitle("Estimated Marginal Effects of Neighborhood Poverty") +
			scale_y_continuous(
				name="Standard Deviations",
				limits=c(-0.50, 0.50),
				breaks=round(seq(-0.50,0.50,0.10), 2)) +
			scale_x_discrete(name=" ") +
		    	theme(plot.title=element_text(face="bold"),
      		      legend.position="none",
				panel.grid.major=element_blank(), 
				panel.grid.minor=element_blank(),
				panel.background=element_blank(), 
				axis.line=element_line(colour="black")) +
			geom_hline(yintercept=0, colour="black", linetype=2)

ggsave("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\figures\\figure_4.tiff", height=9, width=9)

print(startTime)
print(Sys.time())

sink()

save(output.w1,file="C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\_TEMP\\tox_est.RData")

registerDoSEQ()
