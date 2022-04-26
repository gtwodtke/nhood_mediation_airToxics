################################################
################################################
##                                            ##
## PROGRAM NAME: 12_create_figure_5           ##
## AUTHOR: GW                                 ##
## DESCRIPTION:                               ##
##                                            ##
##  sensitivity analysis                      ##
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
	"gridExtra",
	"metR")

for(package.i in list.of.packages) {
	suppressPackageStartupMessages(library(package.i, character.only = TRUE))
	}

startTime<-Sys.time()

##### LOAD EFFECT ESTIMATES #####
load("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\data\\_TEMP\\efx_est.RData")

##### COMPUTE BIAS-ADJUSTED ESTIMATES #####
sens.grid<-expand.grid(v1=seq(-0.25,0.25,0.01),v2=seq(-0.25,0.25,0.01))
adj.grid<-cbind(sens.grid,
	nde.rd.adj=output.rd.w1["NDE","estimate"]-(sens.grid$v1*sens.grid$v2),
	nde.mt.adj=output.mt.w1["NDE","estimate"]-(sens.grid$v1*sens.grid$v2),
	nie.rd.adj=output.rd.w1["NIE","estimate"]+(sens.grid$v1*sens.grid$v2),
	nie.mt.adj=output.mt.w1["NIE","estimate"]+(sens.grid$v1*sens.grid$v2))

##### CONTOUR PLOTS #####
p.nde.rd<-ggplot(adj.grid, aes(x=v1,y=v2,z=nde.rd.adj,colour=stat(level))) +  
		geom_contour(breaks=seq(round(min(adj.grid$nde.rd.adj),2),round(max(adj.grid$nde.rd.adj),2),0.01),show.legend=FALSE) +
		scale_colour_distiller(palette="Greys",direction=-1) +
		xlab(expression(gamma)) +
		ylab(expression(eta)) +
		scale_x_continuous(breaks=seq(-0.25,0.25,0.05)) +
		scale_y_continuous(breaks=seq(-0.25,0.25,0.05)) +
		ggtitle("A. Bias-adjusted NDE on Reading Scores") +
		theme_bw(base_size=11) +
		theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
		geom_text_contour(
			breaks=seq(round(min(adj.grid$nde.rd.adj),2),round(max(adj.grid$nde.rd.adj),2),0.01),
			stroke=0.3,
			size=2.5,
			skip=0,
			color="black")

p.nde.mt<-ggplot(adj.grid, aes(x=v1,y=v2,z=nde.mt.adj,colour=stat(level))) +  
		geom_contour(breaks=seq(round(min(adj.grid$nde.mt.adj),2),round(max(adj.grid$nde.mt.adj),2),0.01),show.legend=FALSE) +
		scale_colour_distiller(palette="Greys",direction=-1) +
		xlab(expression(gamma)) +
		ylab(expression(eta)) +
		scale_x_continuous(breaks=seq(-0.25,0.25,0.05)) +
		scale_y_continuous(breaks=seq(-0.25,0.25,0.05)) +
		ggtitle("B. Bias-adjusted NDE on Math Scores") +
		theme_bw(base_size=11) +
		theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
		geom_text_contour(
			breaks=seq(round(min(adj.grid$nde.mt.adj),2),round(max(adj.grid$nde.mt.adj),2),0.01),
			stroke=0.3,
			size=2.5,
			skip=0,
			color="black")

p.nie.rd<-ggplot(adj.grid, aes(x=v1,y=v2,z=nie.rd.adj,colour=stat(level))) +  
		geom_contour(breaks=seq(round(min(adj.grid$nie.rd.adj),2),round(max(adj.grid$nie.rd.adj),2),0.01),show.legend=FALSE) +
		scale_colour_distiller(palette="Greys",direction=-1) +
		xlab(expression(gamma)) +
		ylab(expression(eta)) +
		scale_x_continuous(breaks=seq(-0.25,0.25,0.05)) +
		scale_y_continuous(breaks=seq(-0.25,0.25,0.05)) +
		ggtitle("C. Bias-adjusted NIE on Reading Scores") +
		theme_bw(base_size=11) +
		theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
		geom_text_contour(
			breaks=seq(round(min(adj.grid$nie.rd.adj),2),round(max(adj.grid$nie.rd.adj),2),0.01),
			stroke=0.3,
			size=2.5,
			skip=0,
			color="black")

p.nie.mt<-ggplot(adj.grid, aes(x=v1,y=v2,z=nie.mt.adj,colour=stat(level))) +  
		geom_contour(breaks=seq(round(min(adj.grid$nie.mt.adj),2),round(max(adj.grid$nie.mt.adj),2),0.01),show.legend=FALSE) +
		scale_colour_distiller(palette="Greys",direction=-1) +
		xlab(expression(gamma)) +
		ylab(expression(eta)) +
		scale_x_continuous(breaks=seq(-0.25,0.25,0.05)) +
		scale_y_continuous(breaks=seq(-0.25,0.25,0.05)) +
		ggtitle("D. Bias-adjusted NIE on Math Scores") +
		theme_bw(base_size=11) +
		theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
		geom_text_contour(
			breaks=seq(round(min(adj.grid$nie.mt.adj),2),round(max(adj.grid$nie.mt.adj),2),0.01),
			stroke=0.3,
			size=2.5,
			skip=0,
			color="black")

tiff("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\figures\\figure_5.tiff",
	height=9,
	width=9,
	units='in',
	res=600)

grid.arrange(
	p.nde.rd,p.nde.mt,
	p.nie.rd,p.nie.mt,
	ncol=2)

dev.off()

sink("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\programs\\_LOGS\\12_create_figure_5_log.txt")

print(startTime)
print(Sys.time())

sink()
