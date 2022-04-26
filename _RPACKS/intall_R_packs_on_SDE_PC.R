## ON NETWORKED COMPUTER:
# GET PACKAGE DEPENDENCIES
getDependencies <- function(packs) {
  dependencyNames <- unlist(
  tools::package_dependencies(packages = packs, 
    db = available.packages(), 
    which = c("Depends", "Imports"),
    recursive = TRUE))
  packageNames <- union(packs, dependencyNames)
  packageNames
  }
packages <- getDependencies(c("haven","ebal","dplyr","tidyr","CBPS","IDPmisc",
                              "ggplot2","survey","mgcv","foreign","devtools",
                              "ellipsis","MASS","nnet","glmnet","withr","Hmisc",
                              "whisker","boot","survival","geeM","geepack",
                              "maptools","rpart","xgboost","caret","BART","ranger",
                              "SuperLearner","dotwhisker","mice","RfEmpImp","glmnet",
					"randomForest","shapr","fastshap","doParallel",
					"ggrepel","Rcpp","metR","foreach","palmerpenguins",
					"tidyverse","kableExtra","doRNG","sys","e1071",
					"SHAPforxgboost","bartMachine","arsenal","rlang"))

# DOWNLOAD ALL PACKAGES
setwd("C:\\Users\\Geoffrey Wodtke\\Dropbox\\D\\projects\\nhood_mediation_toxins\\programs\\_RPACKS\\packages")
pkgInfo <- download.packages(pkgs = packages, destdir = getwd(), type = "win.binary")
write.csv(file = "pkgFilenames.csv", basename(pkgInfo[, 2]), row.names = FALSE)

## MOVE PACKAGES AND CSV FILE TO OFF-NETWORK COMPUTER

## ON OFF-NETWORK COMPUTER:
# READ PACKAGE FILE NAMES AND INSTALL FROM LOCAL DRIVE
setwd("C:\\Users\\wodtke\\Desktop\\projects\\nhood_mediation_toxins\\programs\\_RPACKS\\packages")
pkgFilenames <- read.csv("pkgFilenames.csv", stringsAsFactors = FALSE)[, 1]
install.packages(pkgFilenames, repos = NULL, type = "win.binary")
