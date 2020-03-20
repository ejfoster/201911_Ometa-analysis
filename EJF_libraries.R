### Main libraries used 
### By: Erika Foster
### Created:20180917
### Updated:20191127

###  Objective: Place all of my commonly uploaded libraries in one source document

#x <- c("ggmap", "rgdal", "rgeos", "maptools", "dplyr", "tidyr", "tmap")
   # install.packages(x) # warning: uncommenting this may take a number of minutes
#lapply(x, library, character.only = TRUE) 

#install CRAN packages
#install.packages(c("doBy",  "stringr", "ggthemes", "plyr", lmerTest", 
#                  "corrplot", "gplots", "grid", "ape", "devtools", "rgl", "fdrtool", "qiimer", 
#                  "pheatmap","cluster","NbClust", "RColorBrewer", "corrgram","candisc", "data.table", 
#                   "car", "rgl", "stringi", "MTS", "Mass", "mixOmics", "httpuv",  "indicspecies", 
#                   "htmltools", "lsmeans", "e1071", 
#                   "bios2mds", ##needed for phyloseq
#                   "RVAideMemoire","rgdal", "raster", "HH", "jvm", "ggforce", "outliers",
#                   "metafor", "metaforest", "sp", "mmpf", "rgdal",
#                    "psych", "Cairo", "dichromat", "CCA","boot",
#                   "WriteXLS", "Rmisc", "GGally", "tidyverse", "dlookr", "broom"))  #for CAP pairwise.factorfit #HH for anocva

#install.packages(c("sf", "raster", "spData", "spDataLarge", "tmap", "leaflet", "mapview", "gstat"))
#structuring data, aggregating etc.
library(dlookr); library(broom)
library(raster); library(sp); library(mmpf); library(rgdal) #for get worldclim data
library(psych)
library(metafor); library(metaforest); 
library(WriteXLS)
library(plyr)
library(multcomp)#for anovas

library(tidyverse) #includes  [1] "broom"      "cli"        "crayon"     "dbplyr"     "dplyr"     
#>  [6] "forcats"    "ggplot2"    "haven"      "hms"        "httr"      
#> [11] "jsonlite"   "lubridate"  "magrittr"   "modelr"     "pillar"    
#> [16] "purrr"      "readr"      "readxl"     "reprex"     "rlang"     
#> [21] "rstudioapi" "rvest"      "stringr"    "tibble"     "tidyr"     
#> [26] "xml2"       "tidyverse" 
library(vegan) #for CCA, for bray curtis

library(lme4) #anova models
library(lmerTest) #changes anova ()
library(lsmeans)#anova models

library(MASS) # for NMS


#install.packages("sourcetools", type = "source")

#install.packages("outliers")
library(tidyr)
library(stringr)

library(outliers)
library(HH)
library(RVAideMemoire)
library(sourcetools) 
library(boot); library(cluster); library (foreign); library(lattice); library (Matrix); library(mgcv);
library(rpart); library(survival) 
library(doBy) #summaryBy function
library(Rmisc) #summarySE for summaries

library(dplyr) #wrangling data, select
library(plyr) #rename() 
library(Hmisc) #correlation
library(corrgram) #correlograms
library(corrplot); library(GGally)
library(psych) #for corr.test to work on matrix data

#graphing and exporting
library(reshape)
library(reshape2) #for "melt ()" ggplot2 graphing
library(Cairo) #exporting plots
 library(corrplot); library(pheatmap); #plotting
library(gplots) #fine tune some of the other packages
library(grid);library(ggthemes) #for graphics theme_few()
library(RColorBrewer)
library(dichromat)
library(ggforce)
library(dichromat)

#install.packages(c("stringi"),configure.args=c("--disable-cxx11"),repos="https://cran.rstudio.com")
#install.packages(c("stringr"),configure.args=c("--disable-cxx11"),repos="https://cran.rstudio.com")

library(stringi); 
library(CCA)
#load packages
library(car)
library(htmltools)

#all in tidyverse
#library(stringr); library(tidyr); library(stringr); library(ggplot2);

library(plyr) # for revalue of factor levels
#library(dplyr) #slect variables
library(lmerTest);  #for lmer mixed effects models
library(lsmeans)
library(emmeans) #for new verions of R 3.5
library(candisc)
library(MTS) # for mulitvariate linear model
library(MASS)

library(devtools)#install older versions of programs
library(rgl)# 3D ordination plots 
library(fdrtool) #false discovery rate (FDR) in data analysis
library(vegan)
library(GGally)


####################
# microbiome packages
###################
#install_version("vegan", version = "2.4-6", repos = "http://cran.us.r-project.org")

#library(phyloseq) #for ordinate() function
#library(ape) #make phylo trees
#library(apTreeshape)

#library(phangorn)
#library(qiimer) #https://cran.r-project.org/web/packages/qiimer/qiimer.pdf
#library(dada2)
#library(httpuv)
#library(mixOmics) #indicator species anlysis????
#library(indicspecies) #indicator species anlysis


####################
#Bioconductor - to install metagenomic functions
###################
   ### Use source to install bioclite.R which will allow you to load packages from Bioconductor
#source("http://bioconductor.org/biocLite.R")
   ###Run biocLite() to install dependencies when you install or when you update R.
#biocLite() #only UPdate some packages (SAY NO TO VEGAN UPDATES)
   ### Run biocLite to install the packages: phyloseq, metagenomSeq and DESeq2 
#biocLite("Rcpp")
#biocLite('phyloseq') #only UPdate some packages (SAY NO TO VEGAN UPDATES)
#biocLite('metagenomeSeq')
#biocLite('DESeq2')
#biocLite('dada2')
#biocLite("ggtree")
#biocLite("biomformat")

#library("phyloseq")
#library("metagenomeSeq")
#library("ggtree")
#library("biomformat")
#library("dada2")
#library("DESeq2")

### newest classifier "idTaxa"
#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DECIPHER")

#library(DECIPHER); packageVersion("DECIPHER")

####################
# Functions created by Erika 
####################

Emaxmin2<-function(df, y, var1, var2){ 
  varnames<-c("var1", "var2") # add others
  stats<-aggregate( y ~ var1+var2, data=df, FUN=function(x) c(max(x), (min(x))))
  head(stats)
  place.holder<- as.data.frame(stats[,ncol(stats)])
  stats.summary <-as.data.frame (cbind(stats[,1:ncol(stats)-1], place.holder))
  names(stats.summary)<-c(varnames, "max", "min")
  return(stats.summary)
}

Emeans2<-function(df, y, var1, var2){ 
  varnames<-c("var1", "var2") # add others
  stats<-aggregate( y ~ var1+var2, data=df, FUN=function(x) c(mean(x), (sd(x)/(sqrt(length(x))))))
  head(stats)
  place.holder<- as.data.frame(stats[,ncol(stats)])
  stats.summary <-as.data.frame (cbind(stats[,1:ncol(stats)-1], place.holder))
  names(stats.summary)<-c(varnames, "y.mean", "y.se")
  return(stats.summary)
}

Emeans4<-function(df, y, varname1, var1, varname2, var2, varname3, var3, varname4, var4){ 
  varnames<-c(varname1, varname2, varname3, varname4) # add others
  stats<-aggregate( y ~ var1+var2+var3+var4, data=df, FUN=function(x) c(mean(x), (sd(x)/(sqrt(length(x))))))
  head(stats)
  place.holder<- as.data.frame(stats[,ncol(stats)])
  stats.summary <-as.data.frame (cbind(stats[,1:ncol(stats)-1], place.holder))
  names(stats.summary)<-c(varnames, "y.mean", "y.se")
  return(stats.summary)
}

Emeans3<-function(df, y, varname1, var1, varname2, var2, varname3, var3){ 
  varnames<-c(varname1, varname2, varname3) # add others
  stats<-aggregate( y ~ var1+var2+var3, data=df, FUN=function(x) c(mean(x), (sd(x)/(sqrt(length(x))))))
  head(stats)
  place.holder<- as.data.frame(stats[,ncol(stats)])
  stats.summary <-as.data.frame (cbind(stats[,1:ncol(stats)-1], place.holder))
  names(stats.summary)<-c(varnames, "y.mean", "y.se")
  return(stats.summary)
}

Emeans2<-function(df, y, varname1, var1, varname2, var2){ 
  varnames<-c(varname1, varname2) # add others
  stats<-aggregate( y ~ var1+var2, data=df, FUN=function(x) c(mean(x), (sd(x)/(sqrt(length(x))))))
  head(stats)
  place.holder<- as.data.frame(stats[,ncol(stats)])
  stats.summary <-as.data.frame (cbind(stats[,1:ncol(stats)-1], place.holder))
  names(stats.summary)<-c(varnames, "y.mean", "y.se")
  return(stats.summary)
}

#pairwise comparisions for permanovas
pairwise.adonis <- function(x,factors, sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  co = combn(unique(factors),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(co[1,elem],co[2,elem]),] ~ factors[factors %in% c(co[1,elem],co[2,elem])] , method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}
###
