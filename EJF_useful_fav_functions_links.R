########################################
# My saved functions and R resources
########################################

### By: Erika Foster
### Created:20180917
### Updated:20191127
     
# Still need to clean up all of these old additions and categorize more logically, remove repeats

########################################
# Online Resources
########################################
#basic resource: Datacamp.com

#begining bayes bayesian stats in R
#https://www.datacamp.com/community/open-courses/beginning-bayes-in-r

#learning R
library(adventr)
learnr::run_tutorial("name_of_tutorial", package = "adventr")
#for chapter 3
learnr::run_tutorial("adventr_03", package = "adventr")

#regression moderndive package
#https://moderndive.com/6-multiple-regression.html
#confidence intervals,using infer package
#https://moderndive.com/8-confidence-intervals.html#ci-build-up

### Intro to R
#wrangling, graphing: https://moderndive.com/3-wrangling.html#filter
#tidy data: https://moderndive.com/4-tidy.html
#summarise: http://rpubs.com/justmarkham/dplyr-tutorial

### Git Hub introduction
#basic functions / pull requests: https://guides.github.com/activities/hello-world/
# connect R and GitHub https://www.youtube.com/watch?v=-c2uNqEE6-c
# get older versions on GitHub: https://www.infoworld.com/video/97367/r-tip-how-to-use-git-and-github-with-r-projects time 6:25
#happygitwithR.com #ch4, 6, 7 for setup

########################################
# Summary statistics
########################################

#summarize data: 
summarize(data, mean_bty_avg = mean(bty_avg), mean_score = mean(score),
          median_bty_avg = median(bty_avg), median_score = median(score)) 
#can use 
flights%>%
  summarize_each(funs(min(..na.rm=TRUE), max(..na.rm=TRUE)),
                 match("Delay")) #n=
                 
flights %>%
  group_by(Month, DayofMonth) %>%
  tally(sort=TRUE)
                 
flights%>%
  groupby(Carrier) >%>
  summarise_each(funs(mean), Cancelled, Diverted) #or 2:4

########################################
# Extracting climate data from lat/long
########################################               
                 
# extract climate data MAP MAT mean annual precipitation mean annual temperature
#http://www.worldclim.org/formats1
                 
#tmean , prec ; 12 data layers 1 for each month
                 #BIO1 = Annual Mean Temperature; BIO2 = Mean Diurnal Range (Mean of monthly (max temp – min temp)); BIO3 = Isothermality (BIO2/BIO7) (* 100)
                 #BIO4 = Temperature Seasonality (standard deviation *100); BIO5 = Max Temperature of Warmest Month, BIO6 = Min Temperature of Coldest Month
                 #BIO7 = Temperature Annual Range (BIO5-BIO6), BIO8 = Mean Temperature of Wettest Quarter, BIO9 = Mean Temperature of Driest Quarter
                 #BIO10 = Mean Temperature of Warmest Quarter, BIO11 = Mean Temperature of Coldest Quarter, BIO12 = Annual Precipitation
                 #BIO13 = Precipitation of Wettest Month, BIO14 = Precipitation of Driest Month, BIO15 = Precipitation Seasonality (Coefficient of Variation)
                 #BIO16 = Precipitation of Wettest Quarter, BIO17 = Precipitation of Driest Quarter, BIO18 = Precipitation of Warmest Quarter, BIO19 = Precipitation of Coldest Quarter
                 
#download resolution of climate data desired: http://worldclim.org/version2 #version 2 = 1970-2000
#or version 1.4, under generic grid format, downlaod bioclim 2.5 from https://www.worldclim.org/current
                 
r <- raster::getData("worldclim",var="bio",res=2.5)
                 
                 
# also possible to get future climate data: 
# https://rdrr.io/cran/raster/man/getData.html
                 
#r <- r[[c(1,12)]] #Bio1 and Bio12 selected
#names(r) <- c("Temp","Prec")
                 
                 #lats <- c(df1$lat.dec.degree) #use spTransform if not in WGS 84 lat/lon (EPSG 4326)
                 #longs <- c(df1$long.dec.degree) 
                 # coords <- data.frame(x=longs,y=lats)
                 
                 # points <- SpatialPoints(coords, proj4string = r@crs)
                 # values <- raster::extract(r,points)
                 
                 #df.clim <- cbind.data.frame(coordinates(points),values); head(df.clim,2)
                 #df.clim$Temp<-df.clim$Temp/10 #WorldCLim data hs a scale factor of 10 for temp
                 #WriteXLS(df.clim, "ometa_climate.xls") #specific file here
                 
########################################
# Start with intall
########################################

#how to install a package
install.packages("Rcmdr")
                 
#install multiple packages at once
wants <- c("AICcmodavg", "lme4", "multcomp", "nlme", "pbkrtest")
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
                 
#setting the working directory
setwd("C:/Users/erika/Dropbox/1 GRADUATE SCHOOL/CSU Classes/2014 Fall Classes/ECOL 592 MultiVariate Statistics")
                 
#load a data set by choosing the file
flux.data<-read.csv(file.choose(),header=T)
                 
#information on R version/session
sessionInfo()
                 
#install latest version of R
  # installing/loading the package:
if(!require(installr)) {
install.packages("installr"); require(installr)} # load / install+load installr
# using the package:
updateR() # this will start the updating process of your R installation.
                 
#Clear console or clear environment, remove list
rm(list = ls()) #clear environment
                 
#reinstall packages install.packages
                 
install.packages(c("knitr", "kableExtra", "gplots", "ggplot2", "tidyr", "multcomp",
               "vegan", "lme4", "lmerTest", "lsmeans", "reshape2", "outliers",
                "Hmisc", "tables", "stargazer", "plyr", "Cairo", "pbkrtest", "geoR",
                "multcompView", "dplyr", "knitcitations"))

########################################
# Basic R functions
########################################

#creating a sequence of numbers, 0 to 100, by 1s
sequence.name<-seq(0,100,1)
combined.name<-cbind(data.name1,data.sd)

#look at means per month, using data, dependent variable (DV) and function mean
with(N.data, tapply(nitrate.DV, month.trt, mean))

#Calculate mean or a function for each column by column
#which applies the function passed as argument to each column of the dataframe.
sapply(training.data.raw, function(x) sum(is.na(x)))

#other
#student t distribution
dt()
#Example for one-sided -> pt(1.708, df=25, lower.tail=FALSE)
pt()
qt()

#t-test
t.test(Lead$Lead,alternative=c("two.sided"),mu=30)
rt()

#Check to see what is your workspace
ls()

########################################
# Exploring Data Iniitially
########################################
                 str(flux.data)
                 head(flux.data); dim(flux.data);length(flux.data); names(flux.data)
                 ncol(data1) #to get the number of columns in the data set
                 
                 #what data type is an object
                 mode(flux.data$trt)
                 
                 #to make factors into integers
                 as.integer(factor(mydata$column.name,c("B","M","C")))
                 
                 #to rename a column etc
                 names(biomass.df)[5:6]<- c("se", "mean") 
                 
                 #Change NULL plots to CONTROL =  change label of a factor = rename a factor
                 levels(data1.CN$trt)[levels(data1.CN$trt)=="N"] <- "C"
                 
                 #take make negative numbers to zeros
                 data1[,8:(ncol(data1))][data1[,8:(ncol(data1))]<0]=0 #change all negatives to zero, col 8:final
                 data1
                 
########################################
# Remove outliers  or drop values 
########################################
                 
                 #drop average ee.C and ee.N subsets
                 EEAratio2<-subset(EEAratio1, variable != "ee.N")
                 EEAratio3<-subset(EEAratio2, variable != "ee.C")
                 #tail(EEAratio3,5)
                 
                 #select specific treatments, print only columns weight through income
                 newdata <- subset(mydata, sex=="m" & age > 25,
                                   select=weight:income)
                 
                 #select ages under 10 and greater than or equal to 20, with ID and weight only
                 newdata <- subset(mydata, age >= 20 | age < 10, select=c(ID, Weight))
                 
                 #need to remove negative variables from dataset first
                 head(subset.M2,3)
                 subset.B2[,13][subset.B2[,13]<0]=NaN #remove values from column 13
                 subset.M2[,13][subset.B2[,13]<0]=NaN 
                 
                 #Remove outliers based on numerical value
                 soilP2<-subset(subP1, value <60,na.omit=T);soilP2
                 data.mbn.corr$value[data.mbn.corr$value > 37] <- NA #change to NA
                 
                 #need to remove outliers from residual plot
                 nrow(ratio1.june)
                 ratio2.june <- ratio1.june$[ratio < 20]
                 nrow(EEAratio)
                 
                 #check for outliers and remove
                 #Outliers? = Remove #Remove Outlier (enzymes: #33june,#10june)
                 outlier1.cb <-outlier(subset1.cb$value,logical=TRUE)
                 find_outlier1.cb <-which(outlier1.cb==TRUE,arr.ind=TRUE)
                 subset1.cb[find_outlier1.cb,]
                 subset.cb <-subset1.cb[-c(find_outlier1.cb),]
                 str(subset.cb)
                 
                 #Drop levels that are no longer present
                 levels(droplevels(subset1.cb$irrig));levels(droplevels(subset1.cb$variable))
                 levels(droplevels(subset1.cb$date)); levels(droplevels(subset1.cb$trt))
                 str(subset1.cb)
                 value1.cb<-subset1.cb$value
                 
                 #Remove a variable, this code leaves NAs at the bottom of the data frame... work on this?
                 subset.noD_noBM<-subset.noD[ ! subset1.no3$trt %in% c("BM"), ]
                 
                 #To remove or eliminate empty rows or columns
                 soil1<-soil1[,1:ncol(data1)]
                 

                 
########################################
# Checking for Normality
########################################
                 
                 #histogram
                 hist(data)
                 hist(data.name,breaks=seq(8,10,0.1),xlab="Title X Axis",main="Title Y Axis")
                 
                 #qqplot
                 qqnorm(data); qqline(data)  
                 
                 #Shapiro-Wilk test for normality
                 shapiro.test(data) #to test for normality
                 
                 #check boxplots for quick analysis
                 par(mfrow=(c(1,3)))
                 plot(no3.logvalue~trt,data=subset.no3, main="yield by trt bu/ac")
                 plot(no3.logvalue~irrig,data=subset.no3, main="yield by irrig bu/ac")
                 plot(no3.logvalue~corn,data=subset.no3, main="no3 by corn bu/ac")
                 
                 #check resdual v fitted values (for equality of variance)
                 
                 ###
                 ###Transformations
                 ###
                 
                 #outliers, remove values over a certain level
                 mySmallNUms <- myNums[myNums <= 523.689]
                 
                 transform.log<-log(data$column)
                 
                 #Change negative values to zeros, for all of the response variables 13:32
                 crop1[,13:((ncol(crop1)))][(crop1[,13:(ncol(crop1))] <0)] <- 0 
                 
                 #replace value with NA
                 df[df == 0] <- NA
                 
                 #replae values with corrected values, apply function to a subset
                 data.mbc2[data.mbc2$trt=="B",]$value<-(data.mbc2[data.mbc2$trt=="B",]$value)/.42
                 data.mbc2
                 data.mbc.corr<-cbind(data.mbc, data.mbc2$value); names(data.mbc.corr)[12]<- c("value.corr") 
                 data.mbc.corr
                 
                 
                 #Nas and negatives to 0, check to see if na's exist
                 is.na(x)
                 #convert na's in data columns 2 and 3
                 data[,2:3][is.na(data[ , 2:3] ) ] = 0 
                 
                 #recreates data x without missing na values as data y
                 data.y <- data.x[!is.na(data.x)]
                 
                 # data x places z in the values of x+1 for all non-missing nas and positive x's
                 (x+1)[(!is.na(x)) & x>0] -> z
                 
                 #change, relabel, rename substitute a factor variable name to something else
                 junk$nm <- as.character(junk$nm)
                 junk$nm[junk$nm == "B"] <- "b" #This didn't work for me last time, try the two below
                 

                 
                 levels(flux.df$trt)[levels(flux.df$trt)=="B"] <- "Biochar"
                 
                 #subsets
                 treatment1=subset(trt,Trt=="Manure")
                 treatment0=subset(trt,Trt=="Control")
                 meadian(treatment1);median(treatment2)
                 combined.name<-cbind(treatment1,treatment2)
                 #OR
                 treatment0 <- rat.data[rat.data$Trt=="Ctrl",]
                 low=subset(rat.data,Trt=="LowChr");low
                 low <- rat.data[rat.data$Trt=="LowChr",]
                 boxplot(control$Enzyme,low$Enzyme, main="Rat Diets Low-Cromium vs Control", xlab="Control")
                 
                 #Convert outliers to NA s ? outlier to na
                 outlier1<-outlier(crop$cob,logical=TRUE)
                 find_outlier1 <-which(outlier1==TRUE,arr.ind=TRUE)
                 find_outlier1
                 crop[find_outlier1,]<-NA 
          
                 #Outlier to NA function
                 Eoutlier<-function(df, y, n.outliers){
                   myColnumber<-as.character(y) #save name of column to match later
                   for(i in 1:n.outliers){
                     outlier[i] <-outlier(df$y,logical=TRUE)
                     find_outlier[i] <-which(outlier[i]==TRUE,arr.ind=TRUE)
                     df[c(find_outlier), match(myColnumber,names(df))]<-NA  #match column name to a number
                   }
                   df.new<-df
                   return(outlier, df.new)
                 } 
                 
                 df<-sample.df
                 
                 Eoutlier(sample.df, sample.df$ProteinSorbedPC)

                 
########################################
# User Defined Functions (My Functions)
########################################
                 Esums<-function(df, y, var1, var2){ 
                   varnames<-c("var1", "var2") # add others
                   stats<-aggregate( y ~ var1+var2, data=df, FUN=function(x) c(sum(x), (sd(x)/(sqrt(length(x))))))
                   head(stats)
                   place.holder<- as.data.frame(stats[,ncol(stats)])
                   stats.summary <-as.data.frame (cbind(stats[,1:ncol(stats)-1], place.holder))
                   names(stats.summary)<-c(varnames, "y.sum", "y.se")
                   return(stats.summary)
                 }
                 
                 Emeans<-function(df, y, var1, var2){ 
                   varnames<-c("var1", "var2") # add others
                   stats<-aggregate( y ~ var1+var2, data=df, FUN=function(x) c(mean(x), (sd(x)/(sqrt(length(x))))))
                   head(stats)
                   place.holder<- as.data.frame(stats[,ncol(stats)])
                   stats.summary <-as.data.frame (cbind(stats[,1:ncol(stats)-1], place.holder))
                   names(stats.summary)<-c(varnames, "y.mean", "y.se")
                   return(stats.summary)
                 }
                 
            
                 

                 #####Functions to Play with#####
                 #substitutions? with gsub???
                 plot(sapply(dataliste, function(x)gsub(",", ".", x)))
                 
                 ##from AVA
                 #my main data.frame is called "DD" 
                 #my categories are "trt" and "geno" (you'd put your corn type, irrigation, date etc column headings
                 #na.rm=T gets rid of missing values
                 #FUN can be changed to sd, etc.
                 #print column of choice if you have a shitload of columns to wade through
                 
                 aggdata1=aggregate(DD,by=list(trt,geno),FUN=mean,na.rm=T);
                 print(aggdata1[54:56])
                 
                 
                 #########################################################################################################
                 ###
                 ### Graphing
                 ###
                 #####################################################################################################
                 
                 #use ggplot stat_smmary to create graph of mean and sd
                 ggplot(ToothGrowth, aes(y = len, x = supp, colour = dose, group = dose)) + 
                   stat_summary(fun.y = mean,
                                fun.ymin = function(x) mean(x) - sd(x), 
                                fun.ymax = function(x) mean(x) + sd(x), 
                                geom = "pointrange") +
                   stat_summary(fun.y = mean,
                                geom = "line") +
                   facet_wrap( ~ F3))



#Place superscript on axis title
#y = ("Potential Enzyme Activity (nmol" ~ g dry soil^{-1} ~ hr^{-1}  ")" )

#this works for sure:
labs( x="Irrigation Level", 
      y= (expression("Gravimetric Water Content (g water" ~ g^{-1} ~ "dry soil)" )),
      fill=("Soil\nAmendment"))

### ggplot2
set.seed(1234)
x <- replicate(8, round(10 * rexp(2000, 10)))
y <- apply(x, 2, function(column) table(factor(column, levels = 0:9)))
y <- as.data.frame(y)
colnames(y) <- paste('A', seq(1,ncol(y),1), sep='')
rownames(y) <- paste('R', seq(1,nrow(y),1), sep='')

#theme(axis.text.x = element_text(angle=90, vjust=1) #change angle of xaxis labels

library(ggplot2)
library(reshape)
y$ID <- rownames(y)
y.melt <- melt(y, id.var = 'ID')

y.melt <- within(y.melt, ID <- factor(ID, 
                                      c('R10','R9','R8','R7','R6','R5','R4','R3','R2','R1'), 
                                      ordered = TRUE))

ggplot(y.melt, aes(x = variable, y = value, fill = ID)) +
  geom_bar(stat = 'identity') +
  xlab("") +
  ylab("My variable") +
  theme(legend.title=element_blank())

#Point shapes for ggplot2
sum <- ggplot()
for(i in 1:20) { 
  sum <- sum + 
    geom_point(data=data.frame(x=c(i)),aes(x=x,y=x),shape=i)
}
sum

?seq( )
dev.off() #helps with "Error in .Call.graphics(C_palette2)...invalid graphics state

# Color chart
# http://research.stowers-institute.org/efg/R/Color/Chart/
# Color for ggplot hexadecimicals: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/hextable.png
###########################################################################################
###
###AIC analysis
###
###########################################################################################

library(AICcmodavg)
AICc(fitF)
[1] 2304.445
aictab(cand.set=list(fitR, fitF),
       modnames=c("restricted", "full"),
       sort=FALSE, second.ord=FALSE)

###################################
# Saved CIG Code - Summer 2014
####################################

#Code from Ann Hess: We can account for repeated measures using:
lmer (Y ~ Irr*Soil*Var*Time + (1|Block) + (1|Block:Irr)
      + (1|Block:Irr:Soil) + (1|Block:Irr:Var) + (1|Block:Irr:Soil:Var))

#ie
lmer (no3 ~ irrig*trt*corn*date + (1|block) + (1|block:irrig)
      + (1|block:irrig:corn) + (1|block:irrig:trt) + (1|block:irrig:corn:trt))


#Lsmeans, slice by a specific factor or level
contrasts.date.no3<-lsmeans(lmer.no3, pairwise ~ date, 
                            at = list(price1 = 50,
                                      price2 = c(40,60), day = c("2","3","4")) );contrasts.date.no3


### NO Assumption of compound symmetry
anova(lme(Y ~ Xb1*Xb2*Xw1*Xw2,
          random=list(id=pdBlocked(list(~1, pdIdent(~Xw1-1), pdIdent(~Xw2-1)))),
          method="ML", data=d2))




###
### Website Links
###

#Code on publication quality tables and graphs in html form:
#http://www.statmethods.net/interface/output.html

#Using Sweave and knitr:
#https://support.rstudio.com/hc/en-us/articles/200552056-Using-Sweave-and-knitr
install(knitr)

#Notes on LaTex combatibility
#eliminates copy and pasting all analyses
library(knitr)

\documentclass[11]{article}

\title{R-no-web extentsion knitr test}
\author(Erika Foster)

\begin {document}

\maketitle

\begin{abstract}
\end{abstract}

\section{}
testing whether tis compiles correctly
\end {document}

###################################################################################################
###
### Markdown Notes
###
###############################################################################################

###Notes for Markdown:

#show code but not results

#If I want this to appear more like a typewriter style use this notation: `finalprojectdata`
#>*Create a block quoate, with stars to italisze*
#  *Bullet 1
#*to tab over the bullets
#*just tab over a few times
#*Bullet 2
#**Bold a word**

#Tables in markdown

> library(knitr)
> kable(head(iris[,1:3]), format = "markdown")
Significance:
  |  Contrast    |Signficance  |  TEST        |
  |-------------:|------------:|-------------:|
  |           5,1|          3,5|           1,4|
  |           4,9|          3,0|           1,4|
  |           4,7|          3,2|           1,3|
  |           4,6|          3,1|           1,5|
  |           5,0|          3,6|           1,4|
  |           5,4|          3,9|           1,7|
  
  library(ascii)
print(ascii(head(iris[,1:3])), type = 'pandoc')

**Sepal.Length**   **Sepal.Width**   **Petal.Length**  
  --- ------------------ ----------------- ------------------
  1   5.10               3.50              1.40              
2   4.90               3.00              1.40              
3   4.70               3.20              1.30              
4   4.60               3.10              1.50              
5   5.00               3.60              1.40              
6   5.40               3.90              1.70              
--- ------------------ ----------------- ------------------


#######################################################################################
####
### Contrasts matrix to play with for lsmeans
###
########################################################################################

####Extract p-value from the lsmeans table of contrasts that are <=.05
contrasts.trt.irrig.vwc<-lsmeans(lmer.vwc, pairwise ~ trt+irrig); contrasts.trt.irrig.vwc
#create table of significant values
vwc.table<-summary(contrasts.trt.irrig.vwc$contrasts)
vwc.table[(vwc.table$p.value<=.05),]

#print only signifcant p-values less than .05
splt <- strsplit(as.character(anova.gravwc2$"Pr(>F)"), ",")
rows <- unlist(lapply(splt, function(a) all(as.numeric(a) <= .05))) #insert alpha value here (currently set to .05)
anova.gravwc2[rows,c(1,6)] #leave off the column section if full ANOVA table desired


?contrasts
contrasts.test<-cbind(c(1, 0, -1, 0, 0, 0),
                      c(1, 0, 0, 0, -1, 0),
                      c(0, 1, 0, -1, 0, 0),
                      c(0, 1, 0, 0, 0, -1),
                      c(0, 0, 1, 0, -1, 0),
                      c(0, 0, 0, 1, 0, -1)); contrasts.test
contrasts(contrasts.test)

###################
###
###Resdiual and Transformation Code for unequal variance over time
###
#########################
# use residuals to check for variance between factors
plot(resid.soilC.log,subset.soilC$date) #variance between dates looks even, maintain as a fixed effect


#trans.mbc<-sqrt(value.mbc^2)
#subset.mbc<-cbind(mbc.data,trans.mbc)
#subset.mbc<- subset.mbc[! subset.mbc$sqrt.mbc %in% c("-Inf"), ]

#boxcox transformation
#boxcox(mbc ~ irrig*trt*corn + (1|date)+(1|block) +(1|block:irrig)+ (1|block:irrig:corn) + (1|block:irrig:trt) + (1|block:irrig:date:trt),data=subset.mbc)

#plot(resid.mbc)

#Lots of outliers
### start with 'gravwc.data'
outlier1.gravwc <-outlier(gravwc.data$value,logical=TRUE)
find_outlier1.gravwc <-which(outlier1.gravwc==TRUE,arr.ind=TRUE)
gravwc.data[find_outlier1.gravwc,]
subset.gravwc1<-gravwc.data[-c(find_outlier1.gravwc),]
str(subset.gravwc1)

#check normality 1
qqnorm(subset.gravwc1$value);qqline(subset.gravwc1$value)
hist(subset.gravwc1$value)# right skew
shapiro.test(subset.gravwc1$value) # 

### start with 'subset.gravwc1'
outlier3.gravwc <-outlier(subset.gravwc1$value,logical=TRUE)
find_outlier3.gravwc <-which(outlier3.gravwc==TRUE,arr.ind=TRUE)
subset.gravwc1[find_outlier3.gravwc,]
subset.gravwc2<-subset.gravwc1[-c(find_outlier3.gravwc),]
str(subset.gravwc2)

#check normality 2
qqnorm(subset.gravwc2$value);qqline(subset.gravwc2$value)
shapiro.test(subset.gravwc2$value) # 

### start with 'subset.gravwc2'
outlier3.gravwc <-outlier(subset.gravwc2$value,logical=TRUE)
find_outlier3.gravwc <-which(outlier3.gravwc==TRUE,arr.ind=TRUE)
subset.gravwc2[find_outlier3.gravwc,]
subset.gravwc3<-subset.gravwc2[-c(find_outlier3.gravwc),]
str(subset.gravwc3)

#check normality 3
qqnorm(subset.gravwc3$value);qqline(subset.gravwc2$value)
shapiro.test(subset.gravwc3$value) 

### start with 'subset.gravwc3'
outlier4.gravwc <-outlier(subset.gravwc3$value,logical=TRUE)
find_outlier4.gravwc <-which(outlier4.gravwc==TRUE,arr.ind=TRUE)
subset.gravwc3[find_outlier4.gravwc,]
subset.gravwc4<-subset.gravwc3[-c(find_outlier4.gravwc),]
str(subset.gravwc4)

#check normality 4
qqnorm(subset.gravwc4$value);qqline(subset.gravwc3$value)
shapiro.test(subset.gravwc4$value) 

### start with 'subset.gravwc4'
outlier5.gravwc <-outlier(subset.gravwc4$value,logical=TRUE)
find_outlier5.gravwc <-which(outlier5.gravwc==TRUE,arr.ind=TRUE)
subset.gravwc4[find_outlier5.gravwc,]
subset.gravwc5<-subset.gravwc4[-c(find_outlier5.gravwc),]
str(subset.gravwc5)

#check normality 5
qqnorm(subset.gravwc5$value);qqline(subset.gravwc5$value)
shapiro.test(subset.gravwc5$value)  #p-value .1422

### start with 'subset.gravwc5'
outlier6.gravwc <-outlier(subset.gravwc5$value,logical=TRUE)
find_outlier6.gravwc <-which(outlier5.gravwc==TRUE,arr.ind=TRUE)
subset.gravwc5[find_outlier6.gravwc,]
subset.gravwc6<-subset.gravwc5[-c(find_outlier6.gravwc),]
str(subset.gravwc6)

#check normality 5
qqnorm(subset.gravwc6$value);qqline(subset.gravwc6$value)
shapiro.test(subset.gravwc6$value)  #p-value .1422 = just as good keeping that 6th outlier in
#REMOVED 5 outliers

#Relable response variable to use
gravwc<-subset.gravwc6$value
subset.gravwc6$irrig<-droplevels(subset.gravwc6$irrig)
length(subset.gravwc6$irrig="F")

#create a subset w/out BM, to drop a factor level
flux.df3<-subset(flux.df2, trt !="BM")

###FROM MULTIVARIATE STATISTICS CLASS:

# replace NAs with zeros
seeds_s[is.na(seeds_s)] <- 0  
head(seeds_s)



###########################################################
## Part 10: Yield
###########################################################
####A. Preliminary examination of dataset, checking normality
# {r echo=FALSE}
head(df.crop1)
levels(crop2$variable)

df.crop1<-(droplevels(crop2))
df.crop1
df.crop2<-na.omit(subset(df.crop1, variable=="yield"))

###Normality
value.crop<-df.crop2$value
value.crop
#look at boxplot
plot(value.crop~trt,data=df.crop2, main="yield by soil treatment")
plot(df.crop1$value~irrig,data=df.crop1, main="yield by irrigation")

#Check qqplot and histogram
qqnorm(value.crop);qqline(value.crop)
hist(value.crop)# normal
shapiro.test(value.crop) #high p-value .88

#Relabel response variable to use
subset.crop<-df.crop2
crop<-df.crop2$value

#Check residuals, fit linear model
lmer.crop<-lmer (crop ~ irrig*trt*corn + (1|block)+(1|block:irrig)+ (1|block:irrig:corn) + (1|block:irrig:trt),data=subset.crop)
resid.crop <- resid(lmer.crop)
plot(resid.crop);abline(0,0) #Looks good

###C. ANOVA Table - Combine Yield
#Run ANOVA
anova.crop<-anova(lmer.crop, ddf="Kenward-Roger")
anova.crop


###
### Group data plyr package
###

#mutate. mutate works like transform but lets you build on columns you build.
ddply(d, "year", mutate, mu = mean(count), sigma = sd(count),
      + cv = sigma/mu)

ratio.group<-ddply(cig.ee, ~ variable+date+plot+trt) #to sort data


###
### Data wrangling
###

df.e1<-inner_join(df.b1, sample.df.PC.mean,
                  by = c("Matrix", "Enzyme", "pH"))  %>% rename( "act" = y.mean.x, "PC.sorb"=y.mean.y,
                                                                 "se.act" = y.se.x, "se.sorb"=y.se.y );

soil2<- melt (soil1, id.var= c("block", "ssplot", "date", "plot", "trt", "irrig", "corn", "summer", "depth"));
head(soil2)
#droplevels(soil2$variable); str(soil2)#drop additional factor levels

####################################################################################
### Check Residuals....
####################################################################################

#Assess normality of log transformed data = > Check qqplot and histogram

#Insert the transformation of choice here for the response variable
#Choose to use log.value, value or log1
X<-subset.EEAlog2$ratio.log

#Check the residuals
lmer1<- lmer( X ~ irrig*trt*corn + (1|block) + (1|block:irrig) + (1|block:irrig:corn) + (1|block:irrig:trt) , data=EEAratio.log)
resid1 <- resid(lmer1)#plot residuals
par(mfrow=c(1,1));plot(resid1);abline(0,0) 

anova.ee1<-anova(lmer1, ddf="Kenward-Roger")#run anova
anova.ee1

splt.ee1 <- strsplit(as.character(anova.ee1$"Pr(>F)"), ",")
rows.ee1 <- unlist(lapply(splt, function(a) all(as.numeric(a) <= .05))) #insert alpha value here (currently set to .05)
anova.ee1[rows,c(1:6)] # Significant: 

#Check for differences between treatments 
contrasts.trt.ee1<-lsmeans(lmer1, pairwise ~ trt);contrasts.trt.ee1
contrasts.irrig.ee1<-lsmeans(lmer1, pairwise ~ trt);contrasts.irrig.ee1


log.cig1<-cbind(cig4,log(cig4$value)); names(log.cig1)[12]<-c("log.value"); head(log.cig1,2)
log.cig2<-cbind(log.cig1,log(log.cig1$value+1));names(log.cig2)[13]<-c("log1");head(log.cig2,2)
str(log.cig2)
#Change -Inf to 0 and NaNs to 0

head(log.cig2)
log.cig2[,11][log.cig2[,11]<0]=0 #change all negatives to zero, for the untransformed data column 11
log.cig2[,11][log.cig2[,11]== NaN]=0 #change all negatives to zero, for the untransformed data column 11

#long format dataset and remove row 19 (ssplot 33, june - MAJOR OUTLIER)
data.bg1 <-subset(log.cig2, variable=="bg");data.bg<-data.bg1[-19,]

####################################################################################
###
### Run Anova for each enzyme, for each date (the variance between dates is unequal)
###
####################################################################################
#Change out the dataset ag, bg, cb, lap, nag, phos, xyl
subset.ee1<-data.cb
subset.ee.date
subset.ee.date<-subset(subset.ee1, date=="20140727")
value<-subset.ee.date$value
log.value<-subset.ee.date$log.value
log1<-subset.ee.date$log1

###
###Check normality and transform data
###
par(mfrow=c(2,3))
#Assess normality of values, log transformed and log(x+1)data => Check qqplot and histogram
qqnorm(value);qqline(value); qqnorm(log.value);qqline(log.value); qqnorm(log1);qqline(log1)
hist(value); hist(log.value);hist(log1) 
shapiro.value<-shapiro.test(value);shapiro.log<-shapiro.test(log.value);shapiro.log1<-shapiro.test(log1)
shapiro.value
shapiro.log
shapiro.log1

#Change negatives to zeros
data1$no3[data1$no3 < 0]= 0 ; 

#create new factor columns, combine columns, combine colnames, merge
sample.df$MatrixpH<- as.factor(with(sample.df, paste0(Matrix, pH)))

#########################################################################################
###
### Formatting Text
###
#########################################################################################

#Superscript in ggplot
y= (expression("Gravimetric Water Content ("~ g^{-1} ~ "dry soil)" )),

#Subscript and two lines of text for an axis
xlab(expression(paste("line1 \n line2 a" [b])))

#Forming faucet (2 panels) : link here =>http://docs.ggplot2.org/0.9.3.1/facet_grid.html
facet_grid(~ measurement, margins = FALSE, scales = "fixed", space = "fixed", shrink = TRUE, 
           labeller = "label_value", as.table = TRUE, drop = TRUE)

#####################################################################################
###
### Code shared from other people
###
#####################################################################################

#Panel plots from Yamina 
library(lattice)

my.settings <- list(
  superpose.polygon=list(col=c("white", "darkgray"), border="black"),
  strip.background=list(col="white"),
  strip.border=list(col="black")
)

barchart(sum.clean1$mean ~ sum.clean1$irrig | sum.clean1$broad_funcgrp, 
         groups = sum.clean1$amend, 
         data = sum.clean1, 
         auto.key = list(columns = 2, space = "top", title = "Amendment", cex.title = 1.25), 
         layout = c(5,1),
         par.settings = my.settings,
         par.strip.text=list(col="black", font=2, cex = 1.25),
         ylim = c(-0.1,2.4),
         xlab = list(label = "Irrigation", fontsize = 20),
         ylab = list(label = "Biomass (g C" ~ m^{-2} ~ ")", fontsize = 20),
         scales = list(alternating = 1, tck = c(1,0), cex=1.25))

#Here's the ggplot code
library(ggplot2)
err.lims <- aes(ymax = mean + se, ymin = mean - se) 
dodge <- position_dodge(width=0.5)

ggplot(data = sum.clean1, aes(x = irrig, y = mean, fill = amend, alpha = irrig)) +
  geom_bar(stat = "identity", width = 0.5, position = "dodge", colour = "black") +
  geom_errorbar(err.lims, position=dodge, width =0.2) +
  facet_grid(~ broad_funcgrp) +
  scale_alpha_manual(values = c(0.2, 1)) +
  xlab("Trophic Group") +
  ylab(expression("Biomass (g C" ~ m^{-2} ~ ")" )) +
  scale_fill_manual(values = c("B" = "gray", "C" = "white")) +
  theme(axis.title = element_text(size=24), axis.text = element_text(size=15), 
        axis.title.y=element_text(vjust=1), legend.text = element_text(size=18), 
        legend.title=element_blank(), panel.background = element_rect(fill="white"), 
        panel.grid.minor = element_blank())

#Export high-res high resoultion graphs into powerpowert - From Yamina

library(Cairo)
setwd("/Users/yamina/Dropbox/Research/Biochar Fauna Research/FWEB Interpretations/C Fluc Figures/Riverplot")
Cairo(file="CF_rivplot_1000dpi.png", 
      bg="white",
      type="png",
      units="in", 
      width=10, 
      height=8, 
      pointsize=12, 
      dpi=1000)

riverplot(riv, yscale = 0.03, default_style = style, nsteps = 500)

dev.off()

### Friedman Test
wb <- aggregate(warpbreaks$breaks,
                by = list(w = warpbreaks$wool,
                          t = warpbreaks$tension),
                FUN = mean)
wb
friedman.test(wb$x, wb$w, wb$t)
friedman.test(x ~ w | t, data = wb)

########################################
# tranformations
########################################
bc<-boxcox(lm(y.act3 ~ Enzyme*Matrix*pH, data = df.act3way), na.action=na.omit) #run boxcox trans on linear model
str(bc)
lambda <- bc$x[which.max(bc$y)]

powerTransform <- function(y, lambda1, lambda2 = NULL, method = "boxcox") {
  boxcoxTrans <- function(x, lam1, lam2 = NULL) {
    # if we set lambda2 to zero, it becomes the one parameter transformation
    lam2 <- ifelse(is.null(lam2), 0, lam2)
    if (lam1 == 0L) {
      log(y + lam2)
    } else {
      (((y + lam2)^lam1) - 1) / lam1 }}
  switch(method
         , boxcox = boxcoxTrans(y, lambda1, lambda2)
         , tukey = y^lambda1) }

# re-run with transformation
mod.act3.bc <- lm(powerTransform(y.act3, lambda) ~ Enzyme*Matrix*pH, data=df.act3way,
                  na.action = na.omit)
# QQ-plot
qqnorm(mod.act3.bc$residuals); qqline(mod.act3.bc$residuals)

#Tukey HSD adjustment, pull out significant values Tukey significant
tuk.activity3<-as.data.frame(tuk3way$`Enzyme:Matrix:pH`[,c(1,4)])
str(tuk.activity3)
names(tuk.activity3)[1:2]<-c("diff", "pvalue"); str(tuk.activity3)
tuk.act3way<-subset(tuk.activity3, pvalue < .05)
tuk.act3way

########################################
# CREATE A TABLE OF MEANS AND SE FOR PUBLICATION
########################################
(library(doBy))
sumfun <- function(x, na.rm=TRUE){
  c(m=mean(x), se=sd(x)/(sqrt(length(x))))}

str(soil2)
sum.stats1<-summaryBy( yield+mass+HI+wue+nue ~ irrig*fert*char*innoc, data=soil2, FUN = sumfun)
str(sum.stats1)
sum.stats1

#Export p-values from ANOVA for Table 3 - means/p-values
library(xlsx)
write.xlsx(sum.stats1, "C:/Users/Erika/Dropbox/RESEARCH/CIGDrought/DataAnalysis/Statistics/Table1_Raw2.xlsx")

#Create a column of treatment combinations as a dataframe
label3<-as.data.frame(row.names(anova.no3))
names(label3)[1]<- c("treatments") 
label3

table3.no3<-cbind("no3", anova.no3[c(3,6)]);names(table3.no3)[1]<-"var";table3.no3 #include NumDF and Pr(>F)
str(table3.no3)

#Extract p-values and df from anova tables
Table3.pval<-cbind(table3.no3,table3.nh4,table3.mbc, table3.mbn, table3.CN, table3.soilC, table3.soilN, table3.soilP2,
                   table3.gravwc, table3.ph)
head(Table3.pval,2)
#colnames(Table3.pval) #use to get column numbers where p-value are listed (3,5,7,9,11,13,15,17,19)

#Make all p-values >.1 not significant
Table3.pval[,c(3,5,7,9,11,13,15,17,19)][Table3.pval[,c(3,5,7,9,11,13,15,17,19)]>.1]= "ns"
Table3.pval

#This code doesn't work becasue the "ns" transforms those columns into characters, not numeric
#Table3.pval[,c(3,5,7,9,11,13,15,17,19)][Table3.pval[,c(3,5,7,9,11,13,15,17,19)]<.001]= "<.001"

#Export to excel
write.xlsx(Table3.pval, "C:/Users/Erika/Dropbox/RESEARCH/CIGDrought/DataAnalysis/Statistics/Table3pval.xlsx")

contrasts.date.irrig.no3<-lsmeans(lmer.no3, pairwise ~irrig*date);

sum.no3<-lsmeans(lmer.no3, "irrig", by = "date")
plot(sum.no3, by="date")

mean(summary.soilC[1:3,3]) #average BC %C
mean(summary.soilC[4:6,3]) #average C = 1.49
mean(summary.soilC[7:9,3]) #average M = 1.60

###D. Least Squares Means and Contrast Comparisons - soilN
### By irrig:trt:date
contrasts.irrig.trt.date.soilN<-lsmeans(lmer.soilN, pairwise ~ date+irrig+trt); contrasts.irrig.trt.date.soilN
soilN.table<-summary(contrasts.irrig.trt.date.soilN$contrasts)
soilN.table[(soilN.table$p.value<=.05),]

sum.soilN<-lsmeans(lmer.soilN, c("irrig","trt"), by = "date");sum.soilN
par(mfrow=c(1,3))
plot(sum.soilN, by="date")
plot(sum.soilN, by="trt")
plot(sum.soilN, by="irrig")

########################################
# Tables
########################################
Tab.EE<- as.data.frame(cbind(
  c("-Glucosidase", "Acid Phosphatase"),
  c("BG", "PHOS"), 
  c("4-5.5", "4-7"), 
  c(4, 5.2), 
  c("110-112", "240")))
#c(expression(paste("110-112", ^{a})), (expression(paste("240",^{b}))))

colnames(Tab.EE)<-c("Enzyme", "Abbreviation", "Optimal pH", "Isoelectric Point (pI)", "Atomic Weight (kDa)")

Table.EE<-kable(Tab.EE, row.names=FALSE) #%>%
#add_footnote(c("(Watanabe et al, 1992)", "(Durmus et al, 1999)"), notation ="alphabet") 

########################################
# installing Java JDK rJava
########################################
#OSX: On newer versions of OSX you need to install the Java Development Kit.  The normal Java
#runtime environment IS NOT enough.  To get this go to java.com, click "Free java download",
#then IGNORE the big red button, and select "See all java downloads", on the next screen select
#"Looking for the JDK?" from the left hand menu and select the link to "JDK downloads" in the
#first paragraph.  You can then click the "Download" button underneath JDK in the page you are
#taken to.  Sorry this is such a pain!
  
  
#If you're not sure whether you have java installed then you can test this from a command
#prompt.  To get a command prompt try:

#Windows: Select Start > Run, and type 'cmd' (no quotes) in the box which appears, press OK

#MaxOSX: Run Applications > Utilities > Terminal

#Linux: From your applications menu look for an application called 'Terminal' or 'Konsole'.
#Either of these will give you a usable shell.

#At the command prompt type 'java -version' and press enter.  You should see something like:

#java version "1.8.0_60"
#Java(TM) SE Runtime Environment (build 1.8.0_60-b27)
#Java HotSpot(TM) 64-Bit Server VM (build 25.60-b23, mixed mode)

#If you get an error then you don't have java installed.  If the version listed on the first
#line is less than 1.6 then you might have problems running FastQC.

##############################
# Useful websites and links
##############################


#basic anova language http://ww2.coastal.edu/kingw/statistics/R-tutorials/formulae.html

#R cheat sheet:
#https://www.rstudio.com/resources/cheatsheets/

#R Colors :
#http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

#R Color Pop UP: http://www.visibone.com/popups/colhpop.htm

#R PhC code for shapes:
#http://www.endmemo.com/program/R/pchsymbols.php


#Pairwise comparision mulitple graphs
#http://www.statmethods.net/graphs/scatterplot.html

#boxcox transformation
#https://stackoverflow.com/questions/33999512/how-to-use-the-box-cox-power-transformation-in-r/34002020





