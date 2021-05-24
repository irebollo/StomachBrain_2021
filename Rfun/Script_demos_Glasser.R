
setwd("Z:/Rfun")

#data <- read.csv("CS_glasserROIs_LR.csv",sep=";",dec=",",header=TRUE)
data <- read.csv("CS_glasserROIs_LR2.csv",sep="\t",dec=".",header=TRUE)

#data <- read.cvs2("demos - Copie.csv")

#str(data)
data$gender <- as.factor(data$gender)
data$group <- as.factor(data$group)
data$subject<- as.factor(data$subject)


summary(data)

library(dplyr)
library(tidyverse)
require(MASS)
library(MASS)
library(GGally)

library(tidyr)
data_long <- gather(data, condition, measurement, A1:A360, factor_key=TRUE)

summary(data_long)




regressions<-do.call("rbind",by(data_long, data_long[,c("condition")], function(x) rlm(formula=measurement~Ratio, data=x)$coefficients ))
labels<-unique(data_long[,c("condition")])
fits<-cbind(labels,as.data.frame(regressions))




dataNonZero<-subset(data_long2, measurement!=0)
labelss<-unique(dataNonZero[,c("condition")])
labelss=as.character(labelss)

areas = labelss
areasNames = sprintf("A%s",areas)





model_BMI  <-rlm(cs~BMI,data)
model_FWD  <-rlm(cs~FWD,data)
model_EGGFreq  <-rlm(cs~EGGFreq,data)
model_EGGPower  <-rlm(cs~EGGPower,data)
model_TOD  <-rlm(cs~TOD,data)
model_stai  <-rlm(cs~stai,data)
model_Lastmeal  <-rlm(cs~LM,data)
model_LFPower  <-rlm(cs~LFPower,data)
model_Hfpower  <-rlm(cs~Hfpower,data)
model_ratio <- rlm(cs~ Ratio,data_ratio)
model_varCycleLength  <-rlm(cs~varCycleLength,data)




BF_BMI <- vector(mode = "list", length = length(labelss))
BF_FWD <- vector(mode = "list", length = length(labelss))
BF_EGGfreq <- vector(mode = "list", length = length(labelss))
BF_EGGPower <- vector(mode = "list", length = length(labelss))
BF_TOD <- vector(mode = "list", length = length(labelss))
BF_stai <- vector(mode = "list", length = length(labelss))
BF_Lastmeal <- vector(mode = "list", length = length(labelss))
BF_LFPower <- vector(mode = "list", length = length(labelss))
BF_Hfpower <- vector(mode = "list", length = length(labelss))
BF_ratio <- vector(mode = "list", length = length(labelss))
BF_varCycleLength <- vector(mode = "list", length = length(labelss))


P_BMI <- vector(mode = "list", length = length(labelss))
P_FWD <- vector(mode = "list", length = length(labelss))
P_EGGfreq <- vector(mode = "list", length = length(labelss))
P_EGGPower <- vector(mode = "list", length = length(labelss))
P_TOD <- vector(mode = "list", length = length(labelss))
P_stai <- vector(mode = "list", length = length(labelss))
P_Lastmeal <- vector(mode = "list", length = length(labelss))
P_LFPower <- vector(mode = "list", length = length(labelss))
P_Hfpower <- vector(mode = "list", length = length(labelss))
P_ratio <- vector(mode = "list", length = length(labelss))
P_varCycleLength <- vector(mode = "list", length = length(labelss))


T_BMI <- vector(mode = "list", length = length(labelss))
T_FWD <- vector(mode = "list", length = length(labelss))
T_EGGfreq <- vector(mode = "list", length = length(labelss))
T_EGGPower <- vector(mode = "list", length = length(labelss))
T_TOD <- vector(mode = "list", length = length(labelss))
T_stai <- vector(mode = "list", length = length(labelss))
T_Lastmeal <- vector(mode = "list", length = length(labelss))
T_LFPower <- vector(mode = "list", length = length(labelss))
T_Hfpower <- vector(mode = "list", length = length(labelss))
T_ratio <- vector(mode = "list", length = length(labelss))
T_varCycleLength <- vector(mode = "list", length = length(labelss))

data_ratio <- subset(dataNonZero, Ratio != 'NA')  # Remove passive control group
data_Lastmeal <- subset(dataNonZero, LM != 'NA')  # Remove passive control group
data_stai <- subset(dataNonZero, stai != 'NA')  # Remove passive control group



for ( i in 1:length(labelss)){
  currentArea = areasNames[i]
  
  

  
  
  
  
  # BMI
  sss<-subset(dataNonZero, condition==currentArea)
  model0= rlm(measurement~1, sss)
  model1= rlm(measurement~BMI, sss)
  BICC0 = BIC(model0)
  BICC1 = BIC(model1)
  BF_BMI[i]= exp((BICC1-BICC0)/2)
  model1_summary <- summary(model1)
  t <- model1_summary$coefficients [2, "t value"]
  p <- 2 * pt(abs(t), df = nobs(model1) - 2, lower.tail = FALSE)
  T_BMI[i]= t
  P_BMI[i]= p
  
  
  # FWD
  sss<-subset(dataNonZero, condition==currentArea)
  model0= rlm(measurement~1, sss)
  model1= rlm(measurement~FWD, sss)
  BICC0 = BIC(model0)
  BICC1 = BIC(model1)
  BF_FWD[i]= exp((BICC1-BICC0)/2)
  model1_summary <- summary(model1)
  t <- model1_summary$coefficients [2, "t value"]
  p <- 2 * pt(abs(t), df = nobs(model1) - 2, lower.tail = FALSE)
  T_FWD[i]= t
  P_FWD[i]= p
  
  
  
  
  
  # EGGfreq
  sss<-subset(dataNonZero, condition==currentArea)
  model0= rlm(measurement~1, sss)
  model1= rlm(measurement~EGGFreq, sss)
  BICC0 = BIC(model0)
  BICC1 = BIC(model1)
  BF_EGGfreq[i]= exp((BICC1-BICC0)/2)
  model1_summary <- summary(model1)
  t <- model1_summary$coefficients [2, "t value"]
  p <- 2 * pt(abs(t), df = nobs(model1) - 2, lower.tail = FALSE)
  T_EGGfreq[i]= t
  P_EGGfreq[i]= p
  
  # EGGPower
  sss<-subset(dataNonZero, condition==currentArea)
  model0= rlm(measurement~1, sss)
  model1= rlm(measurement~EGGPower, sss)
  BICC0 = BIC(model0)
  BICC1 = BIC(model1)
  BF_EGGPower[i]= exp((BICC1-BICC0)/2)
  model1_summary <- summary(model1)
  t <- model1_summary$coefficients [2, "t value"]
  p <- 2 * pt(abs(t), df = nobs(model1) - 2, lower.tail = FALSE)
  
  T_EGGPower[i]= t
  P_EGGPower[i]= p
  
  # TOD
  sss<-subset(dataNonZero, condition==currentArea)
  model0= rlm(measurement~1, sss)
  model1= rlm(measurement~TOD, sss)
  BICC0 = BIC(model0)
  BICC1 = BIC(model1)
  BF_TOD[i]= exp((BICC1-BICC0)/2)
  model1_summary <- summary(model1)
  t <- model1_summary$coefficients [2, "t value"]
  p <- 2 * pt(abs(t), df = nobs(model1) - 2, lower.tail = FALSE)
  
  T_TOD[i]= t
  P_TOD[i]= p
  
  # stai
  sss<-subset(data_stai, condition==currentArea)
  model0= rlm(measurement~1, sss)
  model1= rlm(measurement~stai, sss)
  BICC0 = BIC(model0)
  BICC1 = BIC(model1)
  BF_stai[i]= exp((BICC1-BICC0)/2)
  model1_summary <- summary(model1)
  t <- model1_summary$coefficients [2, "t value"]
  p <- 2 * pt(abs(t), df = nobs(model1) - 2, lower.tail = FALSE)
  
  T_stai[i]= t
  P_stai[i]= p
  
  # LM
  sss<-subset(data_Lastmeal, condition==currentArea)
  model0= rlm(measurement~1, sss)
  model1= rlm(measurement~LM, sss)
  BICC0 = BIC(model0)
  BICC1 = BIC(model1)
  BF_Lastmeal[i]= exp((BICC1-BICC0)/2)
  model1_summary <- summary(model1)
  t <- model1_summary$coefficients [2, "t value"]
  p <- 2 * pt(abs(t), df = nobs(model1) - 2, lower.tail = FALSE)
  
  T_Lastmeal[i]= t
  P_Lastmeal[i]= p
  
  # Hfpower
  sss<-subset(data_ratio, condition==currentArea)
  model0= rlm(measurement~1, sss)
  model1= rlm(measurement~Hfpower, sss)
  BICC0 = BIC(model0)
  BICC1 = BIC(model1)
  BF_Hfpower[i]= exp((BICC1-BICC0)/2)
  
  t <- model1_summary$coefficients [2, "t value"]
  p <- 2 * pt(abs(t), df = nobs(model1) - 2, lower.tail = FALSE)
  
  
  T_Hfpower[i]= t
  P_Hfpower[i]= p
  
  # LFPower
  sss<-subset(data_ratio, condition==currentArea)
  model0= rlm(measurement~1, sss)
  model1= rlm(measurement~LFPower, sss)
  BICC0 = BIC(model0)
  BICC1 = BIC(model1)
  BF_LFPower[i]= exp((BICC1-BICC0)/2)
  model1_summary <- summary(model1)
  t <- model1_summary$coefficients [2, "t value"]
  p <- 2 * pt(abs(t), df = nobs(model1) - 2, lower.tail = FALSE)
  
  T_LFPower[i]= t
  P_LFPower[i]= p
  
  
  # ratio
  sss<-subset(data_ratio, condition==currentArea)
  model0= rlm(measurement~1, sss)
  model1= rlm(measurement~Ratio, sss)
  BICC0 = BIC(model0)
  BICC1 = BIC(model1)
  BF_ratio[i]= exp((BICC1-BICC0)/2)
  model1_summary <- summary(model1)
  t <- model1_summary$coefficients [2, "t value"]
  p <- 2 * pt(abs(t), df = nobs(model1) - 2, lower.tail = FALSE)
  
  T_ratio[i]= t
  P_ratio[i]= p
  

  # varCycleLength
  sss<-subset(dataNonZero, condition==currentArea)
  model0= rlm(measurement~1, sss)
  model1= rlm(measurement~varCycleLength, sss)
  BICC0 = BIC(model0)
  BICC1 = BIC(model1)
  BF_varCycleLength[i]= exp((BICC1-BICC0)/2)
  model1_summary <- summary(model1)
  t <- model1_summary$coefficients [2, "t value"]
  p <- 2 * pt(abs(t), df = nobs(model1) - 2, lower.tail = FALSE)
  T_varCycleLength[i]= t
  P_varCycleLength[i]= p
}



adf <- function(x) {
xxx<-as.data.frame(t(x))
xx <-as.data.frame(t(xxx))
return (xx)
}




BFS <-cbind(areasNames,adf(BF_BMI),adf(T_BMI),adf(P_BMI),
            adf(BF_FWD),adf(T_FWD),adf(P_FWD),
            adf(BF_EGGfreq),adf(T_EGGfreq),adf(P_EGGfreq),
            adf(BF_EGGPower),adf(T_EGGPower),adf(P_EGGPower),
            adf(BF_TOD),adf(T_TOD),adf(P_TOD),
            adf(BF_stai),adf(T_stai),adf(P_stai),
            adf(BF_Lastmeal),adf(T_Lastmeal),adf(P_Lastmeal),
            adf(BF_LFPower),adf(T_LFPower),adf(P_LFPower),
            adf(BF_Hfpower),adf(T_Hfpower),adf(P_Hfpower),
            adf(BF_ratio),adf(T_ratio),adf(P_ratio),
            adf(BF_varCycleLength),adf(T_varCycleLength),adf(P_varCycleLength))

BFS = as.data.frame(BFS)

library(clipr)
clipr::write_clip(BFS,object_type="table")

write.csv2(as.data.frame(BFS),'BFglasser.csv')

install.packages('data.table')
library(data.table)

BFS <-cbind(adf(BF_BMI),adf(T_BMI),adf(P_BMI),
            adf(BF_FWD),adf(T_FWD),adf(P_FWD),
            adf(BF_EGGfreq),adf(T_EGGfreq),adf(P_EGGfreq),
            adf(BF_EGGPower),adf(T_EGGPower),adf(P_EGGPower),
            adf(BF_TOD),adf(T_TOD),adf(P_TOD),
            adf(BF_stai),adf(T_stai),adf(P_stai),
            adf(BF_Lastmeal),adf(T_Lastmeal),adf(P_Lastmeal),
            adf(BF_LFPower),adf(T_LFPower),adf(P_LFPower),
            adf(BF_Hfpower),adf(T_Hfpower),adf(P_Hfpower),
            adf(BF_ratio),adf(T_ratio),adf(P_ratio),
            adf(BF_varCycleLength),adf(T_varCycleLength),adf(P_varCycleLength))
fwrite(BFS, file ="BFglasserNONAMES.csv")




model_FWD<-do.call("rbind",by(data_long, data_long[,c("condition")], function(x) rlm(formula=measurement~FWD, data=x)$coefficients ))
labels<-unique(data_long[,c("condition")])
fits<-cbind(labels,as.data.frame(model_BMI))


model_FWD<-do.call("rbind",by(data_long, data_long[,c("condition")], function(x) rlm(formula=measurement~FWD, data=x)$summary ))
labels<-unique(data_long[,c("condition")])
fits<-cbind(labels,as.data.frame(model_BMI))



long_fits = reshape2::melt(as.data.frame(fits),id.vars = c("labels"))