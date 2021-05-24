

setwd("Z:/Rfun")

data <- read.csv("demos.csv",sep=";",dec=".")
#data <- read.csv2("demos - Copie.csv")

str(data)
data$gender <- as.factor(data$gender)
data$group <- as.factor(data$group)


summary(data)

library(dplyr)
library(tidyverse)
require(MASS)
library(GGally)



############ compute bf for all variables

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

model_intercept <- rlm(cs~ 1,data)
data_Lastmeal <- subset(data, LM != 'NA')  # Remove passive control group
model_intercept_Lastmeal  <- rlm(cs~ 1,data_Lastmeal)
data_stai <- subset(data, stai != 'NA')  # Remove passive control group
model_intercept_stai <- rlm(cs~ 1,data_stai)
data_ratio <- subset(data, Ratio != 'NA')  # Remove passive control group
model_interceptRATIO <- rlm(cs~ 1,data_ratio)

### Compute






bf_BMI = exp((BIC(model_BMI)-BIC(model_intercept))/2)
bf_FWD =exp((BIC(model_FWD)-BIC(model_intercept))/2)
bf_EGGFreq =exp((BIC(model_EGGFreq)-BIC(model_intercept))/2)
bf_EGGPower =exp((BIC(model_EGGPower)-BIC(model_intercept))/2)
bf_TOD =exp((BIC(model_TOD)-BIC(model_intercept))/2)
bf_ratio =exp((BIC(model_ratio)-BIC(model_interceptRATIO))/2)
bf_LFPower =exp((BIC(model_LFPower)-BIC(model_interceptRATIO))/2)
bf_HFPower =exp((BIC(model_Hfpower)-BIC(model_interceptRATIO))/2)
bf_stai =exp((BIC(model_stai)-BIC(model_intercept_stai))/2)
bf_Lastmeal =exp((BIC(model_Lastmeal)-BIC(model_intercept_Lastmeal))/2)
bf_EGGvar =exp((BIC(model_varCycleLength)-BIC(model_intercept))/2)



##### Obtain p and t values for all the variables



summaryBMI<- summary(model_BMI)
t_BMI <- summaryBMI$coefficients [2, "t value"]
p_BMI <- 2 * pt(abs(t_BMI), df = nobs(model_BMI) - 2, lower.tail = FALSE)

summarySTAI<- summary(model_stai)
t_STAI <- summarySTAI$coefficients [2, "t value"]
p_STAI <- 2 * pt(abs(t_STAI), df = nobs(model_stai) - 2, lower.tail = FALSE)

summaryTOD<- summary(model_TOD)
t_TOD <- summaryTOD$coefficients [2, "t value"]
p_TOD <- 2 * pt(abs(t_TOD), df = nobs(model_TOD) - 2, lower.tail = FALSE)

summaryLM<- summary(model_Lastmeal)
t_LM <- summaryLM$coefficients [2, "t value"]
p_LM <- 2 * pt(abs(t_LM), df = nobs(model_Lastmeal) - 2, lower.tail = FALSE)

summaryFWD<- summary(model_FWD)
t_FWD <- summaryFWD$coefficients [2, "t value"]
p_FWD <- 2 * pt(abs(t_FWD), df = nobs(model_FWD) - 2, lower.tail = FALSE)


summaryEGGF<- summary(model_EGGFreq)
t_EGGF <- summaryEGGF$coefficients [2, "t value"]
p_EGGF <- 2 * pt(abs(t_EGGF), df = nobs(model_EGGFreq) - 2, lower.tail = FALSE)

summaryEGGP<- summary(model_EGGPower)
t_EGGP <- summaryEGGP$coefficients [2, "t value"]
p_EGGP <- 2 * pt(abs(t_EGGP), df = nobs(model_EGGPower) - 2, lower.tail = FALSE)

summaryEGGV<- summary(model_varCycleLength)
t_EGGV <- summaryEGGV$coefficients [2, "t value"]
p_EGGV <- 2 * pt(abs(t_EGGV), df = nobs(model_varCycleLength) - 2, lower.tail = FALSE)

summaryLF<- summary(model_LFPower)
t_LF <- summaryLF$coefficients [2, "t value"]
p_LF <- 2 * pt(abs(t_LF), df = nobs(model_LFPower) - 2, lower.tail = FALSE)

summaryHF<- summary(model_Hfpower)
t_HF <- summaryHF$coefficients [2, "t value"]
p_HF <- 2 * pt(abs(t_HF), df = nobs(model_Hfpower) - 2, lower.tail = FALSE)

summaryRatio<- summary(model_ratio)
t_Ratio <- summaryRatio$coefficients [2, "t value"]
p_Ratio <- 2 * pt(abs(t_Ratio), df = nobs(model_ratio) - 2, lower.tail = FALSE)

all_p <- c(p_BMI,p_STAI,p_TOD,p_LM,p_FWD,p_EGGF,p_EGGP,p_EGGV,p_LF,p_HF,p_Ratio)
all_t <- c(t_BMI,t_STAI,t_TOD,t_LM,t_FWD,t_EGGF,t_EGGP,t_EGGV,t_LF,t_HF,t_Ratio)
all_bf<- c(bf_BMI,bf_stai,bf_TOD,bf_Lastmeal,bf_FWD,bf_EGGFreq,bf_EGGPower,bf_EGGvar,bf_LFPower,bf_HFPower,bf_ratio)

d_min <-c(min(data$BMI,na.rm = T),min(data$stai,na.rm = T),min(data$TOD,na.rm = T),min(data$LM,na.rm = T),min(data$FWD,na.rm = T),min(data$EGGFreq,na.rm = T),min(data$EGGPower,na.rm = T),min(data$varCycleLength,na.rm = T),min(data$LFPower,na.rm = T),min(data$Hfpower,na.rm = T),min(data$Ratio,na.rm = T))
d_max <-c(max(data$BMI,na.rm = T),max(data$stai,na.rm = T),max(data$TOD,na.rm = T),max(data$LM,na.rm = T),max(data$FWD,na.rm = T),max(data$EGGFreq,na.rm = T),max(data$EGGPower,na.rm = T),max(data$varCycleLength,na.rm = T),max(data$LFPower,na.rm = T),max(data$Hfpower,na.rm = T),max(data$Ratio,na.rm = T))
d_mean <-c(mean(data$BMI,na.rm = T),mean(data$stai,na.rm = T),mean(data$TOD,na.rm = T),mean(data$LM,na.rm = T),mean(data$FWD,na.rm = T),mean(data$EGGFreq,na.rm = T),mean(data$EGGPower,na.rm = T),mean(data$varCycleLength,na.rm = T),mean(data$LFPower,na.rm = T),mean(data$Hfpower,na.rm = T),mean(data$Ratio,na.rm = T))
d_median <-c(median(data$BMI,na.rm = T),median(data$stai,na.rm = T),median(data$TOD,na.rm = T),median(data$LM,na.rm = T),median(data$FWD,na.rm = T),median(data$EGGFreq,na.rm = T),median(data$EGGPower,na.rm = T),median(data$varCycleLength,na.rm = T),median(data$LFPower,na.rm = T),median(data$Hfpower,na.rm = T),median(data$Ratio,na.rm = T))
d_mean <-c(mean(data$BMI,na.rm = T),mean(data$stai,na.rm = T),mean(data$TOD,na.rm = T),mean(data$LM,na.rm = T),mean(data$FWD,na.rm = T),mean(data$EGGFreq,na.rm = T),mean(data$EGGPower,na.rm = T),mean(data$varCycleLength,na.rm = T),mean(data$LFPower,na.rm = T),mean(data$Hfpower,na.rm = T),mean(data$Ratio,na.rm = T))
d_sd <-c(sd(data$BMI,na.rm = T),sd(data$stai,na.rm = T),sd(data$TOD,na.rm = T),sd(data$LM,na.rm = T),sd(data$FWD,na.rm = T),sd(data$EGGFreq,na.rm = T),sd(data$EGGPower,na.rm = T),sd(data$varCycleLength,na.rm = T),sd(data$LFPower,na.rm = T),sd(data$Hfpower,na.rm = T),sd(data$Ratio,na.rm = T))



ptbf_stats<- data.frame(all_t,all_p,all_bf,d_min,d_max,d_mean,d_median,d_sd)




library(clipr)
clipr::write_clip(ptbf_stats)

#### plot candidate models EGG freq
library(ggplot2)

model1 <-rlm(EGGFreq~LFPower,data)
model1_summary <- summary(model1)
t_1_LFEGG <- model1_summary$coefficients [2, "t value"]
p_1_LFEGG <- 2 * pt(abs(t_1_LFEGG), df = nobs(model1) - 2, lower.tail = FALSE)

dev.new()
ggplot(data, aes(x=LFPower, y=jitter(EGGFreq,4))) + 
  geom_point(shape=20,size=3, color="black")+
  geom_smooth(method=rlm,se =FALSE ,color="black",fullrange=TRUE)+
  labs(title=paste("Robust regression EGG Frequency ~ LF HRV power p ",p_1_LFEGG),
       y="EGG frequency (Hz.)", x = "Low frequency power *10e-4 ms2")

model1 <-rlm(EGGFreq~Hfpower,data)
model1_summary <- summary(model1)
t <- model1_summary$coefficients [2, "t value"]
p <- 2 * pt(abs(t), df = nobs(model1) - 2, lower.tail = FALSE)

dev.new()
ggplot(data, aes(x=Hfpower, y=jitter(EGGFreq,4))) + 
  geom_point(shape=20,size=3, color="black")+
  geom_smooth(method=rlm,se =FALSE ,color="black",fullrange=TRUE)+
  labs(title=paste("Robust regression EGG Frequency ~ HF HRV power p ",p),
       y="EGG frequency (Hz.)", x = "High frequency power *10e-4 ms2")

model1 <-rlm(EGGFreq~Ratio,data)
model1_summary <- summary(model1)
t <- model1_summary$coefficients [2, "t value"]
p <- 2 * pt(abs(t), df = nobs(model1) - 2, lower.tail = FALSE)

dev.new()
ggplot(data, aes(x=Ratio, y=jitter(EGGFreq,4))) + 
  geom_point(shape=20,size=3, color="black")+
  geom_smooth(method=rlm,se =FALSE ,color="black",fullrange=TRUE)+
  labs(title=paste("Robust regression EGG Frequency ~ Ratio HRV power p ",p),
       y="EGG frequency (Hz.)", x = "Ratio LF/HF")


### cs

model1 <-rlm(cs~LFPower,data)
model1_summary <- summary(model1)
t_1_LFEGG <- model1_summary$coefficients [2, "t value"]
p_1_LFEGG <- 2 * pt(abs(t_1_LFEGG), df = nobs(model1) - 2, lower.tail = FALSE)

dev.new()
ggplot(data, aes(x=LFPower, y=jitter(cs,4))) + 
  geom_point(shape=20,size=3, color="black")+
  geom_smooth(method=rlm,se =FALSE ,color="black",fullrange=TRUE)+
  labs(title=paste("Robust regression coupling strenght ~ LF HRV power p ",p_1_LFEGG),
       y="coupling strenght ", x = "Low frequency power *10e-4 ms2")

model1 <-rlm(cs~Hfpower,data)
model1_summary <- summary(model1)
t <- model1_summary$coefficients [2, "t value"]
p <- 2 * pt(abs(t), df = nobs(model1) - 2, lower.tail = FALSE)

dev.new()
ggplot(data, aes(x=Hfpower, y=jitter(cs,4))) + 
  geom_point(shape=20,size=3, color="black")+
  geom_smooth(method=rlm,se =FALSE ,color="black",fullrange=TRUE)+
  labs(title=paste("Robust regression coupling strenght ~ HF HRV power p ",p),
       y="coupling strenght ", x = "High frequency power *10e-4 ms2")

model1 <-rlm(cs~Ratio,data)
model1_summary <- summary(model1)
t <- model1_summary$coefficients [2, "t value"]
p <- 2 * pt(abs(t), df = nobs(model1) - 2, lower.tail = FALSE)

dev.new()
ggplot(data, aes(x=Ratio, y=jitter(cs,4))) + 
  geom_point(shape=20,size=3, color="black")+
  geom_smooth(method=rlm,se =FALSE ,color="black",fullrange=TRUE)+
  labs(title=paste("Robust regression coupling strenght ~ Ratio HRV power p ",p),
       y="coupling strenght ", x = "Ratio LF/HF")


model1 <-rlm(boldSTDgasnet~Ratio,data)
model1_summary <- summary(model1)
t <- model1_summary$coefficients [2, "t value"]
p <- 2 * pt(abs(t), df = nobs(model1) - 2, lower.tail = FALSE)

dev.new()
ggplot(data, aes(x=Ratio, y=jitter(boldSTDgasnet,4))) + 
  geom_point(shape=20,size=3, color="black")+
  geom_smooth(method=rlm,se =FALSE ,color="black",fullrange=TRUE)+
  labs(title=paste("Robust regression BOLDvariance gasnet ~ Ratio HRV power p ",p),
       y="BOLD std in gasnet ", x = "Ratio LF/HF")


### EGG power


model1 <-rlm(EGGPower~LFPower,data)
model1_summary <- summary(model1)
t_1_LFEGG <- model1_summary$coefficients [2, "t value"]
p_1_LFEGG <- 2 * pt(abs(t_1_LFEGG), df = nobs(model1) - 2, lower.tail = FALSE)

dev.new()
ggplot(data, aes(x=LFPower, y=jitter(EGGPower,4))) + 
  geom_point(shape=20,size=3, color="black")+
  geom_smooth(method=rlm,se =FALSE ,color="black",fullrange=TRUE)+
  labs(title=paste("Robust regression EGGPower ~ LF HRV power p ",p_1_LFEGG),
       y="EGGPower ", x = "Low frequency power *10e-4 ms2")

model1 <-rlm(EGGPower~Hfpower,data)
model1_summary <- summary(model1)
t <- model1_summary$coefficients [2, "t value"]
p <- 2 * pt(abs(t), df = nobs(model1) - 2, lower.tail = FALSE)

dev.new()
ggplot(data, aes(x=Hfpower, y=jitter(EGGPower,4))) + 
  geom_point(shape=20,size=3, color="black")+
  geom_smooth(method=rlm,se =FALSE ,color="black",fullrange=TRUE)+
  labs(title=paste("Robust regression EGGPower ~ HF HRV power p ",p),
       y="EGGPower ", x = "High frequency power *10e-4 ms2")

model1 <-rlm(EGGPower~Ratio,data)
model1_summary <- summary(model1)
t <- model1_summary$coefficients [2, "t value"]
p <- 2 * pt(abs(t), df = nobs(model1) - 2, lower.tail = FALSE)

dev.new()
ggplot(data, aes(x=Ratio, y=jitter(EGGPower,4))) + 
  geom_point(shape=20,size=3, color="black")+
  geom_smooth(method=rlm,se =FALSE ,color="black",fullrange=TRUE)+
  labs(title=paste("Robust regression EGGPower ~ Ratio HRV power p ",p),
       y="EGGPower ", x = "Ratio LF/HF")



########
# Compare models EGGfreq



model2 <-rlm(EGGFreq~ Ratio + Hfpower + LFPower,data)
model3 <-rlm(EGGFreq~  Hfpower,data)
model4 <-rlm(EGGFreq~  Hfpower + Ratio,data)
model5 <-rlm(EGGFreq~  LFPower + Ratio,data)
model6 <-rlm(EGGFreq~ Ratio,data)

BIC(model2)
BIC(model3)
BIC(model4)
BIC(model5)
BIC(model6)

# Compare models cs

model1 <-rlm(cs~ Ratio,data)
model2 <-rlm(cs~ Ratio + Hfpower + LFPower,data)
model3 <-rlm(cs~ Ratio + EGGFreq,data)
model4 <-rlm(cs~  Hfpower + Ratio,data)
model5 <-rlm(cs~  LFPower + Ratio,data)
model6 <-rlm(cs ~ Ratio +EGGFreq + Hfpower,data)


BIC(model1)
BIC(model2)
BIC(model3)
BIC(model4)
BIC(model5)
BIC(model6)


# bayes factor EGG frequency and 


