# FAMEs GC-FID Data Import and Analysis
# Jack A Hutchings, jackh@wustl.edu
# 2019.06.12

# Please double check the parameters sheet of the processing template to ensure the variables
# match the values used in your extraction procedure.

###################################   Required Files   #####################################

# Needed to run this script.
library(dplyr)
library(ggplot2)
library(broom)
library(tidyr)
library(readxl)
library(writexl)

#You'll need to set your working directory to wherever your data are!
setwd("C:/Box Sync/Konecky Lab/Data/GCFID Exports/2019_05_02_SapTest_FAMEs")

rm(list=ls())

input <- "GC-FID FAMEs Processing Template.xlsx"
raw_data <- list(gcfid = read_excel(input,sheet=2,na="N/A"),
                 calib = read_excel(input,sheet=4),
                 mix.std = read_excel(input,sheet = 3)[,1:2],
                 parameters = read_excel(input,sheet=3,na="N/A")[,4:5],
                 samples = read_excel(input,sheet=1,na="N/A"))
rm(input)

###Load parameters for use during data processing
parameters <- spread(na.omit(raw_data$parameters),Parameter,Value)

#########  STEP 1  ##############   Calibration Curves   ###################################

calib <- raw_data$calib[,c(1:3)] %>% 
  full_join(raw_data$mix.std,by="comp") %>%  na.omit() %>% mutate(dummy=T) %>% 
  full_join(parameters %>% mutate(dummy=T),by="dummy") %>% select(-dummy) %>% 
  mutate(ug.injected = (conc)* #ug/mL
           (dil.factor^dil)* #unitless, thus still ug/mL
           (1/1000)* #unit conversion, now ug/uL
           (inj.vol) #volume injected, now ug
  ) %>% 
  select(comp,dil,area,ug.injected)

calib.sum <- group_by(calib,comp) %>% do(fit=lm(ug.injected ~ area,.)) %>%
  tidy(fit) %>% select(comp,term,estimate) %>% spread(term,estimate) %>% setNames(.,c("comp","intercept","slope")) %>%
  gather(coef,value,intercept,slope) %>% group_by(coef) %>% mutate(mean = mean(value)) %>% spread(comp,value) %>%
  gather(comp,value,-coef) %>% spread(coef,value)

calib.check <- full_join(calib.sum %>% filter(comp!="mean"),calib,by="comp") %>%
  mutate(residual = (slope*area+intercept) - ug.injected,
         relative.residual = residual/ug.injected*100) %>% select(dil,comp,intercept,slope,residual,relative.residual) 

#Calibration Plot
ggplot(calib,aes(area,ug.injected)) +
  geom_point() +
  stat_smooth(method = "lm", se = F,formula = y ~ x,color="blue") +
  facet_wrap(~comp,scales="free") +
  theme_bw() +
  theme(axis.title = element_text(size=18, vjust=0.25),
        title = element_text(size=16,vjust=0.75),
        axis.text=element_text(size=20), panel.grid=element_blank(),
        strip.text = element_text(size=20), legend.title=element_text(size=17), legend.text=element_text(size=17)) +
  labs(x="GC-FID Response (pA)",y="Mass Injected (ug)")

#Residual Plot
calib.check %>% 
  ggplot(aes(dil,relative.residual)) +
  geom_point() +
  facet_wrap(~comp,scales = "free_y") +
  labs(x = "Dilution #", y = "Relative Residual (%)")

#Clean-up from the calibration section
rm(calib,calib.check,calib.stats)

#######  STEP 2  #################   Data Analysis  ########################################

data <- raw_data$gcfid %>% full_join(calib.sum,by="comp") %>%
  full_join(raw_data$samples,by = c("id1", "id2")) %>%
  mutate(intercept = ifelse(is.na(intercept),intercept[which(comp=="mean")],intercept),
         slope = ifelse(is.na(slope),slope[which(comp=="mean")],slope)) %>%
  filter(comp!="mean") %>% 
  mutate(dummy=T) %>% full_join(parameters %>% mutate(dummy=T),by="dummy") %>% select(-dummy) %>% 
  mutate(ug.injected = slope * area + intercept,
         irs.theoretical = irs / irs.vol * irs.spike / volume * inj.vol) %>% group_by(id1,id2) %>% 
  mutate(ug.total = ug.injected / inj.vol * volume) %>% # Placeholder ug.total calculation for samples without IRS
  # mutate(recovery = ug.injected[which(comp == "EVAL")] / eval.theoretical,
  #        mg.total = ug.injected / recovery / inj.vol * residue.vol) %>%
  # select(sample,rep,comp,recovery,ug.total) %>% 
  mutate(uggs = ug.total / mass_mg * 1000, #micrograms per gram sediment
         uggs = ifelse(is.na(uggs),0,uggs)) %>% 
  select(id1,id2,carbon_toc,comp,uggs) %>% 
  spread(comp,uggs) %>% gather(comp,uggs,-c(id1,id2,carbon_toc)) %>% 
  mutate(uggs = ifelse(is.na(uggs),0,uggs)) %>% spread(comp,uggs) %>% 
  mutate(
    sigmaLCFA = sum(C22:C32),
    lambdaLCFA = sigmaLCFA*100/carbon_toc,
    ACL = (22*C22+24*C24+26*C26+28*C28+30*C30+32*C32+34*C34)/(C22+C24+C26+C28+C30+C32+C34),
    CPI = (C22+C24+C26+C28+C30+C32+C34)/(C21+C23+C25+C27+C29+C31+C33)
  )

# Output into an Excel file! Enjoy... :-)
write_xlsx(data,path="output.xlsx")
