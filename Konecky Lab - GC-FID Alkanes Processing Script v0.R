# Alkanes GC-FID Data Import and Analysis
# V0
# Jack A Hutchings, jackh@wustl.edu
# 2018.12.03

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
setwd("C:/Box Sync/Konecky Lab/Data/GCFID Exports/2019_04_16_SNY18 Alkanes Batch 2")

rm(list=ls())

input <- "Konecky Lab - GC-FID Alkanes Processing Template v0.xlsx"
raw_data <- list(gcfid = read_excel(input,sheet=2,na="N/A"),
                 calib = read_excel(input,sheet=4),
                 mix.std = read_excel(input,sheet = 3)[,1:2],
                 parameters = read_excel(input,sheet=3,na="N/A")[,4:5],
                 samples = read_excel(input,sheet=1,na="N/A"))
rm(input)

###Load parameters for use during data processing
parameters <- spread(na.omit(raw_data$parameters),Parameter,Value)

#########  STEP 1  ##############   Calibration Curves   ###################################

calib <- raw_data$calib[,c(1:3)] %>% setNames(.,c("dil","comp","resp")) %>% 
  full_join(raw_data$mix.std,by="comp") %>%  na.omit() %>% mutate(dummy=T) %>% 
  full_join(parameters %>% mutate(dummy=T),by="dummy") %>% select(-dummy) %>% 
  mutate(ug.injected = (conc)* #ug/mL
                     (dil.factor^dil)* #unitless, thus still ug/mL
                     (1/1000)* #unit conversion, now ug/uL
                     (inj.vol) #volume injected, now ug
         ) %>% 
  select(comp,dil,resp,ug.injected) %>% 
  filter(dil > 0)

calib.sum <- group_by(calib,comp) %>% do(fit=lm(ug.injected ~ resp,.)) %>%
  tidy(fit) %>% select(comp,term,estimate) %>% spread(term,estimate) %>% setNames(.,c("comp","intercept","slope")) %>%
  gather(coef,value,intercept,slope) %>% group_by(coef) %>% mutate(mean = mean(value)) %>% spread(comp,value) %>%
  gather(comp,value,-coef) %>% spread(coef,value)

calib.check <- full_join(calib.sum %>% filter(comp!="mean"),calib,by="comp") %>%
  mutate(residual = (slope*resp+intercept) - ug.injected,
         relative.residual = residual/ug.injected) %>% select(dil,comp,intercept,slope,residual,relative.residual) 

calib.stats <-  calib %>% group_by(comp) %>% arrange(dil) %>% do(fit = lm(ug.injected ~ resp,.)) %>% 
  glance(fit)

#Calibration Plot
  ggplot(calib,aes(resp,ug.injected)) +
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
  geom_bar(stat = "identity",position = "dodge") +
  facet_wrap(~comp,scales = "free_y") +
  labs(x = "Dilution #", y = "Relative Residual")

# Fix Me Later!
  
# #Let's make sure your samples lie within the range of you calibration curve
#   raw_data$gcfid[,1:4] %>% setNames(.,c("id","rep","comp","area")) %>% full_join(calib.sum) %>%
#   mutate(ug.actual = slope * area + intercept) %>%
#     full_join(calib) %>%
#     ggplot(aes(log.mg.actual,log.resp,color = id)) +
#     geom_point() +
#     # geom_smooth(method=lm,se=F) +
#     facet_wrap(~comp,scales="free") +
#     theme_bw() +
#     theme(axis.title = element_text(size=18, vjust=0.25),
#           title = element_text(size=16,vjust=0.75),
#           axis.text=element_text(size=20), panel.grid=element_blank(),
#           strip.text = element_text(size=20), legend.title=element_text(size=17), legend.text=element_text(size=17)) +
#     labs(x="Log of Mass Injected - mg",y="Log of Response Area")
                          
#Clean-up from the calibration section
rm(calib.plot,check.calib.bracketing,residual.plot,calib,calib.check,calib.stats)

#######  STEP 2  #################   Data Analysis  ########################################

data <- raw_data$gcfid %>% setNames(.,c("id1","id2","comp","resp")) %>% 
  spread(comp,resp) %>% gather(comp,resp,-c(id1,id2)) %>% mutate(resp = ifelse(is.na(resp),0,resp)) %>% 
  full_join(calib.sum,by="comp") %>%
  filter(comp!="mean") %>% filter(!is.na(slope)) %>% 
  mutate(dummy=T) %>% full_join(parameters %>% mutate(dummy=T),by="dummy") %>% select(-dummy) %>% 
  mutate(ug.injected = resp * slope + intercept,
         ug.injected = ifelse(resp == 0,0,ug.injected),
         irs.theoretical = irs / irs.vol * irs.spike / residue.vol * inj.vol) %>% group_by(id1,id2) %>% 
  mutate(ug.total = ug.injected / inj.vol * residue.vol) %>% # Placeholder ug.total calculation for samples without IRS
  # mutate(recovery = ug.injected[which(comp == "EVAL")] / eval.theoretical,
  #        mg.total = ug.injected / recovery / inj.vol * residue.vol) %>%
  # select(sample,rep,comp,recovery,ug.total) %>% 
  full_join(raw_data$samples,by = c("id1", "id2")) %>%
  mutate(uggs = ug.total / mass_mg * 1000) %>%  #micrograms per gram sediment
  select(id1,id2,carbon_toc,comp,uggs) %>% 
  spread(comp,uggs) %>%
  mutate(
         sigmaLCA = sum(C20,C21,C22,C23,C24,C25,C26,C27,C28,C29,C30,C31,C32,C33,C34,C35),
         lambdaLCA = sigmaLCA*100/carbon_toc,
         ACL = (21*C21+23*C23+25*C25+27*C27+29*C29+31*C31+33*C33+35*C35)/(C21+C23+C25+C27+C29+C31+C33+C35),
         CPI = (C21+C23+C25+C27+C29+C31+C33+C35)/(C20+C22+C24+C26+C28+C30+C32+C34)
         )

# recoveries <- select(data,sample,rep,recovery)

# data.summary <- gather(data,variable,value, c(4:42)) %>% group_by(sample,variable) %>%
#   summarise(mean = mean(value),sd = sd(value),rsd = sd(value)/mean(value))

# data.final <- select(data.summary,sample,variable,mean) %>% spread(variable,mean)

data.final <- data
  
# Output into a csv file! Enjoy... :-)
write_xlsx(data.final,path="output.xlsx")
