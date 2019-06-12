# Alkanes GC-FID Data Import and Analysis
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
#This is easiest done using the Session > Set Working Directory option in the RStudio GUI
setwd()

rm(list=ls())

input <- "GC-FID n-Alkanes Processing Template.xlsx"
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
  mutate(ng.injected = (conc)* #ng/uL
           (dil.factor^dil)* #unitless, thus still ng/uL
           (inj.vol) #volume injected, now ng
  ) %>% 
  select(comp,dil,area,ng.injected) %>% 
  filter(dil > 0)

calib.sum <- group_by(calib,comp) %>% do(fit=lm(ng.injected ~ area,.)) %>%
  tidy(fit) %>% select(comp,term,estimate) %>% spread(term,estimate) %>% setNames(.,c("comp","intercept","slope")) %>%
  gather(coef,value,intercept,slope) %>% group_by(coef) %>% mutate(mean = mean(value)) %>% spread(comp,value) %>%
  gather(comp,value,-coef) %>% spread(coef,value)

calib.check <- full_join(calib.sum %>% filter(comp!="mean"),calib,by="comp") %>%
  mutate(residual = (slope*area+intercept) - ng.injected,
         relative.residual = residual/ng.injected*100) %>% select(dil,comp,intercept,slope,residual,relative.residual) 


#Calibration Plot
ggplot(calib,aes(area,ng.injected)) +
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

#Clean-up from the calibration section
rm(calib,calib.check)

#######  STEP 2  #################   Data Analysis  ########################################

data <- raw_data$gcfid %>% 
  spread(comp,area) %>% gather(comp,area,-c(id1,id2)) %>% mutate(area = ifelse(is.na(area),0,area)) %>% 
  full_join(calib.sum,by="comp") %>%
  filter(comp!="mean") %>% filter(!is.na(slope)) %>% 
  mutate(dummy=T) %>% full_join(parameters %>% mutate(dummy=T),by="dummy") %>% select(-dummy) %>% 
  mutate(ng.injected = area * slope + intercept,
         ng.injected = ifelse(area == 0,0,ng.injected),
         irs.theoretical = irs / irs.vol * irs.spike / residue.vol * inj.vol) %>% group_by(id1,id2) %>% 
  mutate(ug.total = ng.injected / inj.vol * residue.vol / 1000) %>% # Placeholder ug.total calculation for samples without IRS
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


# Output into an Excel file! Enjoy... :-)
write_xlsx(data,path="output.xlsx")
