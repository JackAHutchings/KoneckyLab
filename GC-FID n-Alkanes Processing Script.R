# Alkanes GC-FID Data Import and Analysis
# Jack A Hutchings, jackh@wustl.edu
# 2019.06.12

# Needed to run this script.
library(dplyr)
library(ggplot2)
library(broom)
library(tidyr)
library(readxl)
library(writexl)

rm(list=ls())
options(scipen=999)

setwd(dirname(rstudioapi::getActiveDocumentContext()$ path))

####################################################################################################

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

calib <- raw_data$calib %>% 
  full_join(raw_data$mix.std,by="comp") %>%  na.omit() %>% mutate(dummy=T) %>% 
  full_join(parameters %>% mutate(dummy=T),by="dummy") %>% select(-dummy) %>% 
  mutate(ng.injected = (conc)* #ng/uL
           (dil.factor^dil)* #unitless, thus still ng/uL
           (inj.vol) #volume injected, now ng
  ) %>% 
  select(batch,comp,dil,area,ng.injected) %>% 
  mutate(log_inj = log(ng.injected),
         log_area = log(area))

calib.sum <- group_by(calib,batch,comp) %>% do(fit=lm(log_inj ~ log_area,.)) %>%
  tidy(fit) %>% select(batch,comp,term,estimate) %>% spread(term,estimate) %>% setNames(.,c("batch","comp","intercept","slope")) %>%
  gather(coef,value,intercept,slope) %>% group_by(batch,coef) %>% mutate(mean = mean(value)) %>% spread(comp,value) %>%
  gather(comp,value,-c(coef,batch)) %>% spread(coef,value)

calib.check <- full_join(calib.sum %>% filter(comp!="mean"),calib,by=c("batch","comp")) %>%
  mutate(residual = exp((slope*log_area+intercept)) - ng.injected,
         relative.residual = residual/ng.injected*100) %>% select(batch,dil,comp,intercept,slope,residual,relative.residual) 


#Calibration Plot
ggplot(calib,aes(log_area,log_inj,color=batch)) +
  geom_point() +
  stat_smooth(method = "lm", se = F,formula = y ~ x) +
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



#######  STEP 2  #################   Data Analysis  ########################################

data <- raw_data$gcfid %>% 
  spread(comp,area) %>% gather(comp,area,-c(id1,id2)) %>% mutate(area = ifelse(is.na(area),0,area)) %>% 
  full_join(raw_data$samples,by = c("id1", "id2")) %>%
  full_join(calib.sum,by=c("batch","comp")) %>%
  filter(comp!="mean") %>% filter(!is.na(slope)) %>% 
  mutate(dummy=T) %>% full_join(parameters %>% mutate(dummy=T),by="dummy") %>% select(-dummy) %>% 
  mutate(ng.injected = exp(log(area) * slope + intercept),
         ng.injected = ifelse(area == 0,0,ng.injected),
         irs.theoretical = irs / irs.vol * irs.spike / volume * inj.vol) %>% group_by(id1,id2) %>% 
  mutate(ug.total = ng.injected / inj.vol * volume / 1000) %>% 
  mutate(uggs = ug.total / mass_mg * 1000) %>%  #micrograms per gram sediment
  select(id1,id2,mass_mg,carbon_toc,comp,uggs) %>% 
  spread(comp,uggs) %>%
  mutate(
    sigmaLCA = sum(C20,C21,C22,C23,C24,C25,C26,C27,C28,C29,C30,C31,C32,C33,C34,C35),
    lambdaLCA = sigmaLCA*100/carbon_toc,
    ACL = (21*C21+23*C23+25*C25+27*C27+29*C29+31*C31+33*C33+35*C35)/(C21+C23+C25+C27+C29+C31+C33+C35),
    CPI = (C21+C23+C25+C27+C29+C31+C33+C35)/(C20+C22+C24+C26+C28+C30+C32+C34),
    C29_ng_total = C29 * mass_mg
  )

readme <- data.frame(readme = c("Units for all individaul compounds and sigmaLCA are:  \u03bcg / g sediment",
                                "simgaLCA is the sum concentration of all saturated alkanes between C20 and C35",
                                "lambdaLCA is the same sum (C20 through C35) reported relative to organic carbon in units of \u03bcg / g OC",
                                "ACL, or average chain length, uses odd-chain alkanes from C21 through C35 and is mass-based (as opposed to molar)",
                                "CPI, or carbon preference index, is the unitless, mass-based ratio of the sum of odd over even alkanes from C20 through C35",
                                "C29_ng_total is the ng total amount of C29 present in the extract. Useful for calculating GC-IRMS volumes."))

# Output into an Excel file! Enjoy... :-)
write_xlsx(list(readme=readme,sample_report=data,calibration_summary=calib.sum,calibration_data=calib),path="output.xlsx")
