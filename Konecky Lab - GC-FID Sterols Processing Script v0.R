# Sterol Data Import and Analysis
# V0.0 {Actually These Are Steroid Alcohols}
# Jack A Hutchings, jahutch2@gmail.com
# 2019.05.02

# Please double check the parameters sheet of the processing template to ensure the variables
# match the values used in your extraction procedure.


###################################   Required Files   #####################################

# Needed to run this script.
library(dplyr)
library(ggplot2)
library(broom)
library(tidyr)
library(readxl)

#You'll need to set your working directory to wherever your data are!
setwd("C:/Box Sync/Konecky Lab/Data/GCFID Exports/2019_04_30_SapTest_NeutralPolars")

rm(list=ls())

input <- "Konecky Lab - GC-FID Sterols Processing Template v0.xlsx"
raw_data <- list(gcfid = read_excel(input,sheet=2,na="N/A"),
                 calib = read_excel(input,sheet=4,na=""),
                 mix.std = read_excel(input,sheet = 3)[,1:3],
                 parameters = read_excel(input,sheet=3,na="N/A")[,4:5],
                 samples = read_excel(input,sheet=1,na="N/A"))
rm(input)

###Load parameters for use during data processing
parameters <- spread(na.omit(raw_data$parameters),Parameter,Value)

#########  STEP 1  ##############   Calibration Curves   ###################################

calib <- raw_data$calib %>%
  full_join(raw_data$mix.std,by="comp") %>% filter(!is.na(comp)) %>%  mutate(dummy=T) %>% 
  full_join(parameters %>% mutate(dummy=T),by="dummy") %>% select(-dummy) %>% 
  group_by(dil,batch) %>% 
  mutate(area = ifelse(is.na(area),0,area),
         rrf = area / area[which(comp=="RFS")]) %>% 
  filter(comp!="RFS") %>% 
  mutate(ug = (mass)* #mg
              (1/mix.vol)* #mg/mL BUT ALSO ug/uL
              (dil.factor^dil)* #ug/uL
              (transfer.vol)* #ug
              (1/(transfer.vol+der.vol))* #ug/uL
              (inj.vol) #ug
         )

# can_curves_be_combined <- calib %>% 
#   select(comp,dil,curve_id,ug,rrf) %>% 
#   group_by(comp) %>% 
#   do(combined = lm(ug~rrf+I(rrf^2)+I(rrf^3),data=.),
#      separate = lm(ug~rrf+I(rrf^2)+I(rrf^3)*curve_id,data=.)) %>% 
#   gather(condition,model,combined,separate) %>% 
#   group_by()
# 
# 
# #Checks the relative residual of each calibration and flags points with a value of >0.25 (i.e., 25%)
# #These points are then removed if and only if they are at the beginning or end of the calibration.
# #Removal of these can reduce the calibration range from the high or low end, but improves the overall quality.
# #This excercise is repeated until no points have a relative residual of >25%. If this process reduces the calibration
# #to less than 5 points, then the entire compound is flagged as bad and should be re-calibrated.
# 
# calib.check.and.expel <- data.frame(relative.residual=1)
# removed.points <- ungroup(calib) %>%  select(curve,dil,comp) %>% filter(!is.numeric(dil))
# while(max(abs(calib.check.and.expel$relative.residual))>0.25){
# calib.check.and.expel <- group_by(calib,comp) %>%
#   anti_join(removed.points,by=c("curve","dil","comp")) %>% 
#   do(fit=lm(ug ~ poly(rrf,3,raw=T),.)) %>%
#   tidy(fit) %>% select(comp,term,estimate) %>% spread(term,estimate) %>% setNames(.,c("comp","int","c","b","a")) %>%
#   full_join(anti_join(calib,removed.points,by=c("curve","dil","comp")),by="comp") %>%
#   group_by(curve,comp) %>%
#   mutate(residual = exp(a * rrf^3 + b * rrf^2 + c * rrf + int) - exp(ug),
#          relative.residual = residual/exp(ug),
#          pass.fail = ifelse((dil==max(dil)|dil==min(dil))&abs(relative.residual)>0.25,"fail","pass")) %>%
#   select(curve,dil,comp,int,a,b,c,residual,relative.residual,pass.fail)
#   removed.points <- full_join(removed.points,filter(calib.check.and.expel,pass.fail=="fail"),by=c("curve","dil","comp")) %>% 
#     select(curve,dil,comp)
# }
# 
# calib.check.and.expel <- group_by(calib.check.and.expel,comp) %>% 
#   mutate(cal.fail = ifelse(length(unique(dil))<5,"fail","pass"))

#Take a look at the removed points to see what was removed.


# calib.sum <- calib.check.and.expel %>% 
#   group_by(comp) %>% 
#   mutate(points = length(unique(dil))) %>% 
#   select(comp:c,points) %>% 
#   distinct() %>% 
#   mutate(bad.cal = ifelse(points<5,"BAD","GOOD"))

calib.sum <- group_by(calib,batch,comp) %>% do(fit=lm(ug ~ rrf,.)) %>% 
  tidy(fit) %>% select(batch,comp,term,estimate) %>% spread(term,estimate) %>% setNames(.,c("batch","comp","intercept","slope")) %>% 
  gather(coef,value,intercept,slope) %>% group_by(batch,coef) %>% mutate(epiCOP = value[which(comp=="COP")],`24eCOP` = value[which(comp=="COP")]) %>% 
  spread(comp,value) %>% 
  gather(comp,value,-c(batch,coef)) %>% spread(coef,value)

#Calibration Plot
  ggplot(calib,aes(rrf,ug)) +
    geom_point() +
    stat_smooth(method = "lm", se = F,formula = y ~ (x),color="blue") +
    facet_wrap(~comp,scales="free") +
    theme_bw() +
    theme(axis.title = element_text(size=18, vjust=0.25),
          title = element_text(size=16,vjust=0.75),
          axis.text=element_text(size=20), panel.grid=element_blank(),
          strip.text = element_text(size=20), legend.title=element_text(size=17), legend.text=element_text(size=17)) +
    labs(x="Response",y="Mass Injected - ug")

# #Residual Plot
  # calib.check %>%
  # ggplot(aes(dil,relative.residual)) +
    # geom_point() +
    # facet_wrap(~comp,scales = "free_y") +
    # labs(x = "Dilution #", y = "Relative Residual")

# #Let's make sure your samples lie within the range of you calibration curve
#   raw_data$gcfid[,1:4] %>% setNames(.,c("id1","id2","comp","resp")) %>%
#     group_by(id1,id2) %>%
#     mutate(rrf = (resp / resp[which(comp=="RFS")])) %>%
#     full_join(calib.sum) %>%
#     mutate(ug = rrf*slope+intercept) %>%
#     full_join(calib) %>%
#     ggplot(aes(ug,rrf,color = id1)) +
#     geom_point() +
#     # geom_smooth(method=lm,se=F) +
#     facet_wrap(~comp,scales="free") +
#     theme_bw() +
#     theme(axis.title = element_text(size=18, vjust=0.25),
#           title = element_text(size=16,vjust=0.75),
#           axis.text=element_text(size=20), panel.grid=element_blank(),
#           strip.text = element_text(size=20), legend.title=element_text(size=17), legend.text=element_text(size=17)) +
#     labs(x="Log of Mass Injected - ug",y="Log of Response Area")
  
  #Table indicating if any of your samples/compounds lie outside the range of you calibration
  # bad.comps <- raw_data$gcms %>% 
  #   group_by(sample,rep) %>% 
  #   mutate(sam.rrf = log(resp / resp[which(comp=="ARS")])) %>% 
  #   filter(comp!="ARS") %>%
  #   select(sample,rep,comp,sam.rrf) %>% 
  #   full_join(calib,by="comp") %>% 
  #   group_by(sample,rep,comp,sam.rrf) %>%
  #   summarise(max.calib.log.resp = max(rrf),
  #             min.calib.log.resp = min(rrf)) %>% 
  #   ungroup() %>%
  #   mutate(sample_status = ifelse(sam.rrf > max.calib.log.resp,"High",
  #                                 ifelse(sam.rrf < min.calib.log.resp,"Low",
  #                                        "Perfect"))) %>%
  #   select(sample,rep,comp,sample_status) %>% 
  #   group_by(sample,rep) %>%
  #   spread(comp,sample_status) %>% 
  #   select(-c(EVAL,MCAD))
                          
#Clean-up from the calibration section
rm(calib.stats,calib,calib.check)

#######  STEP 2  #################   Data Analysis  ########################################

data <- raw_data$gcfid %>% 
  group_by(id1,id2) %>% 
  mutate(rrf = (area / area[which(comp=="RFS")])) %>% 
  filter(comp!="RFS") %>%
  full_join(raw_data$samples,by=c("id1","id2")) %>% 
  full_join(calib.sum,by=c("batch","comp")) %>%
  mutate(dummy=T) %>% full_join(parameters %>% mutate(dummy=T),by="dummy") %>% select(-dummy) %>% 
  mutate(ug.injected = rrf*slope+intercept,
         ug.total.uncorrected = ug.injected / inj.vol * (transfer.vol + der.vol) / transfer.vol * volume,
         irs.expected = with(parameters,irs/irs.vol*irs.spike),
         recovery = ug.total.uncorrected[which(comp=="IRS")]/irs.expected,
         ug.total = ug.total.uncorrected/recovery) %>% 
  filter(comp!="IRS") %>% 
  select(id1,id2,mass_mg,carbon_toc,comp,recovery,ug.total) %>%
  mutate(uggs = ug.total / (mass_mg/1000)) %>% 
  select(id1,id2,comp,recovery,carbon_toc,uggs) %>% 
  spread(comp,uggs) %>% 
  mutate(sigma_sterols = CAMP+CHO+SITO+STIG,
         lamda_sterols = sigma_sterols*100/carbon_toc,
         veg_to_man = COP/`24eCOP`,
         COP_to_epiCOP = COP/epiCOP)

# Output into a csv file! Enjoy... :-)
write.csv(data,file="sterols.csv")

# Bonus Plots :) ----------------------------------------------------------

summary.plot <- lignin.summary %>% filter(variable == "vsc"|variable == "adal.v"|variable == "adal.s"|variable == "c.v"|variable == "s.v") %>%
  ggplot(aes(sample,mean,fill = sample)) +
  geom_col() +
  geom_errorbar(aes(ymax = mean + sd,ymin = mean - sd)) +
  theme_gdocs() + theme(legend.position = "none") +
  facet_wrap(~variable,scales = "free_y")
summary.plot