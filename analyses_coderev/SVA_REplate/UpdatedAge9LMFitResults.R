library(ewastools)
library(limma)
library(readr)
library(tidyverse)

#Load in the SER+random effect results (lmfitblock.R)
age9_RE_full <- readRDS(file="/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERfullblock.rds")
age9_RE_sig <- readRDS(file="/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERsig247block.rds")

#Load in the SER + SV surrogate variables

age9SER_SVfull <- readRDS("/home/alreiner/Projects/ffcw/output/age9SERSVfull.rds")
age9SER_SVsig <- readRDS("/home/alreiner/Projects/ffcw/output/age9SERSVsig247.rds")

overlap_nums <- function(age9, age15){
  union <- length(union(age9, age15)) # elements in union (all probes found sig between the two ages)
  #intersect=returns the elements that are in both sets (the sig probes found in both age9 and age15)
  intersect(age9, age15) # 
  #make list of overlapping probes
  overlap <- length(intersect(age9, age15))
  inter <- overlap/union
  return(inter) #concordance of ALL sig probes
}

overlap_nums(rownames(age9_RE_sig),rownames(age9SER_SVsig)) #83.3% overlap


comp1 <- merge(age9SER_SVfull %>% select(probe_id,P.Value,logFC),
               age9_RE_full%>%select(probe_id,P.Value,logFC), by="probe_id", all=FALSE)
comp1 <- comp1 %>% rename("SV_PVal"="P.Value.x",
                          "SV_logFC"="logFC.x",
                        "RE_PVal"="P.Value.y",
                        "RE_logFC"="logFC.y")

#calculate correlation between p values
hist(comp1$SV_PVal)
cor(comp1$SV_PVal, comp1$RE_PVal, method="pearson") #0.7543
cor(comp1$SV_PVal, comp1$RE_PVal, method="spearman")#0.8116 #reporting this because variablels
#are not normally ddistributed

#calculate correlation between effect estimates for sex
hist(comp1$SV_logFC,breaks=100)
hist(comp1$RE_logFC)
cor(comp1$SV_logFC, comp1$RE_logFC, method="pearson") #0.9867
cor(comp1$SV_logFC, comp1$RE_logFC, method="spearman")#0.9307

sum(comp1$SV_logFC > 0.1)
sum(comp1$SV_logFC < -0.1)


compsig <- ggplot(data=comp1, aes(x=-log10(SV_PVal),y=-log10(RE_PVal)))+
  #geom_jitter(size=.5,width = 0.1)+
  geom_point(size=.5)+
  theme_classic()+
  theme(axis.title =element_text(size=16),
        axis.text = element_text(size=20))+
  xlab("SV Sensitivity Model -log10P-Value")+
  ylab("Random Plate Effect -log10P-Value")+
  geom_abline(intercept=0,slope=1)+
  ggtitle("P-Value Comparison (Random Effect for Plate vs SVA Sensitivity Model)")

png("/home/alreiner/Projects/ffcw/output/BlockProbeResults/SERrandVSERsva.png", width = 12.5, height = 6, units = "in", res = 800)
compsig
dev.off()

#Effect estimate comparison plot
compEE <- ggplot(data=comp1, aes(x=SV_logFC,y=RE_logFC))+
  #geom_jitter(size=.5,width = 0.1)+
  geom_point(size=.5)+
  theme_classic()+
  theme(axis.title =element_text(size=16),
        axis.text = element_text(size=20))+
  xlab("Adjusted % difference in DNA methylation between \nfemale and male children in sensitivity model")+
  ylab("Adjusted % difference in DNA methylation between \nfemale and male children in main model")+
  geom_abline(intercept=0,slope=1)+
  #ggtitle("P-Value Comparison (Random Effect for Plate vs SVA Sensitivity Model)")

png("/home/alreiner/Projects/ffcw/output/BlockProbeResults/SERrandVSERsvaEE.png", width = 8, height = 6, units = "in", res = 800)
#png("/home/alreiner/Projects/ffcw/output/BlockProbeResults/SERrandVSERsvaEE_lr.png", width = 8, height = 6, units = "in", res = 150)
compEE
dev.off()







##########################
#Compare main model with SAND indicator to main model
cityfull <- readRDS(file="/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERfullblockCITY.rds")
citysig <- readRDS(file="/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERsig247blockCITY.rds")

comp1_city <- merge(cityfull %>% select(probe_id,P.Value,logFC),
               age9_RE_full%>%select(probe_id,P.Value,logFC), by="probe_id", all=FALSE)
comp1_city <- comp1_city %>% rename("city_PVal"="P.Value.x",
                          "city_logFC"="logFC.x",
                          "RE_PVal"="P.Value.y",
                          "RE_logFC"="logFC.y")

#calculate correlation between p values
hist(comp1_city$city_PVal) #non-normal
cor(comp1_city$city_PVal, comp1_city$RE_PVal, method="pearson") #0.7543
cor(comp1_city$city_PVal, comp1_city$RE_PVal, method="spearman")#0.8116 #reporting this because variablels
#are not normally ddistributed

#calculate correlation between effect estimates for sex
hist(comp1$SV_logFC,breaks=100)
hist(comp1$RE_logFC)
cor(comp1_city$city_logFC, comp1_city$RE_logFC, method="pearson") #0.9867
cor(comp1_city$city_logFC, comp1_city$RE_logFC, method="spearman")#0.9307

sum(comp1$SV_logFC > 0.1)
sum(comp1$SV_logFC < -0.1)


compsig_city <- ggplot(data=comp1_city, aes(x=-log10(city_PVal),y=-log10(RE_PVal)))+
  #geom_jitter(size=.5,width = 0.1)+
  geom_point(size=.5)+
  theme_classic()+
  theme(axis.title =element_text(size=16),
        axis.text = element_text(size=20))+
  xlab("SV Sensitivity Model -log10P-Value")+
  ylab("Random Plate Effect -log10P-Value")+
  geom_abline(intercept=0,slope=1)+
  ggtitle("P-Value Comparison (Random Effect for Plate vs SVA Sensitivity Model)")

png("/home/alreiner/Projects/ffcw/output/BlockProbeResults/SERrandVSERsva.png", width = 12.5, height = 6, units = "in", res = 800)
compsig_city
dev.off()


