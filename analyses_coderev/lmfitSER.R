library(ewastools)
library(limma)
library(readr)
library(tidyverse)

########################Goal
# Test out whether the block paramater of lmfit function from limma will allow us to account for 
# the random effect of plate 

###############################Load in needed data
#age 9 phenotype data
analysisdat <- readRDS("/home/alreiner/Projects/ffcw/data/age9final_816.rds")

#age 15 data
analysisdat15 <- readRDS("/home/alreiner/Projects/ffcw/data/age15final_796.rds")

#filter age 9 data to matching age 15 data
analysisdat <- analysisdat %>% filter(id %in% c(analysisdat15$idnum))

#pheno data
pdqc <- read_csv("/nfs/turbo/bakulski1/People/alreiner/analytic/pd_analytic_freeze1.csv") #1684

#full pheno data
fullphen <- readRDS("/nfs/turbo/bakulski1/People/alreiner/data/fullpheno.rds")

#Add additional variables of interest to age 15 data
#childs age in months: ck6yagem (youth age at time of home visit)
analysisdat15 <- left_join(analysisdat15, pdqc %>% select(c("MethID","ck6yagem","cp6povco")), by="MethID")

# #betas without gap or sex formed in: beta_formation.R
# #beta age 9, no gaps, no XY (n=787)
betanogapnoXY_9 <- readRDS("/home/alreiner/Projects/ffcw/data/age9beta_nogapnoXY_796.rds")
betanogapnoXY_9 <- betanogapnoXY_9[,match(analysisdat$MethID,colnames(betanogapnoXY_9))]
#
# #beta age 15, no gaps, no XY (n=767)
# betanogapnoXY_15 <- readRDS("/home/alreiner/Projects/ffcw/data/beta_nogapnoXY_15_796.rds")
# betanogapnoXY_15 <- betanogapnoXY_15[,match(analysisdat15$MethID,colnames(betanogapnoXY_15))]

###############################################

##Age 9 chr 21 results
#########################
# # 1. subset beta to only chr 21 probes for test
# #get chr 21 probes list
# mani <- ewastools:::manifest_450K
# chr21probes_mani <- mani %>% select(probe_id,chr)%>%
#   filter(chr=="21")
# chr21probes <- chr21probes_mani$probe_id
# 
# #subset beta age 9
# chr21beta9 <- betanogapnoXY_9[rownames(betanogapnoXY_9)%in% chr21probes,]

#test that code is doing what i think
# tt <- c("A","B")
# matty <- matrix(data=c(0.2,0.3,0.4,0.5,0.6,0.8),nrow=3,ncol=2)
# rownames(matty)<-c("A","B","C")
# matty
# matty[rownames(matty)%in%tt,]


# Suppose each plate is a block
# blocky <- as.factor(analysisdat$Sample_Plate)
# blocky

#design matrix
# designchr21 <- model.matrix(~as.factor(cm1bsex)+ Epithelial.cells+as.factor(cm1ethrace), data=analysisdat)
# rownames(designchr21) <- colnames(chr21beta9)
# 
# # #calculate correlation
# # dupcor <- duplicateCorrelation(chr21beta9,design=designchr21,block=blocky)
# # #time how long it took to run
# # system.time(duplicateCorrelation(chr21beta9,design=designchr21,block=blocky))
# # #took 96.38 seconds to run for 4243 probes (1min30 secs)
# # # expect 391980 probes to take 8904 seconds (148 minutes or 2.47 hours)
# # dupcor$consensus.correlation
# #
# # #get sig hits by sex for random effects model for Plate
# fitchr21 <- lmFit(chr21beta9,design=designchr21)
# fitchr21 <- eBayes(fitchr21)
# #
# # #top 10 most sig probes
# chr21topres <- topTable(fitchr21,coef=2,adjust.method = "none",p.value = 2.4E-7) #top 10 most sig probes
# chr21topres$probe_id <-rownames(chr21topres)
# 
# #results for all probes tested
# chr21allres <- topTable(fitchr21,coef=2,number=nrow(chr21beta9),adjust.method = "none")
# chr21allres$probe_id <-rownames(chr21allres)
# #
# 
# #read in rand effect for plate results
# randeff <- readRDS("/home/alreiner/Projects/ffcw/output/BlockProbeResults/chr21res.rds")
# #
# # # Make a plot comparing the magnitude of effect estimates for all significant probes
# # # using the original age 9 results (model=sex+batch+race+epicellprop)
# # # and using the new random effect for plate for just chr 21 (model=sex+race+epicellprop+randeffectPlate)
# # #merge data
# comp1 <- merge(randeff %>% select(probe_id,P.Value),
#               chr21allres%>%select(probe_id,P.Value), by="probe_id", all=FALSE)
# comp1 <- comp1 %>% rename("RandPVal"="P.Value.x",
#                         "SERPVal"="P.Value.y")
# 
# #check how many
# sum(comp1$RandPVal > comp1$SERPVal) #2059
# sum(comp1$SERPVal > comp1$RandPVal) #1525
# 
# sum(comp1$RandPVal < 2.4E-7) #74
# sum(comp1$SERPVal < 2.4E-7) #63
# 
# compsig <- ggplot(data=comp1, aes(x=-log10(RandPVal),y=-log10(SERPVal)))+
#   #geom_jitter(size=.5,width = 0.1)+
#   geom_point(size=.5)+
#   theme_classic()+
#   theme(axis.title =element_text(size=16),
#         axis.text = element_text(size=20))+
#   xlab("Random Plate Effect -log10P-Value")+
#   ylab("SER Model -log10P-Value")+
#   geom_abline(intercept=0,slope=1)+
#   ggtitle("P-Value Comparison (Random Effect for Plate vs Not Controlling for Batch)")
# compsig
# 
# compsigzoom <- ggplot(data=comp1, aes(x=-log10(RandPVal),y=-log10(SERPVal)))+
#   #geom_jitter(size=.5,width = 0.1)+
#   geom_point(size=.5)+
#   theme_classic()+
#   theme(axis.title =element_text(size=16),
#         axis.text = element_text(size=20))+
#   xlab("Random Plate Effect -log10P-Value")+
#   ylab("SER Model -log10P-Value")+
#   xlim(0,30)+
#   ylim(0,30)+
#   geom_abline(intercept=0,slope=1)+
#   geom_hline(yintercept=6.62, linetype="dashed", 
#              color = "red", size=0.5)+
#   geom_vline(xintercept=6.62, linetype="dashed", 
#              color = "red", size=0.5)+
#   ggtitle("P-Value Comparison (Random Effect for Plate vs Not Controlling for Batch)")
# compsigzoom
# 
# compsigzoom2 <- ggplot(data=comp1, aes(x=-log10(RandPVal),y=-log10(SERPVal)))+
#   #geom_jitter(size=.5,width = 0.1)+
#   geom_point(size=.5)+
#   theme_classic()+
#   theme(axis.title =element_text(size=16),
#         axis.text = element_text(size=20))+
#   xlab("Random Plate Effect -log10P-Value")+
#   ylab("SER Model -log10P-Value")+
#   xlim(0,8)+
#   ylim(0,8)+
#   geom_abline(intercept=0,slope=1)+
#   geom_hline(yintercept=6.62, linetype="dashed", 
#              color = "red", size=0.5)+
#   geom_vline(xintercept=6.62, linetype="dashed", 
#              color = "red", size=0.5)+
#   ggtitle("P-Value Comparison (Random Effect for Plate vs Not Controlling for Batch)")
# compsigzoom2
# 
# # 
# # compsig1 <- ggplot(data=comp1, aes(x=RandPVal,y=batchPVal))+
# #   geom_jitter(size=.5,width = 0.1)+
# #   theme_classic()+
# #   theme(axis.title =element_text(size=16),
# #         axis.text = element_text(size=20))+
# #   ylim(-.10,0.10)+
# #   xlab("Random Plate Effect P-Value")+
# #   ylab("Batch Effect P-Value")+
# #   geom_abline(intercept=0,slope=1)+
# #   ggtitle("P-Value Comparison (Random Effect for Plate vs Controlling for Batch)")
# # compsig1
# 
# 
# 
# 
# 
# ####################################Same process but FULL beta matrix age 9
# 
# lmfittySER <- function(beta, pd){
#   design<- model.matrix(~as.factor(cm1bsex)+ Epithelial.cells+as.factor(cm1ethrace), data=pd)
#   rownames(design) <- colnames(beta) #change design to modSv
#   #fit models 
#   fit <- lmFit(beta,design=design)
#   fit <- eBayes(fit)
#   
#   #top 10 most sig probes
#   topres <- topTable(fit,coef=2,number=nrow(beta),adjust.method = "none",p.value = 2.4E-7) #top 10 most sig probes
#   topres$probe_id <-rownames(topres)
#   
#   #results for all probes tested
#   allres <- topTable(fit,coef=2,number=nrow(beta),adjust.method = "none")
#   allres$probe_id <-rownames(allres)
#   
#   #decideTests prints direction and number of significant probes
#   print(summary(decideTests(fit, adjust.method="none",p.value = 2.4E-7)))
#   return(list(allres, topres))
# }
# 
# SER9noplate <- lmfittySER(beta=betanogapnoXY_9,pd=analysisdat)
# SER15noplate <- lmfittySER(beta=betanogapnoXY_15,pd=analysisdat15)
# 
# saveRDS(SER9noplate[[1]], file="/home/alreiner/Projects/ffcw/output/age9SERfull.rds")
# saveRDS(SER9noplate[[2]], file="/home/alreiner/Projects/ffcw/output/age9SERsig247.rds")
# 
# saveRDS(SER15noplate[[1]], file="/home/alreiner/Projects/ffcw/output/age15SERfull.rds")
# saveRDS(SER15noplate[[2]], file="/home/alreiner/Projects/ffcw/output/age15SERsig247.rds")


########################################################
#SENSITIVITY RESULTS
#########################################################
SVage9 <- read_csv(file="/home/alreiner/Projects/ffcw/data/age9dataSV.csv")
betanogapnoXY_9 <- betanogapnoXY_9[,match(SVage9$MethID,colnames(betanogapnoXY_9))]


#Run LMFIT but adding 
lmfittySERsv <- function(beta, pd){
  design<- model.matrix(~as.factor(cm1bsex)+ Epithelial.cells+as.factor(cm1ethrace)+sv1+sv2+sv3+sv4+sv5+sv6+sv7+sv8+sv9+sv10, data=pd)
  rownames(design) <- colnames(beta) #change design to modSv
  #fit models 
  fit <- lmFit(beta,design=design)
  fit <- eBayes(fit)
  
  #top 10 most sig probes
  topres <- topTable(fit,coef=2,number=nrow(beta),adjust.method = "none",p.value = 2.4E-7) #top 10 most sig probes
  topres$probe_id <-rownames(topres)
  
  #results for all probes tested
  allres <- topTable(fit,coef=2,number=nrow(beta),adjust.method = "none")
  allres$probe_id <-rownames(allres)
  
  #decideTests prints direction and number of significant probes
  print(summary(decideTests(fit, adjust.method="none",p.value = 2.4E-7)))
  return(list(allres, topres))
}

SER9SV <- lmfittySER(beta=betanogapnoXY_9,pd=SVage9)
#SER15noplate <- lmfittySER(beta=betanogapnoXY_15,pd=analysisdat15)

saveRDS(SER9SV[[1]], file="/home/alreiner/Projects/ffcw/output/age9SERSVfull.rds")
saveRDS(SER9SV[[2]], file="/home/alreiner/Projects/ffcw/output/age9SERSVsig247.rds")

# saveRDS(SER15noplate[[1]], file="/home/alreiner/Projects/ffcw/output/age15SERfull.rds")
# saveRDS(SER15noplate[[2]], file="/home/alreiner/Projects/ffcw/output/age15SERsig247.rds")


