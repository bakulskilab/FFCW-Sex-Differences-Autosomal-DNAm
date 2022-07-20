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
# #
# # #subset beta age 9
#  chr21beta9 <- betanogapnoXY_9[rownames(betanogapnoXY_9)%in% chr21probes,]

#test that code is doing what i think
# tt <- c("A","B")
# matty <- matrix(data=c(0.2,0.3,0.4,0.5,0.6,0.8),nrow=3,ncol=2)
# rownames(matty)<-c("A","B","C")
# matty
# matty[rownames(matty)%in%tt,]


# # Suppose each plate is a block
# blocky <- as.factor(analysisdat$Sample_Plate)
# # blocky
# 
# #design matrix
# designchr21 <- model.matrix(~as.factor(cm1bsex)+ Epithelial.cells+as.factor(cm1ethrace), data=analysisdat)
# rownames(designchr21) <- colnames(chr21beta9)
# 
# #calculate correlation
# dupcor <- duplicateCorrelation(chr21beta9,design=designchr21,block=blocky)
# #time how long it took to run
# #system.time(duplicateCorrelation(chr21beta9,design=designchr21,block=blocky))
# #took 96.38 seconds to run for 4243 probes (1min30 secs)
# # expect 391980 probes to take 8904 seconds (148 minutes or 2.47 hours)
# dupcor$consensus.correlation
# 
# #get sig hits by sex for random effects model for Plate
# fitchr21 <- lmFit(chr21beta9,design=designchr21,block=blocky,correlation=dupcor$consensus)
# fitchr21 <- eBayes(fitchr21)
# #
# # #top 10 most sig probes
# chr21topres <- topTable(fitchr21,coef=2,adjust.method = "none",p.value = 2.4E-7) #top 10 most sig probes
# chr21topres$probe_id <-rownames(chr21topres)
# 
# #results for all probes tested
# chr21allres <- topTable(fitchr21,coef=2,number=nrow(chr21beta9),adjust.method = "none")
# chr21allres$probe_id <-rownames(chr21allres)

#
# #compare results to the chr 21 probes when did main model adjusting for batch!!
# #loaad in old results
# age9totalres <- readRDS("/home/alreiner/Projects/ffcw/output/SigProbes/age9_SEBR_247_8324.rds")
#
# #filter down to chr 21 probes
# age9totalres$probe_id <- rownames(age9totalres)
#
# #get only sig chr 21 probes from original main model results
# age9sig21 <- age9totalres %>% filter(probe_id %in% chr21probes) #69 sig probes are on chr 21
# sum(chr21allres$P.Value < 2.4E-7) #74 sig probes are found on chr 21 using new blocked method
#
# #get sig probes from new rand effect model for Plate
# sigblockchr21 <- chr21allres %>% filter(P.Value < 2.4E-7)
#
# sum(age9sig21$probe_id %in% sigblockchr21$probe_id) #63 probes in common
# #save results
# #saveRDS(,"/home/alreiner/Projects/ffcw/output/BlockProbeResults/chr21res.rds)


#saveRDS(chr21allres,"/home/alreiner/Projects/ffcw/output/BlockProbeResults/chr21res.rds")

# # Make a plot comparing the magnitude of effect estimates for all significant probes 
# # using the original age 9 results (model=sex+batch+race+epicellprop)
# # and using the new random effect for plate for just chr 21 (model=sex+race+epicellprop+randeffectPlate)
# #merge data
# comp1 <- merge(sigblockchr21 %>% select(probe_id,P.Value), 
#               age9sig21%>%select(probe_id,P.Value), by="probe_id", all=FALSE)
# comp1 <- comp1 %>% rename("RandPVal"="P.Value.x",
#                         "batchPVal"="P.Value.y")
# 
# 
# compsig <- ggplot(data=comp1, aes(x=-log10(RandPVal),y=-log10(batchPVal)))+
#   #geom_jitter(size=.5,width = 0.1)+
#   geom_point(size=.5)+
#   theme_classic()+
#   theme(axis.title =element_text(size=16),
#         axis.text = element_text(size=20))+
#   xlab("Random Plate Effect -log10P-Value")+
#   ylab("Batch Effect -log10P-Value")+
#   geom_abline(intercept=0,slope=1)+
#   ggtitle("P-Value Comparison (Random Effect for Plate vs Controlling for Batch)")
# compsig
# 
# compsig1 <- ggplot(data=comp1, aes(x=RandPVal,y=batchPVal))+
#   geom_jitter(size=.5,width = 0.1)+
#   theme_classic()+
#   theme(axis.title =element_text(size=16),
#         axis.text = element_text(size=20))+
#   ylim(-.10,0.10)+
#   xlab("Random Plate Effect P-Value")+
#   ylab("Batch Effect P-Value")+
#   geom_abline(intercept=0,slope=1)+
#   ggtitle("P-Value Comparison (Random Effect for Plate vs Controlling for Batch)")
# compsig1





####################################Same process but FULL beta matrix age 9

lmfittySER_block <- function(beta, pd){
  design<- model.matrix(~as.factor(cm1bsex)+ Epithelial.cells+as.factor(cm1ethrace), data=pd)
  rownames(design) <- colnames(beta) #change design to modSv

  blocky <- as.factor(pd$Sample_Plate)


  #get correlation between samples on same plate
  dupcor <- duplicateCorrelation(beta,design=design,block=blocky)

  #fit models
  fit <- lmFit(beta,design=design,block=blocky,correlation=dupcor$consensus)
  fitty <- eBayes(fit)
  
  #top 10 most sig probes
  # topres <- topTable(fit,coef=2,number=nrow(beta),adjust.method = "none",p.value = 2.4E-7) #top 10 most sig probes
  # topres$probe_id <-rownames(topres)
  # 
  # #results for all probes tested
  # allres <- topTable(fit,coef=2,number=nrow(beta),adjust.method = "none")
  # allres$probe_id <-rownames(allres)

  #decideTests prints direction and number of significant probes
  #print(summary(decideTests(fit, adjust.method="none",p.value = 2.4E-7)))
  #return(list(allres, topres))
  return(list(fit,fitty))
}

SER9block <- lmfittySER_block(beta=betanogapnoXY_9,pd=analysisdat)
# SER15block <- lmfittySER_block(beta=betanogapnoXY_15,pd=analysisdat15)
# 
saveRDS(SER9block[[1]], file="/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERlmfitoutputBLOCK.rds")
saveRDS(SER9block[[2]], file="/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SEReBayesoutputBLOCK.rds")


# saveRDS(SER9block[[1]], file="/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERfullblock.rds")
# saveRDS(SER9block[[2]], file="/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERsig247block.rds")
# 
# saveRDS(SER15block[[1]], file="/home/alreiner/Projects/ffcw/output/BlockProbeResults/age15SERfullblock.rds")
# saveRDS(SER15block[[2]], file="/home/alreiner/Projects/ffcw/output/BlockProbeResults/age15SERsig247block.rds")


#Check the df.prior
lmfitout <- readRDS("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERlmfitoutputBLOCK.rds")
ebayesout <- readRDS("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SEReBayesoutputBLOCK.rds")
dgdfs <- lmfitout$df.residual
dodfs <- ebayesout$df.prior
ebayesout$df.total

2*pt(-102,df=791)

?pt

pt(102.39293,df=790+dodfs)
