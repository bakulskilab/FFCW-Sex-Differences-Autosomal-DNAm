library(ewastools)
library(limma)
library(readr)
library(tidyverse)
library(dplyr)

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

SANDind <- read.csv("/nfs/turbo/bakulski1/People/alreiner/insand.csv")
SANDidnums <- as.character(SANDind$idnum)

#Age 9 and 15 add city indicator
analysisdat <- analysisdat %>%
  dplyr::mutate(m1city_DetToled = ifelse(m1city==4 |m1city==17,"Yes","No"),
         SAND=ifelse(id %in% c(SANDidnums), "Yes","No"))

analysisdat15 <- analysisdat15 %>%
  dplyr::mutate(m1city_DetToled = ifelse(m1city==4 |m1city==17,"Yes","No"),
         SAND=ifelse(idnum %in% c(SANDidnums), "Yes","No"))

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
betanogapnoXY_15 <- readRDS("/home/alreiner/Projects/ffcw/data/beta_nogapnoXY_15_796.rds")
betanogapnoXY_15 <- betanogapnoXY_15[,match(analysisdat15$MethID,colnames(betanogapnoXY_15))]

###############################################



####################################Same process but FULL beta matrix age 9

lmfittySERcity_block <- function(beta, pd){
  design<- model.matrix(~as.factor(cm1bsex)+ Epithelial.cells+as.factor(cm1ethrace)+factor(SAND), data=pd)
  rownames(design) <- colnames(beta) #change design to modSv
  
  blocky <- as.factor(pd$Sample_Plate)
  
  
  #get correlation between samples on same plate
  dupcor <- duplicateCorrelation(beta,design=design,block=blocky)
  
  #fit models 
  fit <- lmFit(beta,design=design,block=blocky,correlation=dupcor$consensus)
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

SER9block <- lmfittySERcity_block(beta=betanogapnoXY_9,pd=analysisdat)
SER15block <- lmfittySERcity_block(beta=betanogapnoXY_15,pd=analysisdat15)

saveRDS(SER9block[[1]], file="/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERfullblockCITY.rds")
saveRDS(SER9block[[2]], file="/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERsig247blockCITY.rds")

saveRDS(SER15block[[1]], file="/home/alreiner/Projects/ffcw/output/BlockProbeResults/age15SERfullblockCITY.rds")
saveRDS(SER15block[[2]], file="/home/alreiner/Projects/ffcw/output/BlockProbeResults/age15SERsig247blockCITY.rds")


