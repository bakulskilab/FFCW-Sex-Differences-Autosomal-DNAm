# BETA MATRICES formulation without gaps or sex chr probes
library(readr)
library(limma)
library(ewastools)
library(tidyverse)

#load og beta matrix
#old beta : readRDS("/nfs/turbo/bakulski1/People/alreiner/data/betaqc.rds")
beta <- readRDS("/nfs/turbo/bakulski1/People/alreiner/analytic/beta_analytic_freeze1.rds")
#437588 probes

#remove cross reactive probes from beta matrix
crossy <- load("/nfs/turbo/bakulski1/cross.probes.info.450k.rda") #29233
crossylist <- as.character(cross.probes.info$TargetID) 

#count num of these probes that are actually in the beta matrix
sum(crossylist %in% rownames(beta)) #27141 probes are in the beta
beta <- beta[!(rownames(beta) %in% crossylist),] #expected count= 410447

#load most recent age 9 data
age9dat <- readRDS("/home/alreiner/Projects/ffcw/data/age9final_816.rds")
#load most datarecent age 15 
age15dat <- readRDS("/home/alreiner/Projects/ffcw/data/age15final_796.rds")

#filter the age9 data to only kids that also have age 15 data
age9dat <- age9dat %>% filter(id %in% c(age15dat$idnum))

#load illumina 450k probe info to get sex chr probes
mani <- ewastools:::manifest_450K
#data(IlluminaHumanMethylation450kmanifest)

# get a list of cpg sites that are on sex chromosomes
sexprobes <- mani %>%
  filter(chr %in% c("X","Y"))%>%
  select(probe_id)

sexprobesX <- mani %>%
  filter(chr %in% c("X"))%>%
  select(probe_id)

xprobes <- c(sexprobesX$probe_id)


#plot indiivdual probes
betax = beta[rownames(beta)%in%xprobes,match(age9dat$MethID,colnames(beta))]


lmfittySER_block <- function(beta, pd){
  design<- model.matrix(~as.factor(cm1bsex)+ Epithelial.cells+as.factor(cm1ethrace), data=pd)
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

SER9blockX <- lmfittySER_block(beta=betax,pd=age9dat)

saveRDS(SER9blockX[[1]], file="/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERfullblock_Xchr.rds")
saveRDS(SER9blockX[[2]], file="/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERsig247block_Xchr.rds")


xchrfull <- readRDS("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERfullblock_Xchr.rds")
sum(xchrfull$P.Value < 2.4E-7)/nrow(xchrfull) #91.6% of probes differentially methylated by sex

#Prettify version for paper
xchrfull_paper <- xchrfull %>% mutate(logFCpct=logFC*100)%>% select(probe_id,logFC,logFCpct,AveExpr,t,P.Value)
write.csv(xchrfull_paper,"/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERfullblock_Xchr_manu.csv")
#############################
#EDA#
#############################
#pick single x chr probe
# x1_betas <-betax[1,]
# age9dat$x1betas = x1_betas
# age9dat$x1500betas = betax[1500,]
# age9dat$x8000betas = betax[8000,]
# 
# age9dat$cm1bsex <- factor(age9dat$cm1bsex)
# levels(age9dat$cm1bsex) <- c("Male","Female")
# 
# #side by side boxplots of females vs males
# pr1 <- ggplot(data=age9dat, aes(x=factor(cm1bsex), y=x1betas))+
#   geom_boxplot()+
#   theme(text = element_text(size=20))+
#   xlab("Sex")
# pr1
# #side by side boxplots of females vs males
# pr1500 <- ggplot(data=age9dat, aes(x=factor(cm1bsex), y=x1500betas))+
#   geom_boxplot()+
#   theme(text = element_text(size=20))+
#   xlab("Sex")
# 
# pr8000 <- ggplot(data=age9dat, aes(x=factor(cm1bsex), y=x8000betas))+
#   geom_boxplot()+
#   theme(text = element_text(size=20))+
#   xlab("Sex")
# 
# require(gridExtra)
# grid.arrange(pr1, pr1500, pr8000,ncol=1)
# 
# #check how often we see females have beta values dip below 0.5
# 
# #go through each proobe in x chr and check
# #if at cleast 1 beta value for females goes below 0.5
# #subset beta to gals only 
# methidgirls <- age9dat %>% filter(cm1bsex=="Female")%>%select(MethID)
# galmethid <- methidgirls$MethID
# 
# m[1]
# m = lapply(1:8896, function(x) betax[x,match(methidgirls$MethID,colnames(betax))]<0.5)
# summs = lapply(m,function(x) sum(x,na.rm=T)) #number of people with beta < 0.5 for that probe
# #for probe 1 on x chr, all gals have a value less than 0.5
# nums <- unlist(summs)
# sum(nums==0, na.rm=T) #num of TRUES for which 0 females have less than 0.5
# nums
# la <- c(TRUE,TRUE,FALSE)
# sum(la, na.rm=T)
# #subset beta to only girls
# 
# #loop through each probe and determine if any people have less than 0.5


