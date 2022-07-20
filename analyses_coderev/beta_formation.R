# BETA MATRICES formulation without gaps or sex chr probes
library(readr)
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


#count numb of autosomal vs x chr sites in the 
sum(rownames(beta) %in% sexprobesX$probe_id) #x chr
sum(!rownames(beta) %in% sexprobes$probe_id) #autosomes 401,545

#load in gapprobes: made in gappys.R
gaps9 <- readRDS("/home/alreiner/Projects/ffcw/output/gapprobes9.rds")
gaps15 <- readRDS("/home/alreiner/Projects/ffcw/output/gapprobes796.rds")

# Age 9 ####################remove sex and gap probes
beta_nogap9 <- beta[!rownames(beta) %in% gaps9,match(age9dat$MethID,colnames(beta))] #left with 406,148 probes
beta_nogap_noXY9_fin <- beta_nogap9[!rownames(beta_nogap9) %in% sexprobes$probe_id,match(age9dat$MethID,colnames(beta_nogap9))] #left with 406,148 probes
#only 405,527 probes with the 15 year old gapes removed and sex probes
saveRDS(beta_nogap_noXY9_fin, file=paste0("/home/alreiner/Projects/ffcw/data/age9beta_nogapnoXY_",nrow(age9dat),".rds"))


#Age 15####################remove sex and gap probes
#no gap 15
beta_nogap15 <- beta[!rownames(beta) %in% gaps15,match(age15dat$MethID,colnames(beta))] #left with 406,148 probes
sum(sexprobes$probe_id %in% rownames(beta_nogap15))

beta_noXYmani15_794 <- beta[!rownames(beta) %in% sexprobes$probe_id,match(age15dat$MethID,colnames(beta))] #414,378 probes left
sum(sexprobes$probe_id %in% rownames(beta_noXYmani15_794)) #no sex probes in our beta_noXY
#beta_noXYmani <- beta[!rownames(beta) %in% sexprobes$probe_id,match(age9fin816$MethID,colnames(beta))] #414,378 probes left
beta_nogap_noXY15_794<- beta_noXYmani15_794[!rownames(beta_noXYmani15_794) %in% gaps15,match(age15dat$MethID,colnames(beta_noXYmani15_794))] # probes left
#now 794 samples, 405,492probes total (matches expected)
saveRDS(beta_nogap_noXY15_794, file=paste0("/home/alreiner/Projects/ffcw/data/beta_nogapnoXY_15_",nrow(age15dat),".rds"))



###############################Remove only sex chr probes

# Age 9 ####################remove sex probes
beta_noXY9_fin <- beta[!rownames(beta) %in% sexprobes$probe_id,match(age9dat$MethID,colnames(beta))] 

saveRDS(beta_noXY9_fin, file=paste0("/home/alreiner/Projects/ffcw/data/age9beta_noXY_",nrow(age9dat),".rds"))


#Age 15####################remove sex probes
beta_noXY15_fin <- beta[!rownames(beta) %in% sexprobes$probe_id,match(age15dat$MethID,colnames(beta))]
saveRDS(beta_noXY15_fin, file=paste0("/home/alreiner/Projects/ffcw/data/beta_noXY_15_",nrow(age15dat),".rds"))


#check overlap of gap and sex chr probes
#length(intersect(sexprobes$probe_id, gapprobes15)) #[1] 3019 --confirmed that this is number of probes that were both gap and sex chr probes
#length(intersect(sexprobes$probe_id, rownames(beta))) #9290
