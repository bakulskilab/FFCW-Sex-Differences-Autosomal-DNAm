#Gap probe identification
library(minfi)
library(readr)
library(tidyverse)


#load beta
beta <- readRDS("/nfs/turbo/bakulski1/People/alreiner/analytic/beta_analytic_freeze1.rds")
#remove cross reactive probes from beta matrix
crossy <- load("/nfs/turbo/bakulski1/cross.probes.info.450k.rda") #29233
crossylist <- as.character(cross.probes.info$TargetID) 
#count num of these probes that are actually in the beta matrix
sum(crossylist %in% rownames(beta)) #27141 probes are in the beta
beta <- beta[!(rownames(beta) %in% crossylist),] #expected count= 410447

#load age 9 data
age9fin <- readRDS("/home/alreiner/Projects/ffcw/data/age9final_816.rds")

#load most datarecent age 15 
age15fin <- readRDS("/home/alreiner/Projects/ffcw/data/age15final_796.rds")

#filter the age9 data to only kids that also have age 15 data
age9fin <- age9fin %>% filter(id %in% c(age15fin$idnum))


#Age 9 gaphunter
#n=787 as of 1/3/22
gaps <- gaphunter(beta[,match(age9fin$MethID,colnames(beta))], threshold=0.05, keepOutliers = FALSE, outCutoff=0.01, verbose=TRUE)
#threshold=The difference in consecutive, ordered beta values that defines the presence of a gap signal. Defaults to 5 percent.
#keepOutliers = FALSE means that outlier -drive gap signals are not included

gapprobes9 <- c(rownames(gaps$proberesults))
saveRDS(gapprobes9, file="/home/alreiner/Projects/ffcw/output/gapprobes9.rds")

####gap hunter 15 year olds
#n=767 as of 1/3/22
assign(paste0("gaps", nrow(age15fin)), 
       gaphunter(beta[,match(age15fin$MethID,colnames(beta))], threshold=0.05, keepOutliers = FALSE, outCutoff=0.01, verbose=TRUE))

assign(paste0("gapprobes", nrow(age15fin)), c(rownames(get(paste0("gaps", nrow(age15fin)))$proberesults)))
saveRDS(get(paste0("gapprobes", nrow(age15fin))), file=paste0("/home/alreiner/Projects/ffcw/output/gapprobes",nrow(age15fin),".rds"))
