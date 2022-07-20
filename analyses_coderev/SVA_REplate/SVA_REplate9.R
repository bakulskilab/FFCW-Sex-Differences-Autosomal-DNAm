library(tidyverse)
library(limma)
library(sva)

#age 9 phenotype data
analysisdat <- readRDS("/home/alreiner/Projects/ffcw/data/age9final_816.rds")

#age 15 data
analysisdat15 <- readRDS("/home/alreiner/Projects/ffcw/data/age15final_796.rds")

#filter age 9 data to matching age 15 data
analysisdat <- analysisdat %>% filter(id %in% c(analysisdat15$idnum))

# #betas without gap or sex formed in: beta_formation.R
# #beta age 9, no gaps, no XY (n=787)
betanogapnoXY_9 <- readRDS("/home/alreiner/Projects/ffcw/data/age9beta_nogapnoXY_796.rds")
betanogapnoXY_9 <- betanogapnoXY_9[,match(analysisdat$MethID,colnames(betanogapnoXY_9))]
#
# #beta age 15, no gaps, no XY (n=767)
betanogapnoXY_15 <- readRDS("/home/alreiner/Projects/ffcw/data/beta_nogapnoXY_15_796.rds")
betanogapnoXY_15 <- betanogapnoXY_15[,match(analysisdat15$MethID,colnames(betanogapnoXY_15))]

#Try without n.sv input again
sva.fit <- function(beta,pd){
  nullmod <- model.matrix(~1,data=pd)
  fullmod<- model.matrix(~as.factor(cm1bsex)+ Epithelial.cells+as.factor(cm1ethrace), data=pd)
  sva.fit <- sva(beta, fullmod, nullmod)
  sv <- data.frame(sva.fit$sv)
  colnames(sv) <- paste0('sv',1:ncol(sv))
  pd.sv <- cbind(pd,sv)
  return(pd.sv)
}

#Try
age9sva <- sva.fit(betanogapnoXY_9,analysisdat)
write.csv(age9sva, file="/home/alreiner/Projects/ffcw/data/age9dataSV.csv")

age15sva <- sva.fit(betanogapnoXY_15,analysisdat15)
write.csv(age15sva, file="/home/alreiner/Projects/ffcw/data/age15dataSV.csv")

######
