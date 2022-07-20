#Get vars from full pheno
#Get vars from pdqc
#Merge
#Turn the variable into NAs
#Whymissing thing
#filter to only cases where whymissing is not NA
library(readr)
library(tidyverse)
library(minfi)
library(ewastools)
library(EpiDISH)
#beta <- readRDS("/nfs/turbo/bakulski1/People/alreiner/data/betaqc.rds") old
beta <- readRDS("/nfs/turbo/bakulski1/People/alreiner/analytic/beta_analytic_freeze1.rds")
fullphen <- readRDS("/nfs/turbo/bakulski1/People/alreiner/data/fullpheno.rds")
pdqc <- read_csv("/nfs/turbo/bakulski1/People/alreiner/analytic/pd_analytic_freeze1.csv") #1684
pdqc_up <- read_csv("/nfs/turbo/bakulski1/People/alreiner/pd_3-1.csv")

####################################SAMPLE CHECK (3/1/22)
#Check how much cm1bsex values change between old pd data and new pd data
pd.qc <- readRDS("/nfs/turbo/bakulski1/People/alreiner/pd.qc.rds") #1824 old

c <- pd.qc %>% filter(childteen=="C") #887 children
t <- pd.qc %>% filter(childteen=="T") #925 teewns
m <- pd.qc %>% filter(childteen=="M") #12 moms

pdqc <- pdqc %>% arrange(as.numeric(idnum))
pdqcid = as.character(pdqc$idnum)
# pdqcmethid = as.character(pdqc$MethID)
# checky <- pd.qc %>% filter(MethID %in% pdqcmethid) %>%
#   arrange(as.numeric(idnum))
# #check how similar sex values are
# cor(checky$cm1bsex,pdqc$cm1bsex) #97.9% correlated sex values (mostly stay the same)
# #how many values are different
# sum(checky$cm1bsex!=pdqc$cm1bsex) #18 people have different sex values between old and new data

checky2 <- pdqc_up %>% filter(methid %in% pdqcmethid) %>%
   arrange(as.numeric(idnum))

# #check how similar sex values are
cor(checky2$cm1bsex,pdqc$cm1bsex) #updated file is concordant with old value
# #how many values are different
sum(checky2$cm1bsex!=pdqc$cm1bsex) #0 people have different sex values between old and new data

#load in updated sample list (1/3/22)
#updatedsamps <- read_csv("/nfs/turbo/bakulski1/People/alreiner/samps_to_keep_12-13.csv")
########################################################################


# get rid of the moms and teens in pdqc
pdqc1 <- pdqc %>%
  dplyr::filter(childteen != "T") %>%
  filter(childteen != "M")%>% 
  mutate(flag_sex = predicted_sex != sex)%>%
  select(MethID,childteen,idnum,m1city,cm1bsex,cm1ethrace,cm1inpov,probe_fail_pct,flag_sex)

#crossdupe = duplicated(paste0(pdqc$childteen,pdqc$idnum)
                       
#n=866 (1 case is not in wave)

#create flag variable for people with discordant sex
idnums_sex_dis <- pdqc1 %>% filter(flag_sex==T)%>%select(idnum) #0 kids
idnums_sex <- c(idnums_sex_dis$idnum)

pdqc1 <- pdqc1 %>%
  group_by(idnum)%>%
  mutate(min = min(probe_fail_pct))%>%
  arrange(probe_fail_pct)%>%
  mutate(crossdupety = duplicated(paste0(childteen,idnum))) #%>%# returns logical T/F if any duplicate C and idnum combinations (paste0--combines C/T and idnum)
  #filter(idnum %in% idnumsdup) #to check that the sample where crossupety=T (which will be removed is the higher probe_pct_fail value)

pqr <- pdqc1 %>%
  filter(crossdupety==T)

# save idnums of these technical duplicates
idnumsdup <- c(pqr$idnum)

#view probe fail pct for these idnums
#lmno <- pdqc1 %>%
  #filter(idnum %in% idnumsdup)%>%
  #select(idnum, probe_fail_pct, crossdupety,childteen)

# select only variables of interest
age9pd <- pdqc1 %>%
  select(MethID,childteen,idnum,m1city,cm1bsex,cm1ethrace,cm1inpov,probe_fail_pct,flag_sex,crossdupety)

# add cell props
# calculate cell proportions
#library(EpiDISH)
#beta with gaps and sex chrom probes still in 
beta_9 <- beta[,match(age9pd$MethID,colnames(beta))] 
#rm(beta)

#Freida way
# beta_9try <- beta[, colnames(beta) %in% age9pd$MethID]
# age9pdT <- age9pd[match(colnames(beta_9try), age9pd$MethID), ]
# all.equal(beta_9,beta_9try) #THESE GIVE DIFFERENT RESULTS!! That's ok--Freida did it in different but equal way
# if(identical(colnames(beta_9), age9pd$MethID)==FALSE){print('STOP: Sample names do not match')}
# identical(colnames(beta_9), age9pd$MethID) #True
# identical(colnames(beta_9try), age9pdT$MethID) #True
# if(dim(beta_9)[2]!=dim(age9pd)[1]){print('STOP: Sample #s do not match')}

#Calculate cell type proportions per sample with EWAS tools
#gaps and sex chr are still included, just rearranged beta cols to match age9pd ordering
cells <- estimateLC(meth=beta_9, ref='saliva',constrained=T)
  # meth= matrix of beta-values for only age 9 kids including gaps and sex chr

#prior CTP alg: 
#cellsold <- epidish(beta.m=beta_9, ref.m=centEpiFibIC.m)

#bivariate plot comparing old cell type alg to new cell type
#plot(cells$Leukocytes, cellsold$estF[,"IC"], xlab="IC Proportion Updated Algorithm",ylab="IC Proportion Old Algorithm")
#plot(cells$Epithelial.cells, cellsold$estF[,"Epi"],xlab="Epithelial Cell Proportion Updated Algorithm",ylab="Epithelial Cell Proportion Old Algorithm")


# merge these new cell proportions with pd
#Assuming: cell type proportion estimates in each row are in the same order as rows in age9pd
test <- cbind(age9pd$MethID, cells)
colnames(test)[1] <- "MethID"
colnames(test)[1]
#colnames(test)

#merge cell type props with pheno data
age9pd_cells <- merge(age9pd, test, by="MethID")


#data.frame(age9pd,cells$estF%>% data.frame) old way

# pdqc vars that we want to merge by
pd <-c("m1city","cm1bsex","cm1ethrace","cm1inpov") #"cm5povco","ck6ethrace"

# full phen vars of interest
age9fp <- fullphen %>%
  select("id","cf1ethrace","cm1edu","cf1edu","hv5_cbmi","hv5_mbmi","IH5mombmi", "m1g1", "m1g4", "k5conf2","hv5_agem", pd)
# add cm1inpov see what happens--merge data has cminpov.x and cominpov.y


# need to keep all variables
mainvars9 <- merge(age9fp,age9pd_cells, by.x=c("id",pd),by.y=c("idnum",pd), all=T) 

# checking for duplicates (all say 25)
#sum(duplicated(mv_9_try$id)) # 25 duplicate ids
#(sum(duplicated(mainvars9_a$id)))
#sum(mainvars9_a$crossdupe, na.rm=TRUE)


# turn all negative #s to NAs (except the k5conf2 variable because we only eliminate cases with -9 value -- not -6 or -3)
mainvars9_a <- mainvars9 %>% mutate_at(vars(-k5conf2),funs(replace(.,.<0,NA)))
#table(mainvars9$k5conf2)

# this might have actually worked
mainvars9_a1 <- mainvars9_a %>%
  mutate(whymissing=ifelse(k5conf2==-9,"no wave5",
                           ifelse(is.na(MethID), "no methyl",
                                  ifelse(crossdupety==TRUE, "duplicate",
                                         ifelse(is.na(cm1bsex),"no sex",
                                                ifelse(flag_sex==T,"sex mismatch",
                                                ifelse(is.na(hv5_agem),"no age",
                                                                               #ifelse(is.na(hv5_ppvtss),"no PPVT",
                                                                               ifelse(is.na(cm1ethrace), "no mom race",
                                                                                      ifelse(is.na(cm1edu),"no mom edu",
                                                                                             ifelse(is.na(hv5_cbmi),"no child BMI", 
                                                                                                    ifelse(is.na(cm1inpov),"no povratio", 
                                                                                                           ifelse(is.na(m1g1),"no mom health", 
                                                                                                                  ifelse(is.na(m1g4), "no mom smoke", NA)))))))))))))

table(mainvars9_a1$whymissing)

age9final <- mainvars9_a1 %>%
  filter(is.na(whymissing)) #now n=816 (without removing kids with no PPVT score)


# loop <- fullsend %>%
#   mutate(whymissing=ifelse(k5conf2==-9,"no wave5",
#                            #ifelse(is.na(MethID), "no methyl",
#                                   #ifelse(crossdupety==TRUE, "duplicate",
#                                          ifelse(is.na(cm1bsex),"no sex",
#                                                 #ifelse(flag_sex==T,"sex mismatch",
#                                                        ifelse(is.na(hv5_agem),"no age",
#                                                               #ifelse(is.na(hv5_ppvtss),"no PPVT",
#                                                               ifelse(is.na(cm1ethrace), "no mom race",
#                                                                      ifelse(is.na(cm1edu),"no mom edu",
#                                                                             ifelse(is.na(hv5_cbmi),"no child BMI", 
#                                                                                    ifelse(is.na(cm1inpov),"no povratio", 
#                                                                                           ifelse(is.na(m1g1),"no mom health", 
#                                                                                                  ifelse(is.na(m1g4), "no mom smoke", NA))))))))))
# 
# table(loop$whymissing)
# 
# looptest <- loop %>%
#   filter(is.na(whymissing)) 
# 
# #REMOVE SAMPLES BASED ON JONAH UPDATE
# #load in updated sample list (1/3/22)
#updatedsamps <- read_csv("/nfs/turbo/bakulski1/People/alreiner/samps_to_keep_12-13.csv")

#methIDs_keep <- c(updatedsamps$MethID2)

#age9fin <- age9final %>% filter(MethID %in% methIDs_keep)

#need to get rid of NAs for epithelial cells since we know they are 0
age9fin <- age9final %>% 
  mutate(Epithelial.cells = ifelse(is.na(Epithelial.cells),
                                   0,Epithelial.cells))

#check that no NAs are present in cell type variables : Leukocytes, Epithelial.cells
sum(is.na(age9fin$Leukocytes)) #no NAs
sum(is.na(age9fin$Epithelial.cells)) #no NAs

#write.csv(age9final, file="/home/alreiner/Projects/ffcw/data/age9final816.csv")
#saveRDS(age9fin, file="/home/alreiner/Projects/ffcw/data/age9fin810_feb2022.rds") #better to use saveRDS to write.csv b/c doesn't add extra col
#age9final <- readRDS("/home/alreiner/Projects/ffcw/data/age9fin787_jan2022.rds")


####################Redo PCA analysis age_9
#NOTE: This will be done in BatchEff.Rmd in analyses folder

######Check how different globalDNAm is between samples with IC (Leukocyte prop) < 0.4
#relative to most other samples with IC prop near 100%
#Hypothesis: will be very different
tt <- age9final %>% filter(Leukocytes < 0.4)
lowICprop <- tt$MethID
colMeans(beta_9[,lowICprop]) #global DNAm is around 44-45% for the low prop IC samples

summary(colMeans(beta_9)) #mean of ALL samples globalDNAM is 47.5%

#Quick check of assoc between CTP and sex: 1-way ANOVA (categorical v contin)
#Assumps are met: 
# no linearity assump
#independence of observations/errors: no observations are correlated with eachother we assume
# normality of errors: check
# equal variance between groups: need to check

sexCTP <- aov(Epithelial.cells~as.factor(cm1bsex),data=age9fin) #1 way ANOVA test
summary(sexCTP)
summary(lm(Epithelial.cells~as.factor(cm1bsex),data=age9fin)) #same as t-test
#p=0.238
#Conc: No significant difference in proportion of epithelial cells by sex

#same exact results for leukocytes b/c dependent values
sexCTP_L <- aov(Leukocytes~as.factor(cm1bsex),data=age9fin) #1 way ANOVA test
summary(sexCTP_L)

#########################Repeat same process but for age 15 sample (age 9 version is in RepeatGlobalDNAm.R)
#age9fin <- readRDS("/home/alreiner/Projects/ffcw/data/age9fin787_jan2022.rds")

# get rid of the moms and children in pdqc
pdqcTeen <- pdqc %>%
  dplyr::filter(childteen != "C") %>%
  filter(childteen != "M")
#n=856

# returns logical T/F if any duplicate C and idnum combinations (paste0--combines C/T and idnum)
pdqcTeen <- pdqcTeen %>% 
  mutate(flag_sex1 = predicted_sex != sex)%>%
  select(MethID,childteen,idnum,m1city,cm1bsex,cm1ethrace,cm1inpov,probe_fail_pct,flag_sex1)

#get IDs of age-9 cohort
nineyoID<-as.character(c(age9fin$id)) #810 IDs
head(nineyoID)
pdqcTeen_match <- pdqcTeen %>% filter(idnum %in% nineyoID) #791 teens have 9 year old data
#9 year olds that have 15 year old data (788 total)
#(note: could be duplicate entries in this data and some 9 year olds may not have 15 year old data)
nineyo_match <- c(pdqcTeen_match$idnum) #791 IDs

#why is there 1 extra idnum at age 15?
#could be a duplicate entry


#Find if any 9 year olds do not have age 15 data
sum(nineyoID %in% nineyo_match) 
#767 of the nrow(age9fin), which is now 787, 9 year olds have 15 year old data (INCLUDING DUPS)

#check if age-9 ids are exact same as age-15 ids that match the age-9 ids
identical(nineyoID,nineyo_match) #not identical

#for duplicate samples, find which had lower probe_fail_pct and keep that one
pdqcTeenA1 <- pdqcTeen_match %>%
  group_by(idnum)%>%
  mutate(min = min(probe_fail_pct))%>%
  arrange(probe_fail_pct)%>%
  mutate(crossdupety = duplicated(paste0(childteen,idnum))) #%>%
  #select(idnum, probe_fail_pct, crossdupety,childteen,min) #n=788

#count number of 15yo id's are duplicated
dups15 <- pdqcTeenA1 %>%
  filter(crossdupety==T) #21 duplicates 
#(788-21 dups = 766 kids in age-9 cohort also had DNAm data at age 15)
# nrow(pdqcTeen_match) - length(dups15)
#788-21 = 767 total samples

sum(dups15 %in% nineyoID) #0 of the 9 year old samples had duplicate samples for their age 15
sum(nineyoID %in% dups15) #0

sum(!age9fin$id %in% pdqcTeenA1$idnum) # 19 age 9 ids are NOT in age 15 sample (810-19=791)
# save idnums of these technical duplicates (n=length(idnumsdup15)) was 29, now 20
idnumsdup15 <- c(dups15$idnum)

#check how many duplicate C,idnum (n=length(idnumsdup15))
sum(duplicated(pdqcTeenA1$idnum)) #21

#view probe fail pct for these idnums
g <- pdqcTeenA1 %>%
  filter(idnum %in% idnumsdup15)%>%
  select(idnum, probe_fail_pct, crossdupety,childteen)

#keep entries where crossupety=F (smaller probe fail pct)
age15fin <- pdqcTeenA1%>%filter(crossdupety==F)
#final check for no duplicates
sum(duplicated(age15fin$idnum))

#write.csv(pdqcTeen_nodup, file="/home/alreiner/Projects/ffcw/data/pdqcTeen_nodup_794.csv", row.names = FALSE)
#saveRDS(age15fin, file=paste0("/home/alreiner/Projects/ffcw/data/age15fin",nrow(age15fin),"_feb2022.rds")) #better to use saveRDS to write.csv b/c doesn't add extra col
#Read in created age 15 data
#age15fin <- readRDS(paste0("/home/alreiner/Projects/ffcw/data/age15fin",nrow(pdqcTeen_nodup),"_jan2022.rds"))


#GAPHUNTER
# found in gappys.R

# BETA MATRICES
# #remake beta matrix without sex chroms or gap probes
# mani <- ewastools:::manifest_450K
# #data(IlluminaHumanMethylation450kmanifest)
# 
# # get a list of cpg sites that are on sex chromosomes
# sexprobes <- mani %>%
#   filter(chr %in% c("X","Y"))%>%
#   select(probe_id)
# 
# #load in gapprobes: gappys.R
# gaps9 <- readRDS("/home/alreiner/Projects/ffcw/output/gapprobes9.rds")
# gaps15 <- readRDS("/home/alreiner/Projects/ffcw/output/gapprobes767.rds")
# 
# # Age 9 ####################remove sex and gap probes
# beta_nogap9 <- beta[!rownames(beta) %in% gaps9,match(age9fin$MethID,colnames(beta))] #left with 406,148 probes
# beta_nogap_noXY9_fin <- beta_nogap9[!rownames(beta_nogap9) %in% sexprobes$probe_id,match(age9fin$MethID,colnames(beta_nogap9))] #left with 406,148 probes
# #only 405,527 probes with the 15 year old gapes removed and sex probes
# saveRDS(beta_nogap_noXY9_fin, "/home/alreiner/Projects/ffcw/data/age9beta_nogapnoXY_816.rds")
# 
# 
# #Age 15####################remove sex and gap probes
# beta_noXYmani15_794 <- beta[!rownames(beta) %in% sexprobes$probe_id,match(age15fin$MethID,colnames(beta))] #414,378 probes left
# sum(sexprobes$probe_id %in% rownames(beta_noXYmani15_794)) #no sex probes in our beta_noXY
# #beta_noXYmani <- beta[!rownames(beta) %in% sexprobes$probe_id,match(age9fin816$MethID,colnames(beta))] #414,378 probes left
# beta_nogap_noXY15_794<- beta_noXYmani15_794[!rownames(beta_noXYmani15_794) %in% gaps15,match(age15fin$MethID,colnames(beta_noXYmani15_794))] # probes left
# #now 794 samples, 405,492probes total (matches expected)
# saveRDS(beta_nogap_noXY15_794, file="/home/alreiner/Projects/ffcw/data/beta_nogapnoXY_15_767.rds")
# 
# #check overlap of gap and sex chr probes
# length(intersect(sexprobes$probe_id, gaps15)) #[1] 3019 --confirmed that this is number of probes that were both gap and sex chr probes
# length(intersect(sexprobes$probe_id, rownames(beta))) #9290
# 

#####Get batch attached
pd.qc <- readRDS("/nfs/turbo/bakulski1/People/alreiner/pd.qc.rds")
library(tidyverse)
library(data.table)
pd <- data.table(pd.qc)
#Slide isolated from MethID
pd[, slide := .(tstrsplit(MethID, "_", keep = 1) %>% unlist)]

#new code from Jonah for estimating batch
starts <- c("200", "201", "399", "9") #all slides start with one of these number combos
pd[, batch := "bmisc"] #initialize batch values as NA or bmisc
for(s in starts) pd[startsWith(slide, s), batch := .(paste0("b", s))]
# for each slide that starts with 200, make the batch variable equal to b200, then do same for each slide starting with 201,399, then 9

#merge pd.qc and analysis1data
age9fin <- left_join(age9fin, pd %>% select(c("MethID","batch")), by="MethID")


#load age 15 final data
age15fin <- left_join(age15fin, pd %>% select(c("MethID","batch")), by="MethID")

#get beta_15 ready
beta_15 <- beta[,match(age15fin$MethID,colnames(beta))] 

#Compute cell type proportion for age 15 data
#Calculate cell type proportions per sample with EWAS tools
#gaps and sex chr are still included, just rearranged beta cols to match age9pd ordering
cells15 <- estimateLC(meth=beta_15, ref='saliva',constrained=T)
# meth= matrix of beta-values for only age 15 kids including gaps and sex chr

# merge these new cell proportions with pd
#Assuming: cell type proportion estimates in each row are in the same order as rows in age9pd
test15 <- cbind(age15fin$MethID, cells15)
colnames(test15)[1] <- "MethID"
colnames(test)[1]
#colnames(test)

#merge cell type props with pheno data
age15fin <- merge(age15fin, test15, by="MethID")


#REAL BATCH variable
#read.csv(unzip("/nfs/turbo/bakulski1/People/alreiner/batchfiles20210416144155.zip"))
Methylation_450K_array_batch_information <- read_csv("/home/alreiner/Projects/ffcw/data/Methylation_450K_array_batch_information.csv")

age9fin <- left_join(age9fin, Methylation_450K_array_batch_information, by="MethID")
age15fin <- left_join(age15fin, Methylation_450K_array_batch_information, by="MethID")

saveRDS(age9fin, file=paste0("/home/alreiner/Projects/ffcw/data/age9final_",nrow(age9fin),".rds"))

saveRDS(age15fin, file=paste0("/home/alreiner/Projects/ffcw/data/age15final_",nrow(age15fin),".rds"))

#factor necessary variables
# names=c('id','m1city','cm1bsex','cm1ethrace','cf1ethrace','cm1edu','cf1edu','m1g1','m1g4','k5conf2')
# age9fin[,names]<-lapply(age9fin[,names],factor)
# age15fin[,names]<-lapply(age15fin[,names],factor)


# ####Get PCs attached
# #Make PCA plot to show if there are clusters in global DNAm by plate
library(ggfortify)
pcdata <- readRDS("/home/alreiner/Projects/ffcw/data/age9beta_nogapnoXY_796.rds")
pcdata15 <- readRDS("/home/alreiner/Projects/ffcw/data/beta_nogapnoXY_15_796.rds")
#Princ comppnents calculation
pca <- prcomp(t(pcdata),scale. = F,center=T)
pca15 <- prcomp(t(pcdata15),scale. = F,center=T)
#make rownames into column in data
pca$x <- as.data.frame(pca$x) %>% mutate(MethID=as.character(rownames(pca$x[,1:20])))
pca15$x <- as.data.frame(pca15$x) %>% mutate(MethID=as.character(rownames(pca15$x[,1:20])))
#head(pca$x)
#View(as.data.frame(pca$x[,c(1:20,817)]))
#summary(pca)

#calculate variance explained by top 2 PCs age 9
var_explained <- pca$sdev^2/sum(pca$sdev^2) #stdev of the given PC squared / sum of squared sdevs of all PCs
var_explained[1:5] #PC1 explains 17.94%, PC2 explains 13.02%

#calculate variance explained by top 2 PCs age 15
var_explained15 <- pca15$sdev^2/sum(pca15$sdev^2) #stdev of the given PC squared / sum of squared sdevs of all PCs
var_explained15[1:5] #

#merge PCA data
#temp load in age 9 fin data
analysisdat <- readRDS("/home/alreiner/Projects/ffcw/data/age9final_816.rds")
analysisdat15 <- readRDS("/home/alreiner/Projects/ffcw/data/age15final_796.rds")
analysisdat <- analysisdat %>% filter(id %in% c(analysisdat15$idnum))

age9finPC <- merge(analysisdat, as.data.frame(pca$x[,c(1:20)]),by.x="MethID", by.y="row.names")

age15finPC <- merge(analysisdat15, as.data.frame(pca15$x[,c(1:20)]),by.x="MethID", by.y="row.names")

#add global DNAm data **source: RepeatGlobalDNAm.R
#age 9
globby9 <- read_csv("/home/alreiner/Projects/ffcw/data/globalDNAm816.csv")
age9finPC <- left_join(age9finPC, globby9 %>% select(c("MethID","globALL")), by="MethID")

#age 15
globby <- read_csv("/home/alreiner/Projects/ffcw/data/globalDNAm796_15.csv")
age15finPC <- left_join(age15finPC, globby %>% select(c("MethID","globALL15")), by="MethID")

saveRDS(age9finPC, file=paste0("/home/alreiner/Projects/ffcw/data/age9final_PC_",nrow(age9finPC),".rds"))

saveRDS(age15finPC, file=paste0("/home/alreiner/Projects/ffcw/data/age15final_PC_",nrow(age15finPC),".rds"))
