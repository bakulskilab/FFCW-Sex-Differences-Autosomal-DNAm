#Repeat global DNA analysis for verifcation and with updated data

#Load beta without sex probes
library(readr)
library(tidyverse)
library(ewastools)

#load analytic sample
age9fin816PC <- readRDS("/home/alreiner/Projects/ffcw/data/age9final_816.rds")
#age9fin816PC <- age9fin816PC[,-c(1,3)]

#age 15 data
age15pd <-readRDS("/home/alreiner/Projects/ffcw/data/age15final_796.rds")
age9fin816PC <- age9fin816PC %>% filter(id %in% c(age15pd$idnum))


#load beta matrix without sex chroms probes (age 9)
beta_noXY <- readRDS("/home/alreiner/Projects/ffcw/data/age9beta_noXY_816.rds")
beta_noXY <- beta_noXY[,match(age9fin816PC$MethID,colnames(beta_noXY))] #make sure beta in same order as data

#load beta matrix without sex chroms probes (age 15)
beta_noXY15 <- readRDS("/home/alreiner/Projects/ffcw/data/beta_noXY_15_796.rds")
beta_noXY15 <- beta_noXY15[,match(age15pd$MethID,colnames(beta_noXY15))] #make sure beta in same order as data

#Make sure beta is in same order as analysisdata
# beta <- readRDS("/nfs/turbo/bakulski1/People/alreiner/data/betaqc.rds") #423,668 probes
# beta_noXYmani <- beta[!rownames(beta) %in% sexprobes$probe_id,match(age9fin816$MethID,colnames(beta))] #414,378 probes left
# sum(sexprobes$probe_id %in% rownames(beta_noXYmani)) #no sex probes in our beta_noXY
# #beta_noXYmani <- beta[!rownames(beta) %in% sexprobes$probe_id,match(age9fin816$MethID,colnames(beta))] #414,378 probes left
# 
# #sanity check
# bet <- c("al","br","ca")
# nothx <- c("al")
# !bet %in% nothx # get F, T,T (only keep beta probes, not sex probes)
# 
# #now 816 samples, n=414,378 probes (still gaps, no sex chr)
# #expect 412,020 probes....
# 
# #Test hypothesis: Only 9290 probe difference b/w og beta and beta_noXY (even though 11,648 sex probes found in mani)
# # 2358 probes in sexprobes are not found in our og beta (removed in QC already)
# sum(!sexprobes$probe_id %in% rownames(beta)) #=2358 check

# ########Save beta with no gaps no sex chr
# beta_nogap_noXYmani_816 <- beta_noXYmani[!rownames(beta_noXYmani) %in% gapprobes816,match(age9fin816$MethID,colnames(beta_noXYmani))] # probes left
# #now 816 samples, 406,114 probes total (matches expected)
# saveRDS(beta_nogap_noXYmani_816, file="data/beta_nogapnoXY816.rds")
# 
# #check overlap of gap and sex chr probes
# length(intersect(sexprobes$probe_id, gapprobes816)) #[1] 2998 --confirmed that this is number of probes that were both gap and sex chr probes


#calculate globalDNAm with colMeans for each sample (no sex chr, INCLUDING gap)
globALL <- colMeans(beta_noXY)*100
#View(globALL)

#add globalDNAm function 
#every probe site in beta matrix has entry in Islands.UCSC--know region of every site--NOT chrom location
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data("Islands.UCSC") #485,512 probes
#View(Islands.UCSC)
genomic.region <- function(X, region, anno='450k'){
  #x=matrix with cpg as row names
  #region=area of interest
  #anno= 450k or EPIC
  
  #filter probe names of Island to the ones in beta_noXYmani
  Islands.UCSC <- Islands.UCSC[rownames(X),]
  #filter probes to only those where relation to island is the given region
  Islands.UCSC <- Islands.UCSC[Islands.UCSC$Relation_to_Island==region,]
  #return beta matrix for the probes in given region
  return(X[rownames(Islands.UCSC),])
}

unique(Islands.UCSC$Relation_to_Island)
#Check that no empty slots in Islands.UCSC
table(Islands.UCSC$Relation_to_Island) #adds to 485,512
sum(Islands.UCSC$Relation_to_Island==" ")


# # compute means by region
meanDNAm.sea <- colMeans(genomic.region(beta_noXY, 'OpenSea'))*100 #this gives same results from last time!

meanDNAm.shoren <- colMeans(genomic.region(beta_noXY,'N_Shore')) * 100
meanDNAm.shores <- colMeans(genomic.region(beta_noXY,'S_Shore')) * 100
meanDNAm.shelfn <- colMeans(genomic.region(beta_noXY,'N_Shelf')) * 100
meanDNAm.shelfs <- colMeans(genomic.region(beta_noXY,'S_Shelf')) * 100
#meanDNAm.shore <- colMeans(genomic.region(beta,'Shore')) * 100
#meanDNAm.shelf <- colMeans(genomic.region(beta,'Shelf')) * 100
meanDNAm.island <- colMeans(genomic.region(beta_noXY,'Island')) * 100

#make table and merge with analysisdata
globalM <- data.frame(globALL,meanDNAm.sea,meanDNAm.shoren,meanDNAm.shores,
                      meanDNAm.shelfn,meanDNAm.shelfs,meanDNAm.island)
head(rownames(globalM))

#label each globalDNAm calc with correct id names
globalM <- globalM %>% mutate(MethID=rownames(globalM))

# merge pheno & global methylation values (was pdqcA)
GlobMdata <- left_join(age9fin816PC %>% select(MethID,cm1bsex,Leukocytes,Epithelial.cells),globalM,by ="MethID") # only child age 9 pheno & global meth values (including duplicate kids)
write.csv(GlobMdata, file="/home/alreiner/Projects/ffcw/data/globalDNAm816.csv")

# #Read it in
# globby <- read_csv("/home/alreiner/Projects/ffcw/data/globalDNAm816.csv")
# 
# #check if the distrib of males and female global DNAm is normal or not
# malesM <- GlobMdata %>% filter(cm1bsex==1)
# hist(malesM$globALL) # slightly skewed
# boxplot(malesM$globALL)
# qqplot(malesM$globALL)
# #calc avg for males
# mean(malesM$globALL) #0.4700985
# sd(malesM$globALL)
# 
# femsM <- GlobMdata %>% filter(cm1bsex==2)
# hist(femsM$globALL) #slightly skewed
# #calc avg for females
# mean(femsM$globALL) #0.4719025
# sd(femsM$globALL)
# 
# #if normal, use 2 sample t-test to determine if difference between males and females
# #if non-normal, use Mann-Whitney nonpparametric test (verification)
# 
# # independent 2-group Mann-Whitney U Test
# wilcox.test(globALL~factor(cm1bsex), data=globby)
# t.test(globALL~as.factor(cm1bsex), data=GlobMdata)
# # where y is numeric and A is A binary factor
# 
# #Results: W=69,631 p=0.0003116
# 
# #Test if the difference is significant across different t-tests
# hist(malesM$meanDNAm.sea)
# hist(femsM$meanDNAm.sea)
# 
# #mann for globalDNAm diff by sex across other regions
# wilcox.test(meanDNAm.sea~factor(cm1bsex), data=GlobMdata) #p=0.019
# wilcox.test(meanDNAm.island~factor(cm1bsex), data=GlobMdata) #p=3.43 e -5
# wilcox.test(meanDNAm.shoren~factor(cm1bsex), data=GlobMdata) #p=3.13 E -6
# wilcox.test(meanDNAm.sea~factor(cm1bsex), data=GlobMdata) #p=0.019
# 
# #make images
# # mean globalm by sex
# View(globby %>% group_by(cm1bsex)%>%
#        summarise(avgm = mean(globALL),
#                  sdm=sd(globALL),
#                  avgseam = mean(meanDNAm.sea),
#                  sdsea = sd(meanDNAm.sea),
#                  avgislandm = mean(meanDNAm.island),
#                  sdisland = sd(meanDNAm.island),
#                  avgshoren = mean(meanDNAm.shoren),
#                  sdshoren = sd(meanDNAm.shoren),
#                  avgshore_s = mean(meanDNAm.shores),
#                  sdshore_s = sd(meanDNAm.shores),
#                  avgshelfn = mean(meanDNAm.shelfn),
#                  sdshelfn = sd(meanDNAm.shelfn),
#                  avgshelf_s = mean(meanDNAm.shelfs),
#                  sdshelf_s = sd(meanDNAm.shelfs)
#        ))
# 
# #Make new ggplot for updated globalDNAm by sex
# library(ggplot2)
# theme_update(text = element_text(size=15))
# 
# #plot of sex differences in global methylation data
# ggplot(data=GlobMdata, aes(y= globALL, x=factor(cm1bsex), fill=factor(cm1bsex) ))+ 
#   geom_boxplot() + 
#   ylim(40,55)+
#   ylab("Autosomal Global DNAm (%)")+
#   xlab("Sex")+
#   scale_x_discrete(labels=c("Male (n=414)","Female (n=402)"))+
#   theme_bw()+
#   theme(legend.position = "none")+ 
#   theme(axis.text.x = element_text(face = "bold",size = 12))+
#   scale_fill_brewer(palette="BuPu") + 
#   ggtitle("Sex differences in autosomal global DNAm")

#Make old ggplot of the by sex boxplots for all regions, shores, x, y etc



#####################################Age 15 Analysis#######################
#In 15.yo in 15yo folder!

#calculate globalDNAm with colMeans for each sample (no sex chr, INCLUDING gap)
globALL15 <- colMeans(beta_noXY15)*100

# # compute means by region
meanDNAm.sea15 <- colMeans(genomic.region(beta_noXY15, 'OpenSea'))*100 #this gives same results from last time!

meanDNAm.shoren15 <- colMeans(genomic.region(beta_noXY15,'N_Shore')) * 100
meanDNAm.shores15 <- colMeans(genomic.region(beta_noXY15,'S_Shore')) * 100
meanDNAm.shelfn15 <- colMeans(genomic.region(beta_noXY15,'N_Shelf')) * 100
meanDNAm.shelfs15 <- colMeans(genomic.region(beta_noXY15,'S_Shelf')) * 100
#meanDNAm.shore <- colMeans(genomic.region(beta,'Shore')) * 100
#meanDNAm.shelf <- colMeans(genomic.region(beta,'Shelf')) * 100
meanDNAm.island15 <- colMeans(genomic.region(beta_noXY15,'Island')) * 100

#make table and merge with analysisdata
globalM15 <- data.frame(globALL15,meanDNAm.sea15,meanDNAm.shoren15,meanDNAm.shores15,
                        meanDNAm.shelfn15,meanDNAm.shelfs15,meanDNAm.island15)
head(rownames(globalM15))

#label each globalDNAm calc with correct id names
globalM15 <- globalM15 %>% mutate(MethID=rownames(globalM15))

# merge pheno & global methylation values (was pdqcA)
GlobMdata15 <- left_join(age15pd %>% select(MethID,cm1bsex,Leukocytes,Epithelial.cells),globalM15,by ="MethID") # only child age 9 pheno & global meth values (including duplicate kids)
write.csv(GlobMdata15, file="/home/alreiner/Projects/ffcw/data/globalDNAm796_15.csv")

#globby <- read_csv("/home/alreiner/Projects/ffcw/data/globalDNAm794_15.csv",rownames.included=F)

# #how many boys v girls in age 15 data
# table(GlobMdata15$cm1bsex) #1 boy, 2 girl
# #median and iqr globalDNAm by sex
# GlobMdata15 %>% group_by(as.factor(cm1bsex))%>%summarise(mglob = median(globALL15),
#                                          per1 = quantile(globALL15,probs=0.25),
#                                          per2 = quantile(globALL15,probs=0.75)
# )
# #test for diff in means (wilcox test for skewed dist)
# wilcox.test(globALL15~factor(cm1bsex), alternative = "two.sided", data=GlobMdata15) #0.00965
# 
# #Check global methyl diff by sex age 15
# View(GlobMdata15 %>% group_by(cm1bsex)%>%
#        summarise(avgm = mean(globALL15),
#                  sdm=sd(globALL15)))
# 
# 
# # independent 2-group Mann-Whitney U Test
# wilcox.test(globALL15~factor(cm1bsex), data=GlobMdata15)
# t.test(globALL~as.factor(cm1bsex), data=GlobMdata)
# # where y is numeric and A is A binary factor
# 
# 
# 
# #####################Make plot of number of significant probes that map to each 
# library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
# annotation.table = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# data(Locations)
# View(Locations)
# #data(Islands.UCSC)
# data(Other)
# View(Other)
# 
# #pos<-Locations[rownames(betaqc),'pos']
# #chrnames <- as.character(Locations[rownames(betaqc),'chr'])
# 
# add_annotations<-function(data){
#   #add annotations
#   Locations <- Locations[rownames(data), ]
#   Islands.UCSC <- Islands.UCSC[rownames(data),]
#   Other <- Other[rownames(data),]
#   # add columns of interest from annotation to annotated version
#   data.annotated <- data
#   data.annotated$chr <- Locations$chr
#   data.annotated$pos <- Locations$pos
#   data.annotated$island <- Islands.UCSC$Relation_to_Island
#   data.annotated$gene <- Other$UCSC_RefGene_Name
#   data.annotated$gene_group <- Other$UCSC_RefGene_Group
#   data.annotated$regulatory_feature <- Other$Regulatory_Feature_Group
#   return(data.annotated)
# }
# 
# #pull out top hits
# TopHits_SEBR_247 <- readRDS("/home/alreiner/Projects/ffcw/output/SigProbes/age9_SEBR_247_11136.rds")
# TopHits_SEBR_full <- readRDS("/home/alreiner/Projects/ffcw/output/SigProbes/age9_SEBR_full_updatedCTP.rds")
# # Annotate top hits
# age9.SEBR.hits.annotated <- add_annotations(TopHits_SEBR_247)
# write_csv(age9.SEBR.hits.annotated, file="/home/alreiner/Projects/ffcw/data/AnnotatedLmFitHits/age9_SEBR_247results_annot_upCTP.csv")
# obs <- read_csv("/home/alreiner/Projects/ffcw/data/AnnotatedLmFitHits/age9_SEBR_247results_annot_upCTP.csv")
# view(obs %>% filter(chr=="chr17"))
# 
# View(table(age9.SEBR.hits.annotated.full$chr))
# 
# #Make a for loop that makes 2x2 table for chisquare test for each CpG island region
# 
# 
# 
# # Annotate top hits
# age9.SEBR.hits.annotated.full <- add_annotations(TopHits_SEBR_full)
# write_csv(age9.SEBR.hits.annotated.full, file="/home/alreiner/Projects/ffcw/data/AnnotatedLmFitHits/age9_SEBR_fullresults_annot_upCTP.csv")
# age9.SEBR.hits.annotated.full <- age9.SEBR.hits.annotated.full %>% 
#   mutate(IsSig = ifelse(rownames(age9.SEBR.hits.annotated.full) %in% rownames(TopHits_SEBR_247),"Yes","No"),
#          IsIsland = ifelse(island=="Island","Island","Not Island"),
#          IsOpenSea = ifelse(island=="OpenSea","OpenSea","Not OpenSea"),
#          IsNShelf = ifelse(island=="N_Shelf","N_Shelf","Not NShelf"),
#          IsSShelf = ifelse(island=="S_Shelf","S_Shelf","Not SShelf"),
#          IsNShore = ifelse(island=="N_Shore","N_Shore","Not NShore"),
#          IsSShore = ifelse(island=="S_Shore","S_Shore","Not SShore")
#          )
# 
# #Create table of sig hits by island region
# table(age9.SEBR.hits.annotated.full$IsIsland, age9.SEBR.hits.annotated.full$IsSig)
# chiIS <- chisq.test(age9.SEBR.hits.annotated.full$IsIsland, age9.SEBR.hits.annotated.full$IsSig) #chisq=17.6, pval=2.69E-5
# chiIS$observed
# chiIS$expected #sig probe count in islands is LOWER than expected
# 
# chiOS <- chisq.test(age9.SEBR.hits.annotated.full$IsOpenSea, age9.SEBR.hits.annotated.full$IsSig) #chisq=17.6, pval=2.69E-5
# chiOS$observed
# chiOS$expected #sig probe count in OS is HIGHER than expected
# 
# #Iterate through the 5 new formed variables and get chi square values
# factor_variables <- c("IsIsland", "IsOpenSea","IsNShore","IsNShelf","IsSShore","IsNShore")
# 
# chisqRESULTS <- lapply(age9.SEBR.hits.annotated.full[, factor_variables], function(x) 
#   chisq.test(age9.SEBR.hits.annotated.full$IsSig, x))
# chisqRESULTS
# 
# chisqRESULTSob <- lapply(age9.SEBR.hits.annotated.full[, factor_variables], function(x) 
#   chisq.test(age9.SEBR.hits.annotated.full$IsSig, x)$observed)
# chisqRESULTSob
# #try manual chi square test
# #chisqval = 17.71 (hand done)--very close to the above chi square rests of 17.626
# 
# 
# 
# 
# #make freq
# freqdata <- age9.SEBR.hits.annotated.full %>%
#   group_by(IsSig,island)%>%
#   summarise(n=n())%>%
#   mutate(freq=n/sum(n))
# 
# freqfull <- age9.SEBR.hits.annotated.full %>%
#   group_by(island)%>%
#   summarise(n=n())%>%
#   mutate(freq=n/sum(n))
# 
# #merge this data
# freqd9 <- merge(freqdata %>% filter(IsSig=="Yes"),freqfull, by="island") %>%
#   rename(n.x="n_sig",
#          freq.x="freq_sig",
#          n.y="n_tot",
#          freq.y="freq_tot") %>%
#   mutate(freqsigsites = n_sig/n_tot,
#          expcountSig = n_tot*(nrow(TopHits_SEBR_247)/nrow(TopHits_SEBR_full)),
#          dev=n_sig - expcountSig,
#          devperc = (dev/expcountSig)*100)
# 
# freqd9$n_tot-freqd9$n_sig
# #Determine if significant difference
# #table(freqd$)
# sum((freqd9$n_sig-freqd9$expcountSig)^2/freqd9$expcountSig) #2138.392
# pchisq(q=2138.392, df=5, lower.tail=F) # p value = 0 -- highly reject H0
# tabbo9 <- matrix(data=c(freqd9$n_sig,freqd9$n_tot-freqd9$n_sig),nrow=6, ncol=2)
# t9 <- chisq.test(tabbo9, correct=F) #X-squared = 2199, df = 5, p-value 0
# #why is this diffeerecent
# 
# #make bar plot of significant site distribution by CpG island region
# ggplot(data=freqdata, aes(x=island, y=freq, fill=IsSig))+
#   geom_bar(stat="identity", position="dodge")+ #dodge--puts bars side by side
#   theme_classic()+
#   theme(text= element_text(size=18))+
#   ylab("Frequency")+
#   xlab("CpG Island Region")+
#   scale_x_discrete(labels=c("Island" = "Island", "N_Shelf" = "North Shelf",
#                             "N_Shore" = "North Shore", "OpenSea"="Open Sea",
#                             "S_Shelf"="South Shelf","S_Shore"="South Shore"))+
#   guides(fill=guide_legend(title="Significant"))
# 
# #make bar plot of significant site distribution by reg_feature
# #NOTE: most sig hits are missing this entry
# #sum(age9.SEBR.hits.annotated$regulatory_feature != "") = 4217  (4217 of 12,998 have available data)
# #sum(age9.SEBR.hits.annotated$regulatory_feature == "") = 8781 (8781 of 12,998 hits are blank)
# 
# ggplot(data=age9.SEBR.hits.annotated %>% filter(regulatory_feature != ""), 
#        aes(x=regulatory_feature))+
#   geom_bar()+
#   theme_classic()+
#   theme(text= element_text(size=13),
#         axis.text.x = element_text(angle = 90))+
#   ylab("Count")+
#   xlab("Regulatory Feature")+
#   scale_x_discrete(labels=c("Gene_Associated" = "Gene associated", "NonGene_Associated" = "Non-gene associated",
#                             "Promoter_Associated" = "Promoter associated", "Unclassified"="Unclassified",
#                             "Gene_Associated_Cell_type_specific"="Gene associated CTS",
#                             "NonGene_Associated_Cell_type_specific"="Non-gene associated CTS",
#                             "Promoter_Associated_Cell_type_specific"="Promoter associated CTS",
#                             "Unclassified_Cell_type_specific"="Unclassified CTS"))
# 
# #####Are there any sig differences in proportions between sig and not sig sites for each genomic region
# 
# 
# 
# #AGe 15 repeat
# 
# # Annotate top hits
# #TopHits_modSEBR15_247 <- readRDS("/home/alreiner/Projects/ffcw/output/SigProbes/age15sigprobes_SEBR_794_9353.rds")
# age15.SEBR.hits.annotated <- add_annotations(TopHits_modSEBR15_247)
# write_csv(age15.SEBR.hits.annotated, file="/home/alreiner/Projects/ffcw/data/AnnotatedLmFitHits/age15_SEBR_247results_annot.csv")
# 
# 
# #make bar plot of significant site distribution by CpG island region
# ggplot(data=age15.SEBR.hits.annotated, aes(x=island))+
#   geom_bar()+
#   theme_classic()+
#   theme(text= element_text(size=18))+
#   ylab("Count")+
#   xlab("CpG Island Region")+
#   scale_x_discrete(labels=c("Island" = "Island", "N_Shelf" = "North Shelf",
#                             "N_Shore" = "North Shore", "OpenSea"="Open Sea",
#                             "S_Shelf"="South Shelf","S_Shore"="South Shore"))
# 
# # Load full hits SEBR age 15
# TopHits_modSEBR15_full <- readRDS("/home/alreiner/Projects/ffcw/output/SigProbes/age15sigprobes_SEBR_full794_CTPup.rds")
# 
# # Annotate top hits
# age15.SEBR.hits.annotated.full <- add_annotations(TopHits_modSEBR15_full)
# age15.SEBR.hits.annotated.full <- age15.SEBR.hits.annotated.full %>% 
#   mutate(IsSig = ifelse(rownames(age15.SEBR.hits.annotated.full) %in% rownames(TopHits_modSEBR15_247),"Yes","No"))
# 
# write_csv(age15.SEBR.hits.annotated.full, file="/home/alreiner/Projects/ffcw/data/AnnotatedLmFitHits/age15_SEBR_fullresults_annot_upCTP.csv")
# 
# #make freq
# freqdata15 <- age15.SEBR.hits.annotated.full %>%
#   group_by(IsSig,island)%>%
#   summarise(n=n())%>%
#   mutate(freq=n/sum(n))
# 
# freqfull15 <- age15.SEBR.hits.annotated.full %>%
#   group_by(island)%>%
#   summarise(n=n())%>%
#   mutate(freq=n/sum(n))
# 
# #merge this data15
# freqd <- merge(freqdata15 %>% filter(IsSig=="Yes"),freqfull15, by="island") %>%
#   rename(n.x="n_sig",
#          freq.x="freq_sig",
#          n.y="n_tot",
#          freq.y="freq_tot") %>%
#   mutate(freqsigsites = n_sig/n_tot,
#          expcountSig = n_tot*(nrow(TopHits_modSEBR15_247)/nrow(TopHits_modSEBR15_full)),
#          dev=n_sig - expcountSig,
#          devperc = (dev/expcountSig)*100)
# #Determine if significant difference
# #table(freqd$)
# sum((freqd$n_sig-freqd$expcountSig)^2/freqd$expcountSig) #2114.214
# 
# #[1] X2 = 2315.794 (close to automated results)
# pchisq(q=2114.214, df=5, lower.tail=F) # p value = 0 -- highly reject H0
# #cchisq.test(table) is not asking the same question
# tabbo <- matrix(data=c(freqd$n_sig,freqd$n_tot-freqd$n_sig),nrow=6, ncol=2)
# t <- chisq.test(tabbo, correct=F) #X-squared = 2164, df = 5, p-value 0
# #why is this diffeerecent
# tabbo
# 
# #At which island region does the significant difference exist?
# chisq.test(tabbo[c(1,2),], correct = FALSE)
# install.packages("devtools")
# devtools::install_github("ebbertd/chisq.posthoc.test")
# chisq.posthoc.test::chisq.posthoc.test(as.table(tabbo))
# 

