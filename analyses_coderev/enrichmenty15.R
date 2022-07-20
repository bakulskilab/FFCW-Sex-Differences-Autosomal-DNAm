#load libraries
#remotes::install_github("coolbutuseless/ggpattern")
#library(ggpattern)
library(readr)
library(tidyverse)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(missMethyl)
imagesdir <- paste0(here::here(),'/images')


#load in full age 9 results and age 15 results
age15sigSEBR <-readRDS("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age15SERsig247block.rds")

age15fullSEBR <- readRDS("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age15SERfullblock.rds")

age15.SEBR.hits.annotated.full <- readRDS(paste0("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age15_SEBR_full_annot",nrow(age15fullSEBR),".rds"))



add_annotations<-function(data){
  #add annotations
  Locations <- Locations[rownames(data), ]
  Islands.UCSC <- Islands.UCSC[rownames(data),]
  Other <- Other[rownames(data),]
  # add columns of interest from annotation to annotated version
  data.annotated <- data
  data.annotated$chr <- Locations$chr
  data.annotated$pos <- Locations$pos
  data.annotated$island <- Islands.UCSC$Relation_to_Island
  data.annotated$gene <- Other$UCSC_RefGene_Name
  data.annotated$gene_group <- Other$UCSC_RefGene_Group
  data.annotated$regulatory_feature <- Other$Regulatory_Feature_Group
  return(data.annotated)
}

#add annotation to age 9 full data
saveRDS(add_annotations(age15fullSEBR), paste0("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age15_SEBR_full_annot_2",nrow(age15fullSEBR),".rds"))


age15.SEBR.hits.annotated.full <- age15.SEBR.hits.annotated.full %>% 
  mutate(IsSig = ifelse(age15.SEBR.hits.annotated.full$probe_id %in% age15sigSEBR$probe_id,"Yes","No"))

sum(is.na(age15.SEBR.hits.annotated.full$island)) #27,180 sites missing location data weirdly

#check if that site is in age 9 data
testannot_nodat <- age15.SEBR.hits.annotated.full[which(is.na(age15.SEBR.hits.annotated.full$island)),]

#probes with no data for CpG island region
noisland <- c(testannot_nodat$probe_id)
saveRDS(noisland, file="/home/alreiner/Projects/ffcw/output/noisland15.rds")

#island information taken from the age 9 annotated data which actually has info for the island regions
noinfo <- readRDS("/home/alreiner/Projects/ffcw/output/noinfo.rds")

age15.SEBR.hits.annotated.full<- left_join(age15.SEBR.hits.annotated.full, noinfo %>% select(probe_id,island), by="probe_id")
age15.SEBR.hits.annotated.full <- age15.SEBR.hits.annotated.full %>% mutate(islandtest=ifelse(is.na(island.x),island.y,island.x))

#Create overrepresentation plot
#make freq
freqdata <- age15.SEBR.hits.annotated.full %>%
  group_by(IsSig,islandtest)%>%
  summarise(n=n())%>%
  mutate(freq=n/sum(n))

freqfull <- age15.SEBR.hits.annotated.full %>%
  group_by(islandtest)%>%
  summarise(n=n())%>%
  mutate(freq=n/sum(n))

freqd15 <- merge(freqdata %>% filter(IsSig=="Yes"),freqfull, by="islandtest") 

names(freqd15)
freqd15 <- freqd15%>%
  dplyr::rename(n_sig="n.x",
                freq_sig="freq.x",
                n_tot="n.y",
                freq_tot="freq.y")%>%
  pivot_longer(cols=c("freq_sig","freq_tot"),names_to = "freq")

#chisq test for homogeneity
#answers the question of: is the distribution of significant sites across island regions
# the same as the distribution of ALL sites across island regions
freqd15_test <- merge(freqdata %>% filter(IsSig=="Yes"),freqfull, by="islandtest") %>%
  dplyr::rename("n_sig"="n.x",
                freq_sig="freq.x",
                n_tot="n.y",
                freq_tot="freq.y")

sigVall <- rbind(freqd15_test$n_sig,freqd15_test$n_tot)
colnames(sigVall) <- c(freqd15_test$island)
rownames(sigVall) <- c("Sig","All")
sigVall
showtbl <- cbind(Sites=c("Sig","All"),sigVall)
showtbl
gt(as.data.frame(showtbl))
gt(as.data.frame(sigVall))
chisq.test(sigVall) #p<2.2E-16 There is significant difference in distributions of significant sites
# across island regions than overall sites

#Where do the significant differences exist? 
#1 proportion z-test (compare proportion in each island region to a fixed
# proportion )
n_sig <- freqd15_test$n_sig #total num probes sig per region
exp_sig <- freqd15_test$freq_tot #total freq per region
totSig15 <- nrow(age15sigSEBR)
k <- c()
j<-c()
for (i in 1:length(n_sig)){
  k=prop.test(x = n_sig[i], n = totSig15, p = exp_sig[i], 
              correct = FALSE)
  print(k$p.value)
}

#Island: 0.0001213494
# NShelf:  9.798583e-17
# NShore: 8.085702e-114
# OS 9.592863e-194
# [SShelf] 6.166095e-10
# [SShore] 1.332154e-131


#change ordering of factor levels!
freqd15$islandtest <- factor(freqd15$islandtest, levels = c("N_Shelf","N_Shore","Island","S_Shore","S_Shelf","OpenSea"))

#make bar plot of significant site distribution by CpG island region
age15dist <- ggplot(data=freqd15, aes(x=islandtest, y=value, fill=freq, label=scales::percent(value %>% round(3),accuracy=0.1)))+
   geom_bar(stat="identity", position="dodge")+ #dodge--puts bars side by side
   theme_classic()+
   theme(text= element_text(size=15),
         axis.text.y= element_text(size=18),
         axis.text.x= element_text(size=18),
         axis.title = element_text(size=22))+
   ylab("Proportion of DNA Methylation Sites")+
   xlab("CpG Island Region")+
   #scale_fill_discrete(name="Proportion of",labels=c("Significant sites", "Total sites"))+
   scale_fill_manual(values=c("black","grey"),name="Proportion of DNA \n Methylation Sites",labels=c("Significantly different by sex", "All sites"))+
   scale_x_discrete(labels=c("Island" = "Island", "N_Shelf" = "North\nShelf",
                            "N_Shore" = "North\nShore", "OpenSea"="Open\nSea",
                            "S_Shelf"="South\nShelf","S_Shore"="South\nShore")) +
   geom_text(position = position_dodge(width = .9),    # move to center of bars
             vjust = -0.5,    # nudge above top of bar
             size = 5)+guides(fill=guide_legend(title="Proportion of"))
age15dist
# png(paste0(imagesdir,"/enrich/enrichyplot_lr.png"),width=11,height=7,units = "in",res=150)
# png(paste0(imagesdir,"/enrich/enrichyplot.png"),width=11,height=7,units = "in",res=800)
# age9dist
# dev.off()
#
# png(paste0(imagesdir,"/enrich/enrichyplot_lr.png"),width=9,height=7,units = "in",res=150)
# age9dist
# dev.off()