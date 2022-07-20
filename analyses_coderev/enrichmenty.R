#load libraries
#remotes::install_github("coolbutuseless/ggpattern")
#library(ggpattern)
library(readr)
library(tidyverse)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(missMethyl)
imagesdir <- paste0(here::here(),'/images')


#load in full age 9 results and age 15 results
age9sigSEBR <-readRDS("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERsig247block.rds")

age9fullSEBR <- readRDS("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERfullblock.rds")

#age15fullSEBR <- readRDS("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age15SERfullblock.rds")

View(Islands.UCSC$Relation_to_Island)
#create add_annotations function
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

#see <- add_annotations(age9fullSEBR)
#add annotation to age 9 full data
#saveRDS(see, paste0("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9_SEBR_full_annot",nrow(age9fullSEBR),".rds"))

#add annotation to age 9 full data
#saveRDS(add_annotations(age15fullSEBR), paste0("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age15_SEBR_full_annot",nrow(age15fullSEBR),".rds"))


age9.SEBR.hits.annotated.full <- readRDS(paste0("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9_SEBR_full_annot",nrow(age9fullSEBR),".rds"))

age9.SEBR.hits.annotated.full <- age9.SEBR.hits.annotated.full %>% 
  mutate(IsSig = ifelse(rownames(age9.SEBR.hits.annotated.full) %in% rownames(age9sigSEBR),"Yes","No"))

sum(is.na(age9.SEBR.hits.annotated.full$island)) #0 sites missing location data as hoped

#get names of probes that did not have island data at age 15 weirdly
noisland <- readRDS("/home/alreiner/Projects/ffcw/output/noisland15.rds")

#create list of data for probes with missing data at age 15
noinfo <- age9.SEBR.hits.annotated.full %>% filter(probe_id %in% noisland)
saveRDS(noinfo, file="/home/alreiner/Projects/ffcw/output/noinfo.rds")

#Create overrep plot (11/5/21)
#make freq
freqdata <- age9.SEBR.hits.annotated.full %>%
  group_by(IsSig,island)%>%
  summarise(n=n())%>%
  mutate(freq=n/sum(n))

freqfull <- age9.SEBR.hits.annotated.full %>%
  group_by(island)%>%
  summarise(n=n())%>%
  mutate(freq=n/sum(n))

freqd9 <- merge(freqdata %>% filter(IsSig=="Yes"),freqfull, by="island") 
names(freqd9)
freqd9 <- freqd9%>%
  dplyr::rename(n_sig="n.x",
         freq_sig="freq.x",
         n_tot="n.y",
         freq_tot="freq.y")%>%
  pivot_longer(cols=c("freq_sig","freq_tot"),names_to = "freq")

#chisq test for homogeneity
#answers the question of: is the distribution of significant sites across island regions
# the same as the distribution of ALL sites across island regions
freqd9_test <- merge(freqdata %>% filter(IsSig=="Yes"),freqfull, by="island") %>%
  dplyr::rename("n_sig"="n.x",
         freq_sig="freq.x",
         n_tot="n.y",
         freq_tot="freq.y")

sigVall <- rbind(freqd9_test$n_sig,freqd9_test$n_tot)
colnames(sigVall) <- c(freqd9_test$island)
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
n_sig <- freqd9_test$n_sig #total num probes sig per region
exp_sig <- freqd9_test$freq_tot #total freq per region
totSig9 <- nrow(age9sigSEBR)
k <- c()
j<-c()
for (i in 1:length(n_sig)){
  k=prop.test(x = n_sig[i], n = totSig9, p = exp_sig[i], 
              correct = FALSE)
  print(k$p.value)
}

#Island: 2.443419e-05
# NShelf: 8.390827e-17
# NShore: 1.448741e-106
# OS 4.316531e-191
# [SShelf] 2.420602e-08
# [SShore] 1.300244e-125

#check island results manually 
phat=3090/8430
pnot=135117/391980
sampsiz = 8430
zstat=(phat-pnot)/(sqrt((pnot*(1-pnot))/sampsiz))
2*pnorm(zstat, mean = 0, sd = 1, lower.tail = F) 
#2-sided p-value = 2.44E-5 (which matches prop.test results above)

#change ordering of factor levels!
freqd9$island <- factor(freqd9$island, levels = c("N_Shelf","N_Shore","Island","S_Shore","S_Shelf","OpenSea"))

#make bar plot of significant site distribution by CpG island region
age9dist <- ggplot(data=freqd9, aes(x=island, y=value, fill=freq, label=scales::percent(value %>% round(3),accuracy=0.1)))+
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
            size = 5)
#+guides(fill=guide_legend(title="Proportion of"))
age9dist
png(paste0(imagesdir,"/enrich/enrichyplot_lr.png"),width=11,height=7,units = "in",res=150)
png(paste0(imagesdir,"/enrich/enrichyplot.png"),width=11,height=7,units = "in",res=800)
age9dist
dev.off()

png(paste0(imagesdir,"/enrich/enrichyplot_lr.png"),width=9,height=7,units = "in",res=150)
age9dist
dev.off()

# #test if sig difference between each blue and grey bar (observed vs expected props)
# prop.test(x = c(3529, 136239), n = c(11136, 406114)) #p=4.39E-05
# prop.test(x = c(381, 19311), n = c(11136, 406114)) #p<6.78E-11
# prop.test(x = c(2538, 53749), n = c(11136, 406114)) #p<2.2E-16
# prop.test(x = c(2295, 137767), n = c(11136, 406114)) #p<2.2E-16
# prop.test(x = c(345, 17122), n = c(11136, 406114)) #p=7.13E-9
# prop.test(x = c(2048, 41926), n = c(11136, 406114)) #p<2.2E-16
# 
