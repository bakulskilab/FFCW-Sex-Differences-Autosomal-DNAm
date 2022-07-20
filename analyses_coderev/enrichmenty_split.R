#load libraries
library(readr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(missMethyl)
library(tidyverse)
imagesdir <- paste0(here::here(),'/images')


#load in full age 9 results and age 15 results
age9sigSEBR <-readRDS("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERsig247block.rds")

#subset to malebiased probes

#subset to female biased probes
sigprobes_f_SEBR9 <-age9sigSEBR %>% filter(logFC>0)%>%rownames() # female biased (15)
sigprobes_m_SEBR9 <-age9sigSEBR %>% filter(logFC<0)%>%rownames() # male biased (15)


age9fullSEBR <- readRDS("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERfullblock.rds")

age15fullSEBR <- readRDS("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age15SERfullblock.rds")

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


age9.SEBR.hits.annotated.full <- readRDS(paste0("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9_SEBR_full_annot",nrow(age9fullSEBR),".rds"))

sum(is.na(age9.SEBR.hits.annotated.full$island))


age9.SEBR.hits.annotated.full <- age9.SEBR.hits.annotated.full %>% 
  mutate(IsSig = ifelse(rownames(age9.SEBR.hits.annotated.full) %in% rownames(age9sigSEBR),"Yes","No"),
         SexBias=ifelse(rownames(age9.SEBR.hits.annotated.full) %in% sigprobes_f_SEBR9,"Fembiased",
                        ifelse(rownames(age9.SEBR.hits.annotated.full) %in% sigprobes_m_SEBR9,"Malebiased", "neither")))

table(age9.SEBR.hits.annotated.full$SexBias)

#Create overrep plot for female  biased probes
#make freq
freqdataf <- age9.SEBR.hits.annotated.full %>%   #just need to filter out the  male biased probes
  filter(SexBias != "Malebiased")%>%
  group_by(IsSig,island)%>%
  summarise(n=n())%>%
  mutate(freq=n/sum(n))


freqfullf <- age9.SEBR.hits.annotated.full %>%
  #filter(SexBias != "Malebiased")%>% #want all sites
  group_by(island)%>%
  summarise(n=n())%>%
  mutate(freq=n/sum(n))

freqd9f<- merge(freqdataf %>% filter(IsSig=="Yes"),freqfullf, by="island")%>%
  dplyr::rename("n_sig"="n.x",
         "freq_sig"="freq.x",
         "n_tot"="n.y",
         "freq_tot"="freq.y")%>%
  pivot_longer(cols=c("freq_sig","freq_tot"),names_to = "freq")

#chisq test for homogeneity
#answers the question of: is the distribution of significant sites across island regions
# the same as the distribution of ALL sites across island regions
freqd9_test <- merge(freqdata %>% filter(IsSig=="Yes"),freqfull, by="island") %>%
  rename(n.x="n_sig",
         freq.x="freq_sig",
         n.y="n_tot",
         freq.y="freq_tot")

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

freqd9f_test<- merge(freqdataf %>% filter(IsSig=="Yes"),freqfullf, by="island")%>%
  dplyr::rename("n_sig"="n.x",
                "freq_sig"="freq.x",
                "n_tot"="n.y",
                "freq_tot"="freq.y")


#1 proportion z-test (compare proportion in each island region to a fixed
# proportion )
n_sig <- freqd9f_test$n_sig
exp_sig <- freqd9f_test$freq_tot
freqdataf_filt <- freqdataf %>% filter(IsSig=="Yes")
totSig9 <- sum(freqdataf_filt$n)
k <- c()
j<-c()
for (i in 1:length(n_sig)){
  k=prop.test(x = n_sig[i], n = totSig9, p = exp_sig[i], 
              correct = FALSE)
  print(k$p.value)
}
print(k)

#hand check results of Open Sea regions for hypermethylated sites in females
phat_f=1039/6425
pnot_f=129569/391980
sampsiz = 6425
zstat=(phat_f-pnot_f)/(sqrt((pnot_f*(1-pnot_f))/sampsiz))
2*pnorm(zstat, mean = 0, sd = 1, lower.tail = T) 
#1 sided p=2.60 E -182
#2 sided p = 5.20 E -182
prop.test(x = 1039, n = 6425, p = pnot_f, 
          correct = F) #p< 2.2 E -16
#X-squared=827.67 df=1


#Test North Shelf: 1.576323e-46
phat_f=174/6425
pnot_f=18225/391980
sampsiz = 6425
zstat=(phat_f-pnot_f)/(sqrt((pnot_f*(1-pnot_f))/sampsiz))
2*pnorm(zstat, mean = 0, sd = 1, lower.tail = T) 
#1-sided p-value=p=7.32E-14
#2-sided p-value = 1.46 E -13 (which matches prop.test results)

prop.test(x = 174, n = 6425, p = pnot_f, 
          correct = T) #p=1.46 E -13
#X-squared=827.67 df=1

#change ordering of factor levels!
freqd9$island <- factor(freqd9$island, levels = c("N_Shelf","N_Shore","Island","S_Shore","S_Shelf","OpenSea"))

#make bar plot of significant site distribution by CpG island region
age9dist <- ggplot(data=freqd9f, aes(x=island, y=value, fill=freq))+
  geom_bar(stat="identity", position="dodge")+ #dodge--puts bars side by side
  theme_classic()+
  theme(text= element_text(size=15))+#, legend.position = "none"
  ylab("Frequency")+
  xlab("CpG Island Region")+
  #scale_fill_discrete(name="Proportion of",labels=c("Significant sites", "Total sites"))+
  scale_fill_manual(values=c("navy","grey"),name="Proportion of",labels=c("Significant sites", "Total sites"))+
  scale_x_discrete(labels=c("Island" = "Island", "N_Shelf" = "North Shelf",
                            "N_Shore" = "North Shore", "OpenSea"="Open Sea",
                            "S_Shelf"="South Shelf","S_Shore"="South Shore"))+
  ggtitle("Age 9 Female-Biased Enrichment")
#+guides(fill=guide_legend(title="Proportion of"))

png(paste0(imagesdir,"/enrich/enrichyplot.png"),width=9,height=7,units = "in",res=800)
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



########################################
#Male-Biased
########################################

#Create overrep plot for male  biased probes
#make freq
freqdatam <- age9.SEBR.hits.annotated.full %>%   #just need to filter out the  male biased probes
  filter(SexBias != "Fembiased")%>%
  group_by(IsSig,island)%>%
  summarise(n=n())%>%
  mutate(freq=n/sum(n))


freqfullm <- age9.SEBR.hits.annotated.full %>%
  #filter(SexBias != "Fembiased")%>%
  group_by(island)%>%
  summarise(n=n())%>%
  mutate(freq=n/sum(n))

freqd9m<- merge(freqdatam %>% filter(IsSig=="Yes"),freqfullm, by="island")%>%
  rename("n_sig"="n.x",
         "freq_sig"="freq.x",
         "n_tot"="n.y",
         "freq_tot"="freq.y")%>%
  pivot_longer(cols=c("freq_sig","freq_tot"),names_to = "freq")

#change ordering of factor levels!
freqd9m$island <- factor(freqd9$island, levels = c("N_Shelf","N_Shore","Island","S_Shore","S_Shelf","OpenSea"))

#make bar plot of significant site distribution by CpG island region
age9dist_malebias <- ggplot(data=freqd9m, aes(x=island, y=value, fill=freq))+
  geom_bar(stat="identity", position="dodge")+ #dodge--puts bars side by side
  theme_classic()+
  theme(text= element_text(size=15))+#, legend.position = "none"
  ylab("Frequency")+
  xlab("CpG Island Region")+
  #scale_fill_discrete(name="Proportion of",labels=c("Significant sites", "Total sites"))+
  scale_fill_manual(values=c("navy","grey"),name="Proportion of",labels=c("Significant sites", "Total sites"))+
  scale_x_discrete(labels=c("Island" = "Island", "N_Shelf" = "North Shelf",
                            "N_Shore" = "North Shore", "OpenSea"="Open Sea",
                            "S_Shelf"="South Shelf","S_Shore"="South Shore"))+
  ggtitle("Age 9 Male-Biased Enrichment")
#+guides(fill=guide_legend(title="Proportion of"))

png(paste0(imagesdir,"/enrich/enrichyplot.png"),width=9,height=7,units = "in",res=800)
age9dist_malebias
dev.off()

#calc p values
freqd9m_test<- merge(freqdatam %>% filter(IsSig=="Yes"),freqfullm, by="island")%>%
  dplyr::rename("n_sig"="n.x",
         "freq_sig"="freq.x",
         "n_tot"="n.y",
         "freq_tot"="freq.y")

#1 proportion z-test (compare proportion in each island region to a fixed
# proportion )
n_sig <- freqd9m_test$n_sig
exp_sig <- freqd9m_test$freq_tot
freqdatam_filt <- freqdatam %>% filter(IsSig=="Yes")
totSig9 <- sum(freqdatam_filt$n)
k <- c()
j<-c()
for (i in 1:length(n_sig)){
  k=prop.test(x = n_sig[i], n = totSig9, p = exp_sig[i], 
              correct = FALSE)
  print(k$p.value)
}
print(k)
