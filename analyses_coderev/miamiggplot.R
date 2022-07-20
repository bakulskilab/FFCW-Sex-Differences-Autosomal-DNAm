###################Load libraries
library(ewastools)
library(readr)
library(tidyverse)
library(limma)
library(ggplot2)

###################Load in appropriate data

TopHits_SEBR_full <-readRDS("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERfullblock.rds")
TopHits_modSEBR15_full <- readRDS("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age15SERfullblock.rds")
age15siggy <- readRDS("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age15SERsig247block.rds")
age9siggy <- readRDS("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERsig247block.rds")

mani <- ewastools:::manifest_450K
TopHits_sexmod1_full <-TopHits_SEBR_full #load in age 9 full results for ANY model
TopHits_mod015_full <- TopHits_modSEBR15_full
overlappo <- intersect(age9siggy$probe_id, age15siggy$probe_id) #8859 for SEBR model

#Turn into function
#input= tophitsresults9, tophitsresults15
sum(age15siggy$probe_id %in% age9siggy$probe_id) #
sum(!age15siggy$probe_id %in% age9siggy$probe_id) # unique to age 15
uniqueprobes15 <- age15siggy$probe_id[!age15siggy$probe_id %in% age9siggy$probe_id] #2145

sum(!age9sigprobes_12643 %in% sigprobes_10744_age15) 
uniqueprobes9 <- age9siggy$probe_id[!age9siggy$probe_id %in% age15siggy$probe_id] #4139

#Try Manhattan plot
#Get probe, chr, pos, p val into one data frame for age 9
top_age9 <- mani %>% filter(probe_id %in% rownames(TopHits_sexmod1_full))%>%select(probe_id,mapinfo,chr)

#merge age 9 probe p-vals to this dataframe
top_age9 <- merge(top_age9, TopHits_sexmod1_full %>% 
                    mutate(probe_id = rownames(TopHits_sexmod1_full))%>%
                    select(probe_id, P.Value), by="probe_id")
top_age9 <- top_age9 %>% dplyr::rename(SNP = probe_id,
                                       CHR = chr,
                                       POS = mapinfo,
                                       pvalue=P.Value)

#Get probe, chr, pos, p val into one data frame for age 15
bot_age15 <- mani %>% filter(probe_id %in% rownames(TopHits_mod015_full))%>%select(probe_id,mapinfo,chr)

#merge age 9 probe p-vals to this dataframe
bot_age15 <- merge(bot_age15, TopHits_mod015_full %>% 
                     mutate(probe_id = rownames(TopHits_mod015_full))%>%
                     select(probe_id,P.Value), by="probe_id")
bot_age15 <- bot_age15 %>% dplyr::rename(SNP = probe_id,
                                         CHR = chr,
                                         POS = mapinfo,
                                         pvalue=P.Value)

#Highlighting probes based on significance in only age 9, only age 15 or both
#overlap8651 = list of probes sig in both ages
#uniqueprobes9 = list of probes only sig at age 9
#uniqueprobes15 = list of probes only sig at age 15


#make variable for every tested age 9 probe that notes whether the probe was significant in ONLY age 9,
#only in age 15 or in both
top_age9 <- top_age9 %>% mutate(Status=ifelse(SNP %in% uniqueprobes9, "nine_only",
                                              ifelse(SNP %in% uniqueprobes15, "fifteen_only",
                                                     ifelse(SNP %in% overlappo, "both", "notsig"))))

top_age9 <- top_age9 %>% mutate(Coleur=ifelse(SNP %in% uniqueprobes9, "#990066",
                                              ifelse(SNP %in% uniqueprobes15, "#3399FF",
                                                     ifelse(SNP %in% overlappo, "#669900", "grey"))))

bot_age15 <- bot_age15 %>% mutate(Status=ifelse(SNP %in% uniqueprobes9, "nine_only",
                                                ifelse(SNP %in% uniqueprobes15, "fifteen_only",
                                                       ifelse(SNP %in% overlappo, "both", "notsig"))))
bot_age15 <- bot_age15 %>% mutate(Coleur=ifelse(SNP %in% uniqueprobes9, "#990066",
                                                ifelse(SNP %in% uniqueprobes15, "#3399FF",
                                                       ifelse(SNP %in% overlappo, "#669900", "grey"))))

#this is key for ordering the colors
#was as.numeric(CHR) for both top and bot
top_age9 <- top_age9 %>% group_by(CHR)%>%arrange(as.numeric(CHR),POS)
bot_age15 <- bot_age15 %>% group_by(CHR)%>%arrange(CHR,POS)

top_age9 <- top_age9[order(as.numeric(top_age9$CHR)),]
bot_age15 <- bot_age15[order(as.numeric(bot_age15$CHR)),]

saveRDS(top_age9, file="/home/alreiner/Projects/ffcw/data/top_age9_miami_SEBR.rds")
saveRDS(bot_age15, file="/home/alreiner/Projects/ffcw/data/bot_age15_miami_SEBR.rds")

topp09 <- readRDS("/home/alreiner/Projects/ffcw/data/top_age9_miami_SEBR.rds")
#topp09 <- topp09[,-1]
botto15 <- readRDS("/home/alreiner/Projects/ffcw/data/bot_age15_miami_SEBR.rds")
#botto15 <- botto15[,-1]
#topp09 <- top_age9
#botto15 <- bot_age15


#Make the data only the probes that are commonly tested at age 9 and 15
topp09 <- topp09 %>% filter(SNP %in% botto15$SNP)
botto15 <- botto15 %>% filter(SNP %in% topp09$SNP)

####################Try ggplot version of miami
# Prepare the dataset
class(topp09$CHR)
topp09$CHR <- factor(topp09$CHR,levels=c(1:22, "X"))
class(topp09$CHR)


age9try <- topp09 %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(topp09, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(as.numeric(CHR), POS) %>%
  mutate(BPcum=POS+tot) #%>%

# Add highlight and annotation information
#mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
#mutate( is_annotate=ifelse(-log10(P)>4, "yes", "no")) 

# Prepare X axis
class(age9try$CHR)
axisdf <- age9try %>% group_by(CHR)%>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
axisdf$CHR <- as.numeric(axisdf$CHR)
axisdf <- axisdf[order(as.numeric(axisdf$CHR)),]
loo <- age9try %>% group_by(CHR) %>% summarise(minno=min(BPcum),maxxo=max(BPcum))
loo$CHR <- as.numeric(loo$CHR)
loo <- loo[order(as.numeric(loo$CHR)),]

# Make the plot
#ggplot(age9try, aes(x=BPcum, y=-log10(pvalue))) +
ninePLOT <- ggplot(data=age9try) +
  # Show all points
  #geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  #scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  geom_point(aes(x=BPcum, y=-log10(pvalue),color=as.factor(Status)), alpha=0.5, size=1) +
  #scale_color_brewer(palette = "PuBuGn")+
  scale_color_manual(values = c("#990066","#3399FF","#66CC00","#99CCCC")) +
  #facet_grid(.~CHR)+
  geom_rect(data=loo, aes(xmin=minno,
                          xmax = maxxo, ymin=0, ymax=Inf,fill=factor(CHR)),alpha=0.3)+
  scale_fill_manual(values = rep(c("grey", "white"), 22 )) +
  # custom X axis:
  scale_x_continuous( label = as.factor(axisdf$CHR), breaks= axisdf$center,expand=expansion(mult=0,add=1e6)) +
  #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  ylim(0,300)+
  # Add highlighted points
  #geom_point(aes(color=as.factor(Status))) +
  
  # Add label using ggrepel to avoid overlapping
  #geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
  
  # Custom the theme:
  theme_bw() +
  ggtitle("Age 9")+
  xlab("")+
  theme( 
    text=element_text(size=16),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(hjust = 0.5,size=rel(1)),
    axis.title.x = element_text(hjust=0),
    axis.text.x = element_text(vjust=-2.2, size=10), #vjust=-2.2
    axis.text.y = element_text(size=20))+
    #axis.title.x = element_text(vjust=-0.5))
  #)+
  coord_cartesian(clip="off") #+  #in order to plot outside the plot box
  #annotation_custom(grid::textGrob("Chromosome", -0.06, -.14, 
   #                                default.units = "native", just="left",
    #                             gp=grid::gpar(fontsize = 15, col="black", lineheight=0.9, font=c(plain=1))))

#-0.04 too far right
#-0.15 too far down

#-0.06,-0.12

# Prepare the dataset
class(botto15$CHR)
botto15$CHR <- factor(botto15$CHR,levels=c(1:22, "X"))
class(botto15$CHR)
age15try <- botto15 %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(botto15, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(as.numeric(CHR), POS) %>%
  mutate(BPcum=POS+tot) #%>%


# Add highlight and annotation information
#mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
#mutate( is_annotate=ifelse(-log10(P)>4, "yes", "no")) 

# Prepare X axis
class(age15try$CHR)
axisdf15 <- age15try %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
axisdf15$CHR <- as.numeric(axisdf15$CHR)
axisdf15 <- axisdf15[order(as.numeric(axisdf15$CHR)),]
loo15 <- age15try %>% group_by(CHR)%>% summarise(minno=min(BPcum),maxxo=max(BPcum))
loo15$CHR <- as.factor(loo15$CHR)
loo15 <- loo15[order(as.numeric(loo15$CHR)),]

#age15try$Status <- factor(age15try$Status, levels=c("Both","Fifteen Only","Nine Only", "Not significant"))
#levels(age15try$Status)

# Make the plot
#ggplot(age9try, aes(x=BPcum, y=-log10(pvalue))) +
fitePLOT <- ggplot(data=age15try) +
  # Show all points
  #geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  #scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  geom_point(aes(x=BPcum, y=-log10(pvalue),color=as.factor(Status)), alpha=0.5, size=1) +
  #scale_color_brewer(palette = "PuBuGn")+
  scale_color_manual(values = c("#990066","#3399FF","#66CC00","#99CCCC"), labels=c("Both","15 only","9 only","Not significant")) +
  #facet_grid(.~CHR)+
  geom_rect(data=loo15, aes(xmin=minno,
                            xmax = maxxo, ymin=0, ymax=Inf,fill=factor(CHR)),alpha=0.3)+
  scale_fill_manual(values = rep(c("grey", "white"), 22 ), guide="none") +
  # custom X axis:
  scale_x_continuous( label = axisdf15$CHR, breaks= axisdf15$center, 
                      position="top", expand=expansion(mult=0,add=1e6)) +
  #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  #ylim(0,350)+
  # Add highlighted points
  #geom_point(aes(color=as.factor(Status))) +
  
  # Add label using ggrepel to avoid overlapping
  #geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    text=element_text(size=16),
    legend.position="bottom",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=20)#vjust element text to even out axis labels
    
  )+
  scale_y_reverse(limits=c(300,0))+
  labs(caption="Age 15")+ 
  theme(plot.caption = element_text(hjust=0.5, size=rel(1)))+
  guides(color=guide_legend(title="Significance Status"))
  #coord_cartesian(clip="off") + 
  #annotation_custom(grid::textGrob("Hello Down Here", -0.05, -.09, default.units = "native", just="left",
  #                                 gp=grid::gpar(fontsize = 11, col="black", lineheight=0.9, font=c(plain=1))))
  
  
  
fitePLOT
ninePLOT

library(grid)
grid.newpage()
#png(paste0(imagesdir,"/PaperManuscript/miami_lr.png"),width=9,height=7,units = "in",res=150)
png(paste0(imagesdir,"/PaperManuscript/miami_hr.png"),width=9,height=7,units = "in",res=800)
grid.draw(rbind(ggplotGrob(ninePLOT), ggplotGrob(fitePLOT), size = "last"))
dev.off()
#miamitry2_mf = very close--Chr slightly over
#fitePLOT

