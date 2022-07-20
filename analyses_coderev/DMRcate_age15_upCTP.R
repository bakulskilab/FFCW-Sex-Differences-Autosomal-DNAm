library(readr)
library(limma)
library(tidyverse)
library(gt)
library(gtsummary)
library(missMethyl)

#load analytic sample with PC data
analysis2dataPC <- readRDS("/home/alreiner/Projects/ffcw/data/age15final_796.rds")
#analysis2dataPC <-analysis2dataPC[,-c(1,2)]

#get beta matrix without sex chroms or gap probes
#load og beta

#beta <- readRDS("/nfs/turbo/bakulski1/People/alreiner/data/betaqc.rds")

# ########sex chrom removal
# library(ewastools)
# library(data.table)
# mani <- ewastools:::manifest_450K

#data(IlluminaHumanMethylation450kmanifest)

# # get a list of cpg sites that are on sex chromosomes
# sexprobes <- mani %>%
#   filter(chr %in% c("X","Y"))%>%
#   select(probe_id)
# # CHANGE GAPS
# gapprobes9 <- readRDS("/home/alreiner/Projects/ffcw/output/gapprobes9.rds")
# 
# # ####################remove sex and gap probes
# beta_nogap9 <- beta[!rownames(beta) %in% gapprobes9,match(analysis2dataPC$MethID,colnames(beta))] #left with 406,148 probes
# beta_nogap_noXY9_fin <- beta_nogap9[!rownames(beta_nogap9) %in% sexprobes$probe_id,match(analysis2dataPC$MethID,colnames(beta_nogap9))] #left with 406,148 probes
# #only 405,527 probes with the 15 year old gapes removed and sex probes
# #saveRDS(beta_nogap_noXY9_fin, "/home/alreiner/Projects/ffcw/data/age9beta_nogapnoXY_816.rds")


betanogapnoXY_15 <- readRDS("/home/alreiner/Projects/ffcw/data/beta_nogapnoXY_15_796.rds")
betanogapnoXY_15 <- betanogapnoXY_15[,match(analysis2dataPC$MethID,colnames(betanogapnoXY_15))]
#get the SEBR significant probes
#SEBRsigprobes <- readRDS("/home/alreiner/Projects/ffcw/output/SigProbes/age9sigprobes_SEBR_816_11136.rds")

#NEED TO CHANGE IF LMFIT RESULTS CHANGE##############
#line 607 in results_writeup.Rmd
SEBRsigprobes <- readRDS("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age15SERsig247block.rds")
SEBRsigprobes <- rownames(SEBRsigprobes)
######################

library(DMRcate)
design <- model.matrix(~factor(analysis2dataPC$cm1bsex)+
                         analysis2dataPC$Epithelial.cells+factor(analysis2dataPC$cm1ethrace))

dupcor <- duplicateCorrelation(betanogapnoXY_15,design=design,block=factor(analysis2dataPC$Sample_Plate))

myannotation <- cpg.annotate("array", betanogapnoXY_15[rownames(betanogapnoXY_15) %in% SEBRsigprobes,], 
                             arraytype = "450K",annotation=c(array = "IlluminaHumanMethylation450k", 
                                                             annotation = "ilmn12.hg19"),
                             analysis.type="differential", 
                             what="Beta",
                             design=design, 
                             block=factor(analysis2dataPC$Sample_Plate), 
                             correlation=dupcor$consensus.correlation,
                             contrasts=F, 
                             cont.matrix = NULL, 
                             coef=2, 
                             fdr=1) 


myannotationFDR0.05 <- cpg.annotate("array", betanogapnoXY_15, 
                                    arraytype = "450K",annotation=c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"),
                                    analysis.type="differential", 
                                    what="Beta",
                                    design=design, 
                                    block=factor(analysis2dataPC$Sample_Plate), 
                                    correlation=dupcor$consensus.correlation,
                                    contrasts=F, 
                                    cont.matrix = NULL, 
                                    coef=2, 
                                    fdr=0.05) 

#myannotation2 <- cpg.annotate("array", beta_nogap_noXYmani, arraytype = "450K",
#                              analysis.type="differential", design=design, coef=2) 
#returned 45,332 sig probes (make sure to reorddr betanogapyXymani to match analysis2data order!)

#both annotation 1 and 2 give same number of sig diff probes

dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
dmrcoutput2 <- dmrcate(myannotationFDR0.05, lambda=1000, C=2)

#dmrcoutput_247 <- dmrcate(myannotation, lambda=1000, C=2, p.adjust.method="none")

#fix error where can't connect to experimentHub 
library(ExperimentHub)
setExperimentHubOption(arg="URL", value = 'https://experimenthub.bioconductor.org')
#"PROXY", httr::use_proxy(), 

results.ranges_15 <- extractRanges(dmrcoutput, genome = "hg19")
results.ranges_15_FDR0.05 <- extractRanges(dmrcoutput2, genome = "hg19")
#2304 DMRs for sex,epi, batch model (get same count when using what=default and what="Beta")
#2362 DMRs for sex,epi,batch,mom race

#DMR results with overlapping genes
saveRDS(as.data.frame(results.ranges_15), "/home/alreiner/Projects/ffcw/output/DMR/DMR_results_247_SEBR_15.rds")
saveRDS(as.data.frame(results.ranges_15_FDR0.05), "/home/alreiner/Projects/ffcw/output/DMR/DMR_results_247_SEBR_15_fdr05.rds")



# results.ranges_15 <- read_csv("/home/alreiner/Projects/ffcw/output/DMR/DMR_results_247_SEBR_upCTP15.csv")
# results.ranges_15 <- results.ranges_15[,-1]
results.ranges_15G <- makeGRangesFromDataFrame(results.ranges_15)


#Picking DMR sig threshold
#View(as.data.frame(results.ranges))
#View(results.ranges$meandiff)
#View(as.data.frame(results.ranges_09)%>% filter(no.cpgs > 3,abs(meandiff)>0.09))

#using sig thresh of > 3 CpG and abs(meandiff)>0.04
# sigDMRdata_SEBR_9 <- as.data.frame(results.ranges_09) %>% filter(no.cpgs > 3, abs(meandiff)>0.04) #111
# TOPsigDMR_9 <- sigDMRdata_SEBR_9 %>% filter(abs(meandiff)>0.07) %>% select(seqnames,no.cpgs,Fisher,maxdiff,meandiff,overlapping.genes)
# genes9 <- sigDMRdata_SEBR_9$overlapping.genes
# #path to read updated ctp SEBR DMR results
# #DMR_results_816_247_SEBR_upCTP <- read_csv("output/DMR/DMR_results_816_247_SEBR_upCTP.csv")
# 
# #Print out final DMR list for results
# rounded9 <- TOPsigDMR_9 %>% 
#   mutate_at(vars(maxdiff,meandiff), funs(round(., 3)))
# #save txt file of results to input into word
# write.table(rounded9, file="/home/alreiner/Projects/ffcw/output/DMR/age15_DMRresults_upCTP_1_11_22.txt",sep=',')
# 

#Calculate concordance of overlapping genes

#separate double genes into 2 sep elemts in a list
# genes_9_split <- unlist(strsplit(genes9, "[,]"))
# print(genes_9_split)
# #get rid of NAs
# genessplit_noNA <- na.omit(unlist(strsplit(genes9, "[,]")))
# print(genessplit_noNA)
# #trim whitespace
# #genessplit_noNA <- na.omit(unlist(strsplit(genes9, "[,]")))
# genessplit_noNA <- trimws(genessplit_noNA) #112
# print(genessplit_noNA)
# #delete the double genes
# sum(duplicated(genessplit_noNA)) #3
# which(duplicated(genessplit_noNA))
# genessplit_noNA <- genessplit_noNA[!duplicated(genessplit_noNA)] #final list of genes for age 9
# saveRDS(genessplit_noNA, file=paste0("/home/alreiner/Projects/ffcw/output/DMR/DMR_overlappgenelist_SEBR_15_",format(Sys.time(), "%b%d"),".rds"))


#GO term enrichment
#get GO ontology terms associated with the DMRs
enrichment_GO_15 <- goregion(results.ranges_15G, all.cpg = rownames(beta_nogap_noXYmani),
                               collection = "GO", array.type = "450K",sig.genes = T)
enrichment_GO_15 <- enrichment_GO_15[order(enrichment_GO_15$P.DE),]
readRDS(enrichment_GO_15, "/home/alreiner/Projects/ffcw/output/DMR/age15_GOresults_fdr1.rds")

#Filter GO terms to biological processes only
enrichment_GO_15_bp<- enrichment_GO_15[order(enrichment_GO_15$P.DE),]%>% 
  filter(ONTOLOGY=="BP")
#Further filter to FDR < 0.05
enrichment_GO_15sig_filtFDR <- enrichment_GO_15_bp %>%
  filter(FDR<0.05) %>%
  arrange(FDR)%>%
  mutate_at(vars(P.DE, FDR), funs(round(., 3)))
#write.table(enrichment_GO_1717sig_15_filt, file="/home/alreiner/Projects/ffcw/analyses/15yo/DMR_1717_GoTerms_BP_11_24.txt",sep=',')
readRDS(enrichment_GO_15sig_filtFDR, "/home/alreiner/Projects/ffcw/output/DMR/age15GO_results_filtFDR05.rds")
