library(readr)
library(limma)
library(tidyverse)
library(gt)
library(gtsummary)
library(missMethyl)

#load analytic sample with PC data
analysis2dataPC <- readRDS("/home/alreiner/Projects/ffcw/data/age9final_816.rds")

#age 15 data
analysisdat15 <- readRDS("/home/alreiner/Projects/ffcw/data/age15final_796.rds")
#originally: readRDS("age15fin767_jan2022.rds)

analysis2dataPC <- analysis2dataPC %>% filter(id %in% c(analysisdat15$idnum))

# #betas without gap or sex formed in: beta_formation.R
# #beta age 9, no gaps, no XY (n=796)

beta_nogap_noXYmani <- readRDS("/home/alreiner/Projects/ffcw/data/age9beta_nogapnoXY_796.rds")
beta_nogap_noXYmani <- beta_nogap_noXYmani[,match(analysis2dataPC$MethID,colnames(beta_nogap_noXYmani))]

#NEED TO CHANGE IF LMFIT RESULTS CHANGE######
#line 614 in results write_up
SEBRsigprobes <-readRDS("/home/alreiner/Projects/ffcw/output/BlockProbeResults/age9SERsig247block.rds")

##########
library(DMRcate)
design <- model.matrix(~factor(analysis2dataPC$cm1bsex)+
                         analysis2dataPC$Epithelial.cells+factor(analysis2dataPC$cm1ethrace)) 

dupcor <- duplicateCorrelation(beta_nogap_noXYmani,design=design,block=factor(analysis2dataPC$Sample_Plate))

myannotation <- cpg.annotate("array", beta_nogap_noXYmani[rownames(beta_nogap_noXYmani) %in% rownames(SEBRsigprobes),], 
                             arraytype = "450K",annotation=c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"),
                             analysis.type="differential", 
                             what="Beta",
                             design=design, 
                             block=factor(analysis2dataPC$Sample_Plate), 
                             correlation=dupcor$consensus.correlation,
                             contrasts=F, 
                             cont.matrix = NULL, 
                             coef=2, 
                             fdr=1) 

# myannotationFDR0.05 <- cpg.annotate("array", beta_nogap_noXYmani, 
#                              arraytype = "450K",annotation=c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"),
#                              analysis.type="differential", 
#                              what="Beta",
#                              design=design, 
#                              block=factor(analysis2dataPC$Sample_Plate), 
#                              correlation=dupcor$consensus.correlation,
#                              contrasts=F, 
#                              cont.matrix = NULL, 
#                              coef=2, 
#                              fdr=0.05) 


#myannotation2 <- cpg.annotate("array", beta_nogap_noXYmani, arraytype = "450K",
#                              analysis.type="differential", design=design, coef=2) 
#returned 45,332 sig probes (make sure to reorddr betanogapyXymani to match analysis2data order!)

#both annotation 1 and 2 give same number of sig diff probes

dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
#dmrcoutput0.05FDR <- dmrcate(myannotationFDR0.05, lambda=1000, C=2)
#dmrcoutput_247 <- dmrcate(myannotation, lambda=1000, C=2, p.adjust.method="none")
#saveRDS(dmrcoutput, "/home/alreiner/Projects/ffcw/output/DMR/DMRoutput_fdr1.rds")

#fix error where can't connect to experimentHub 
library(ExperimentHub)
setExperimentHubOption(arg="URL", value = 'https://experimenthub.bioconductor.org')
#"PROXY", httr::use_proxy(), 

#main
#results.ranges_09 <- extractRanges(dmrcoutput, genome = "hg19")
#results.ranges_09 


#results.ranges_09_FDR <- extractRanges(dmrcoutput0.05FDR, genome = "hg19")
#2304 DMRs for sex,epi, batch model (get same count when using what=default and what="Beta")
#2362 DMRs for sex,epi,batch,mom race

#DMR results with overlapping genes
#saveRDS(as.data.frame(results.ranges_09), "/home/alreiner/Projects/ffcw/output/DMR/DMR_results_247_SEBP_upCTP.rds")
#saveRDS(as.data.frame(results.ranges_09_FDR), "/home/alreiner/Projects/ffcw/output/DMR/DMR_results_247_SEBP_upCTP_FDR0.05.rds")


# results.ranges_09 <- read_csv("/home/alreiner/Projects/ffcw/output/DMR/DMR_results_247_SEBP_upCTP.csv")
# results.ranges_09 <- results.ranges_09[,-1]
# results.ranges_09 <- makeGRangesFromDataFrame(results.ranges_09)
# 
# 
# #Picking DMR sig threshold
# #View(as.data.frame(results.ranges))
# #View(results.ranges$meandiff)
# #View(as.data.frame(results.ranges_09)%>% filter(no.cpgs > 3,abs(meandiff)>0.09))
# 
# #using sig thresh of > 3 CpG and abs(meandiff)>0.04
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
# write.table(rounded9, file="/home/alreiner/Projects/ffcw/output/DMR/age9_DMRresults_upCTP_11_26.txt",sep=',')
# 
# 
# #Calculate concordance of overlapping genes
# 
# #separate double genes into 2 sep elemts in a list
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
# saveRDS(genessplit_noNA, file=paste0("/home/alreiner/Projects/ffcw/output/DMR/DMR_overlappgenelist_SEBR_upCTP9_",format(Sys.time(), "%b%d"),".rds"))
# 
# 

############################
#Enrichment
#############################

# results.ranges_09 <- readRDS("/home/alreiner/Projects/ffcw/output/DMR/DMR_results_247_SEBP_upCTP.rds")
# results.ranges_09_fdr <-readRDS("/home/alreiner/Projects/ffcw/output/DMR/DMR_results_247_SEBP_upCTP_FDR0.05.rds")
# 
# # #Print out final DMR list for results
# roundedo <- results.ranges_09 %>% 
# mutate_at(vars(maxdiff,meandiff), funs(round(., 3)))
# # #save txt file of results to input into word
# # write.table(rounded9, file="/home/alreiner/Projects/ffcw/output/DMR/age9_DMRresults_upCTP_11_26.txt",sep=',')
# write.table(roundedo[1:15,-c(5,7,8,9)],"/home/alreiner/Projects/ffcw/output/DMR/top15DMRs.txt")
# 
# 
# # #GO term enrichment
# #get GO ontology terms associated with the DMRs
# enrichment_GO_fdr1 <- goregion(makeGRangesFromDataFrame(results.ranges_09), all.cpg = rownames(beta_nogap_noXYmani),
#                               collection = "GO", array.type = "450K",sig.genes = T)
# enrichment_GO_fdr1<- enrichment_GO_fdr1[order(enrichment_GO_fdr1$P.DE),]
# saveRDS(enrichment_GO_fdr1, "/home/alreiner/Projects/ffcw/output/DMR/age9GO_results_fdr1_input.rds")
# 
# #Filter GO terms to biological processes only
# enrichment_GO_fdr1<- enrichment_GO_fdr1[order(enrichment_GO_fdr1$P.DE),]%>%
#   filter(ONTOLOGY=="BP")
# #Further filter to FDR < 0.05
# enrichment_GO_sig_9_filtFDR <- enrichment_GO_fdr1 %>%
#   filter(FDR<0.05) %>%
#   arrange(FDR)%>%
#   mutate_at(vars(P.DE, FDR), funs(round(., 3)))
# #write.table(enrichment_GO_1717sig_15_filt, file="/home/alreiner/Projects/ffcw/analyses/15yo/DMR_1717_GoTerms_BP_11_24.txt",sep=',')
# write.csv(enrichment_GO_2024sig_9_filtFDR, "/home/alreiner/Projects/ffcw/output/DMR/age9GO_results_filtFDR.csv")
