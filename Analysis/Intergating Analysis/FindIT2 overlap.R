library(tidyverse)
library(readr)
library(plyranges)
library(AnnotationHub)
library(ensembldb)
library(FindIT2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(dplyr)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
Txdb<-EnsDb.Hsapiens.v86


#Load annotationhub
hub = AnnotationHub(localHub=TRUE)
qr <- query(hub, c("EnsDb", "sapiens", "108"))
edb <- qr[[1]]

#Reading ATAC file
anno_atac <- read_tsv(file = "AnnotationDB/ATACseq/HUVEC_TNF_consensus_peaks.mLb.clN.annotatePeaks.txt")

#Create Granges object for ATACseq
ATAC_peak_GR <- anno_atac %>%
  dplyr::rename("start" = Start,
                "end" = End,
                "seqnames" = Chr,
                "feature_id" = PeakID)%>%
  mutate(seqnames = str_replace(seqnames, "chr", "")) %>%
  
  mutate_at(c('seqnames'), as.numeric) %>%
  dplyr::filter(!nchar(seqnames) > 6) %>%
  dplyr::select(feature_id, seqnames, start, end) %>%
  as_granges()

ATAC_peak_GR

seqlevelsStyle(Txdb)<-"UCSC"
seqlevelsStyle(ATAC_peak_GR)<-"UCSC"

#Annotate using the nearest mode in FindIT2
Anno_nearestgene<- mm_nearestGene(peak_GR = ATAC_peak_GR,
                                  Txdb = Txdb)
Anno_nearestgene

#plotting the summary of the analysis
plot_annoDistance(mmAnno = Anno_nearestgene)
getAssocPairNumber(mmAnno = Anno_nearestgene,output_type = "feature_id")
plot_peakGeneAlias_summary(Anno_nearestgene)
plot_peakGeneAlias_summary(Anno_nearestgene,output_type = "feature_id")

# load the RNAseq data and filter based on p-value and log2FC
DF_RNA <- read_tsv(file = "RNAseq/RNAseq_H1H2_TNF_res.tsv")

#Create annotation DF for rnaseq data
annotationDF <- ensembldb::select(edb, keys = DF_RNA$GENEID,
                                  keytype = "GENEID",
                                  columns = c("GENEID", "GENENAME","SYMBOL",
                                              "SEQNAME","GENESEQSTART",
                                              "GENESEQEND", "GENEBIOTYPE", "ENTREZID","EXONID"),
                                  multiVals="first")

#Merge the annotation
DF_RNA_anno <- DF_RNA %>% left_join(annotationDF, by= c("GENEID"="GENEID"))

#Create Granges object for RNAseq
DF_RNA_gr <- DF_RNA_anno %>%
  dplyr::rename("start" = GENESEQSTART,
                "end" = GENESEQEND,
                "seqnames" = SEQNAME,
                "gene" = GENEID)%>%
  mutate(seqnames = str_replace(seqnames, "chr", "")) %>%
  mutate(start-50000,
         end+50000) %>%
  mutate_at(c('seqnames'), as.numeric) %>%
  dplyr::filter(!nchar(seqnames) > 6) %>%
  dplyr::select(gene, seqnames, start, end) %>%
  as_granges()

overlap_peaks<-Anno_nearestgene %>%
  dplyr::filter(gene_id %in% DF_RNA_gr$gene)
overlap_peaks
overlap_peaks_df<-as.data.frame(overlap_peaks)

#writing the results
write.table(overlap_peaks_df,"Overlaps.tsv",sep = "\t",row.names = FALSE,col.names = TRUE)


#Integrating the results by correlating the RNA and ATAC expression values
#extract the normalize DESeq2 result object with count and metadata
gene_exp<-counts(dds,normalized=TRUE)
peak_exp<-counts(dds1,normalized=TRUE)

#Take only the overlapped peaks and gene
gene_exp<-gene_exp[rownames(gene_exp) %in% overlap_peaks$gene_id,]
peak_exp<-peak_exp[rownames(peak_exp) %in% overlap_peaks$feature_id,]

#make sure samples of RNA and ATAC are in same order and follows same names
#changing the samples names of RNA to match with ATAC
new_colNames<-c("CTRL1_EGM_REP5" ,"CTRL1_EGM_REP2" ,"CTRL1_EGM_REP4" ,"CTRL1_EGM_REP1" ,"CTRL1_EGM_REP3",
                "HIT1_IFN_REP4" , "HIT2_IFN_REP2" , "CTRL2_IFN_REP1", "CTRL2_IFN_REP2" ,"HIT2_IFN_REP5",
                "HIT1_IFN_REP3" , "CTRL2_IFN_REP3" ,"HIT1_IFN_REP2" , "HIT2_IFN_REP1" , "HIT2_IFN_REP3",
                "CTRL2_IFN_REP5", "HIT1_IFN_REP1" , "HIT2_IFN_REP4" , "CTRL2_IFN_REP4", "HIT1_IFN_REP5")
colnames(gene_exp)<- new_colNames

#running correlation analysis
Overlap_cor<-peakGeneCor(mmAnno = overlap_peaks,geneScoreMt = gene_exp,peakScoreMt = peak_exp,
                         parallel = FALSE,verbose = TRUE)
Overlap_cor_df<-as.data.frame(Overlap_cor)

#Adding the correlation degree score to the data frame
Overlap_cor_df<- Overlap_cor_df %>%
  mutate(Correlation_degree = case_when(
    cor > 0.1 & cor < 0.39 ~ "Weak+",
    cor > 0.39 & cor < 0.69 ~ "Moderate+",
    cor > 0.69 & cor < 0.9 ~ "Strong+",
    cor > 0.9 ~ "perfect+",
    cor < (-0.1) & cor > (-0.39) ~ "Weak-",
    cor < (-0.39) & cor > (-0.69) ~ "Moderate-",
    cor < (-0.69) & cor > (-0.9) ~ "Strong-",
    cor < (-0.9) ~ "perfect-",
    TRUE ~ "No correlation"
  ))

#Saving the file
write.table(Overlap_cor_df,"overlap_cor_table.tsv",sep = "\t",col.names = TRUE,row.names = FALSE) 

#creating correlation plot for genes and peaks
#extracting the samples names to label names
cols <- dput(colnames(gene_exp))
newc <- sub("_.*", "", cols)

#creating sublist of highly correlation score
subset(Overlap_cor,cor > 0.69) %>%
    getAssocPairNumber()

Plot_corr <- plot_peakGeneCor(mmAnnoCor =  subset(Overlap_cor,cor > 0.69),
                       select_gene = "GENENAME")+
  geom_point(aes(color=c(newc,newc)))+
  geom_smooth(method = "lm",formula = y~x)+
  labs(color="Sample")

Plot_corr


#plotting the count for correlation degree
plot<-ggplot(Overlap_cor_df,aes(x=Correlation_degree,fill=Correlation_degree))+
  geom_bar()+
  labs(x="Correlation degree",
       y="Number of Peaks")+
  geom_text(stat = "count",aes(label=..count..),vjust=-0.5)+
  scale_x_discrete(labels=c("<-0.39 & >-0.69",">0.39 & <0.69",">-0.1 & <0.1",">0.9","<-0.69 & >-0.9",">0.69 & <0.9","<-0.1 & >-0.39",">0.1 & <0.39"))+
  theme_minimal()
plot
ggsave("Correlation_degree_barplot.png",plot = plot,dpi = 600,width = 10,height = 8)



