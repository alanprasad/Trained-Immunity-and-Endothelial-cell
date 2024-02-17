library(tidyverse)
library(readr)
library(plyranges)
library(AnnotationHub)
library(ensembldb)
library(FindIT2)

#Load annotationhub
hub = AnnotationHub(localHub=TRUE)
qr <- query(hub, c("EnsDb", "sapiens", "108"))
edb <- qr[[1]]

#devtools::install_version("dbplyr", version = "2.3.4")

# load the RNAseq data and filter based on p-value and log2FC
DF_RNA <- read_tsv(file = "RNAseq/RNAseq_H1H2_TNF_res.tsv")
#DF_RNA <- DF_RNA %>%
#dplyr::filter(padj < 0.05 & abs(log2FoldChange > 0.5))

#Create annotation DF for rnaseq data
annotationDF <- ensembldb::select(edb, keys = DF_RNA$GENEID,
                                  keytype = "GENEID",
                                  columns = c("GENEID", "GENENAME","SYMBOL",
                                              "SEQNAME","GENESEQSTART",
                                              "GENESEQEND", "GENEBIOTYPE", "ENTREZID"),
                                  multiVals="first")

#Merge the annotation
DF_RNA_anno <- DF_RNA %>% left_join(annotationDF, by= c("GENEID"="GENEID"))


# Read the DE results of atacseq
DF_atac <- read.table("ATACseq/ATAC_H1C1_IFN_Res.tsv")
new_atac<- DF_atac %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange > 1)|abs(log2FoldChange<(-1)))
new_atac <- rownames_to_column(new_atac, var = "PeakID")
#Load the annotattion for the atacseq
anno_atac <- read_tsv(file = "AnnotationDB/ATACseq/HUVEC_TNF_consensus_peaks.mLb.clN.annotatePeaks.txt")

#Merge the atacseq with annotation 
DF_atac_anno <- new_atac %>% 
  left_join(anno_atac, by = "PeakID") 


#Create Granges object for ATACseq
DF_atac_gr <- anno_atac %>%
  dplyr::rename("start" = Start,
                "end" = End,
                "seqnames" = Chr,
                "feature_id" = PeakID)%>%
  mutate(seqnames = str_replace(seqnames, "chr", "")) %>%
  
  mutate_at(c('seqnames'), as.numeric) %>%
  dplyr::filter(!nchar(seqnames) > 6) %>%
  dplyr::select(feature_id, seqnames, start, end) %>%
  as_granges()

DF_atac_gr

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

DF_RNA_gr

# find peaks that overlap with gene postions.
plyranges::find_overlaps_directed(DF_RNA_gr, DF_atac_gr )
overlaps<-plyranges::find_overlaps_directed(DF_RNA_gr, DF_atac_gr )
seqlevelsStyle(overlaps)<-"UCSC"
overlaps_df<-as.data.frame(overlaps)