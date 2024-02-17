if(!require('remotes')) {
  install.packages('remotes')
  library('remotes')
}
BiocManager::install("hfang-bristol/XGR", dependencies=T, force=T)
install.packages("XGR", repos="http://R-Forge.R-project.org")
install.packages("XGR", repos="https://cran.rstudio.com")
install.packages("XGR")


setwd("C:/Users/Alan")

# Load the library
library(biomaRt)
library(XGR)
library(readr)
library(ggplot2)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
snps<-SNPlocs.Hsapiens.dbSNP144.GRCh38

#reading RNA data
RNA <- read_tsv(file = "RNAseq/RNAseq_res.tsv")

#reading backgorund genes
RNA_BG <- read_tsv(file = "RNAseq/RNAseq_res.tsv")

#Rename the rownames as GENEID
RNA <- rownames_to_column(RNA, var = "GENEID")
BG_RNA <- rownames_to_column(RNA_BG, var = "GENEID")

#Create annotation DF for rnaseq data
annotationDF <- ensembldb::select(edb, keys = RNA$GENEID,
                                  keytype = "GENEID",
                                  columns = c("GENEID", "GENENAME","SYMBOL",
                                              "SEQNAME","GENESEQSTART",
                                              "GENESEQEND", "GENEBIOTYPE", "ENTREZID"),
                                  multiVals="first")

#Merge the annotation
TNF_RNA_anno <- RNA %>% left_join(annotationDF, by= c("GENEID"="GENEID"))

#reading atac-seq data
#Load the annotattion for the atacseq
anno_atac <- read_tsv(file = "AnnotationDB/ATACseq/HUVEC_IFN_consensus_peaks.mLb.clN.annotatePeaks.txt")

#Creating RNA-seq object 
TNF_RNA_gr <- TNF_RNA_anno %>%
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

#Create Granges object for ATACseq
TNF_atac_gr <- anno_atac %>%
  dplyr::rename("start" = Start,
                "end" = End,
                "seqnames" = Chr,
                "feature_id" = PeakID)%>%
  mutate(seqnames = str_replace(seqnames, "chr", "")) %>%
  
  mutate_at(c('seqnames'), as.numeric) %>%
  dplyr::filter(!nchar(seqnames) > 6) %>%
  dplyr::select(feature_id, seqnames, start, end) %>%
  as_granges()

# find peaks that overlap with gene postions.
plyranges::find_overlaps_directed(TNF_RNA_gr, TNF_atac_gr )
TNF_overlaps<-plyranges::find_overlaps_directed(TNF_RNA_gr, TNF_atac_gr )
#seqlevelsStyle(overlaps)<-"UCSC"
TNF_overlaps_df<-as.data.frame(TNF_overlaps)

#creating the overlapping peaks Grange
TNF_atac_gr<- TNF_atac_gr %>%
  dplyr::filter(feature_id %in% TNF_overlaps_df$feature_id)


snpcount(snps)
colnames(snps)

#extracting the SNPs of the overlapped peaks
TNF_overlaps_snps<-snpsByOverlaps(snps,TNF_atac_gr)
TNF_overlaps_snps_df<-as.data.frame(TNF_overlaps_snps)
write.table(TNF_overlaps_snps_df,"Overlap_SNPs.tsv",sep = "\t",row.names = FALSE,col.names = TRUE)


#reading interesing SNPs data
df_snps<-read_tsv("Overlap_SNPs.tsv")
df_snps_bg<-read_tsv("Overlap_background_SNPs.tsv")

#list of snps
overlap_snps<-df_snps$RefSNP_id
overlap_snps_bg<-df_snps_bg$RefSNP_id

# SNP-based enrichment analysis using GWAS Catalog traits (mapped to EF)
#perform enrichment analysis
eTerm <- xEnricherSNPs(data=overlap_snps, ontology="EF",background = overlap_snps_bg,
                       path.mode=c("all_paths"),p.adjust.method = "BH",size.range = c(1,2000),min.overlap = 1)

#view enrichment results for the top significant terms
table_eTerm<-xEnrichViewer(eTerm,top_num = 22)

#save enrichment results to the file called 'EF_enrichments.txt'
res <- xEnrichViewer(eTerm, top_num=length(eTerm$adjp), sortBy="adjp",
                     details=TRUE)
output <- data.frame(term=rownames(res), res)
utils::write.table(output, file="EF_enrichments_IFN.tsv", sep="\t",
                   row.names=FALSE)

# barplot of significant enrichment results
bp <- xEnrichBarplot(eTerm, top_num=15, displayBy="adjp")
bp
ggsave(filename = "XGR_ontology_H2H1_RNA_IFN.png",plot = bp,height = 6,width = 8,dpi = 600)

# visualise the top significant terms in the ontology hierarchy
# color-code terms according to the adjust p-values (taking the form of 10-based negative logarithm)
oh<-xEnrichDAGplot(eTerm, top_num=15, displayBy="adjp",
                   node.info=c("full_term_name"))
# color-code terms according to the z-scores
xEnrichDAGplot(eTerm, top_num=15, displayBy="zscore",
               node.info=c("full_term_name"))
xEnrichGGraph(eTerm)