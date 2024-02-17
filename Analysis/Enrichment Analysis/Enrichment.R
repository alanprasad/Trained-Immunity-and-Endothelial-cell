library(ReactomePA)
library(reactome.db)
library(org.Hs.eg.db)

#read the results table from deseq which contains the details of the genes
RNA_res<-read_tsv("results_table")

columns(org.Hs.eg.db)
GENE_ANNO <- ensembldb::select(edb, keys = RNA_res$GENENAME,
                               keytype = "GENENAME",
                               columns = c("GENEID", "GENENAME","SYMBOL",
                                           "SEQNAME","GENESEQSTART",
                                           "GENESEQEND", "GENEBIOTYPE", "ENTREZID"),
                               multiVals="first")

#enrichment analysis of GO, Reactome and KEGG
pathway.reac<- enrichPathway(as.data.frame(RNA_res)$GENEID,universe = RNA_res$GENEID )
head(pathway.reac)

pathway.GO <- enrichGO(as.data.frame(RNA_res)$GENEID,org.Hs.eg.db,ont = "BP",universe = RNA_res$GENEID)
head(pathway.GO)

pathway.KEGG <- enrichKEGG(gene = RNA_res$GENEID,
                           organism = "hsa",
                           pvalueCutoff = 0.05,
                           pAdjustMethod="BH",
                           universe = RNA_res$GENEID)

#Using dotplot to visualize the enrichment result
dotplot(pathway.GO,showCategory=30)
dotplot(pathway.reac,showCategory=25)
dotplot(pathway.KEGG,showCategory=25)

#alternative methods
#enrichment plot
columns(org.Hs.eg.db)
RNA_res$symbol = mapIds(org.Hs.eg.db,
                             keys = row.names(RNA_res),
                             column = "SYMBOL",
                             keytype = "ENSEMBL",
                             multiVals = "first")
RNA_res$entrez =  mapIds(org.Hs.eg.db,
                                    keys = row.names(RNA_res),
                                    column = "ENTREZID",
                                    keytype = "ENSEMBL",
                                    multiVals = "first")

RNA_res$name = mapIds(org.Hs.eg.db,
                           keys = row.names(RNA_res),
                           column = "GENENAME",
                           keytype = "ENSEMBL",
                           multiVals = "first")
all_genes <- as.character(rownames(RNA_res))
signif_res<- RNA_res[RNA_res$padj<0.05 & !is.na(RNA_res$padj), ]
signif_gene <-as.character(rownames(signif_res))

ego<- enrichGO(gene = signif_gene,
               universe = all_genes,
               keyType = "ENSEMBL",
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               pAdjustMethod = "BH",
               qvalueCutoff = 0.5,
               readable = TRUE)

dotplot(ego,showCategory = 12)
enrich_go<- enrichplot::pairwise_termsim(ego)
emapplot(enrich_go)


#creating barplot from enrichment result table using ggplot
library(ggplot2)

#reading data
data<- read.table("Enrichment.txt",header = TRUE,sep = "\t")
#filtering significant terms
filtered_data<-data[data$p.value<=0.05,]

#top enriched terms
top_20_data <- head(filtered_data[order(filtered_data$p.value), ], 20)
p<-ggplot(top_20_data, aes(x = -log10(p.value), y = reorder(gene_set_name, -p.value))) +
  geom_bar(stat = "identity", fill = "darkorange") +
  #geom_text(aes(label = round(p.value, 3), hjust = -0.2), size = 3, color = "black") +
  labs(x = "-log10(p-value)", y = "Enrichment") +
  ggtitle("Enrichment results ") +
  theme_bw() +
  theme(axis.text.y = element_text(hjust = 1,face = "bold",size = 14))
ggsave("Enrichment.png",plot = p,dpi = 300,height = 8,width = 16)
