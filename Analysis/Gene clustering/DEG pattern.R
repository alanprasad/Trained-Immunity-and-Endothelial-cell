library(DESeq2)
library(dplyr)
library(readr)

# Load data
vein_data <- readRDS("RNAseq/HUVEV_IFN.rds")
AnnotationDB <- read.delim("AnnotationDB/annotationDF2.tsv", header = TRUE)

# Prepare metadata
Metadata <- colData(vein_data)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = assay(vein_data),
                              colData = Metadata,
                              design = ~ donor + condition)

# Perform pre-filtering
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]

#Reading the result table of DESeq2
RNA<-read_tsv("results_DESeq2.tsv")

#variance stabilizing transformation and the rlog
vsd<- vst(dds,blind = FALSE)
rld<- rlog(dds,blind = FALSE)

#Extracting the count table of significant genes
Count<-assay(rld)
Count<- Count[rownames(ma) %in% RNA$GENEID,]

#creating metadata
des<-as.data.frame(colData(rld))

#Relevel the conditions
des$condition<-factor(des$condition,levels = c("CTRL1","1hit","CTRL2","2hit"))

#running DEG pattern analysis
res_pattern<- degPatterns(ma,des,time = "condition",col="type",reduce = TRUE,minc = 5,scale = TRUE,plot = TRUE,pattern = NULL)

#extrenal plotting of DEG using ggplot
ggplot(res_pattern[["normalized"]],
       aes(condition, value,color=cluster)) +
  #geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.5))+
  theme_classic()+
  # change the method to make it smoother
  geom_smooth(aes(group=cluster))

#Heatmap analysis
#Normalized count table of DESeq object
select <- order(rowMeans(counts(dds,normalized = TRUE)),decreasing = TRUE)

#selecting the needed label columns
df<- as.data.frame(colData(dds)[,c("condition","donor")])

#selecting the significant genes count values
table<-assay(dds)[select,]
table<-table[rownames(table) %in% RNA$GENEID,]

#order the table according to the stimulation condition
table<- table[,c(1,5,12,16,17,3,7,9,14,20,2,6,13,18,19,4,8,10,11,15)]

pheatmap(table,cluster_rows = FALSE,show_rownames = FALSE,
         cluster_cols = FALSE,annotation_col = df,scale = "row")