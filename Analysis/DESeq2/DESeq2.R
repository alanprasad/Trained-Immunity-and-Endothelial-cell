library(DESeq2)

#loading data
vein_data <- readRDS("RNAseq/HUVEV_IFN.rds")
head(vein_data)

#reading annotation file
AnnotationDB <- read.delim("AnnotationDB/annotationDF2.tsv",header = TRUE)
head(AnnotationDB)

#row name in colData matches to column names in count data
all(colnames(vein_data)%in%rownames(colData(vein_data)))

#are they in the same order?
all(colnames(vein_data)==rownames(colData(vein_data)))


#creating data frame for count
count_data<- assay(vein_data)
Metadata <- data.frame(colData(vein_data))

#Isolating different gene types
filtered_gene<- AnnotationDB[AnnotationDB$GENEBIOTYPE == "lncRNA", ]
miRNA<- AnnotationDB[AnnotationDB$GENEBIOTYPE == "miRNA", ]
protein_coding<- AnnotationDB[AnnotationDB$GENEBIOTYPE == "protein_coding", ]

#Filtering the long non-coding gene
count_TNF<-count_data[rownames(count_data)%in%filtered_gene$GENEID,]
#optional steps for other gene type
count_miRNA<-count_data[rownames(count_data)%in%miRNA$GENEID,]
count_protein<-count_data[rownames(count_data)%in%protein_coding$GENEID,]

#Adding the missing fields in metadata
Metadata<- Metadata %>%
  mutate(REP = case_when(
    donor=="SD125"~ "REP1",
    donor=="SD126"~ "REP2",
    donor=="SD263"~ "REP3",
    donor=="SD098"~ "REP4",
    donor=="SD450"~ "REP5",
  ))


#creating DEseq2 object
dds<- DESeqDataSetFromMatrix(countData = count_TNF,
                             colData = Metadata,
                             design = ~ donor + condition)
dds

#perform pre-filtering of the data
dds<- dds[rowSums(counts(dds)>=10)>=3,]

#setting the levels of design to compare the level with rest of the levels
dds$condition<- relevel(dds$condition,ref = "2hit")

#Differential expression analysis
dds<- DESeq(dds)

#Extracting the results on desired pvalue threshold
res<- results(dds,alpha = 0.05)
resultsNames(dds)

#specified the coefficient or contrast
#Repaeting the above step with selected condition camparison extracting results
res05_2hvs1h<- results(dds,contrast = c("condition","2hit","1hit"),alpha = 0.05)
res05_2hvs1h_df<- subset(res05_2hvs1h,padj<0.05)
res05_2hvsctrl2<- results(dds,contrast = c("condition","2hit","CTRL2"),alpha = 0.05)
res05_2hvsctrl2_df<- subset(res05_2hvsctrl2,padj<0.05)
res05_1hvsctrl1<- results(dds,contrast = c("condition","1hit","CTRL1"),alpha = 0.05)
res05_1hvsctrl1_df<-subset(res05_1hvsctrl1,padj<0.05)

#changing the datatype of the results and adding the geneid to rownames
res05_2hvsctrl2_df<-as.data.frame(res05_2hvsctrl2_df)
res05_2hvsctrl2_df <- rownames_to_column(res05_2hvsctrl2_df, var = "GENEID")
res05_1hvsctrl1_df<-as.data.frame(res05_1hvsctrl1_df)
res05_1hvsctrl1_df <- rownames_to_column(res05_1hvsctrl1_df, var = "GENEID")
res05_2hvs1h_df<-as.data.frame(res05_2hvs1h_df)
res05_2hvs1h_df <- rownames_to_column(res05_2hvs1h_df, var = "GENEID")

#saving the results in tsv format
write.table(res05_2hvsctrl2_df,file = "RNAseq/RNAseq_H2C2_TNF_res.tsv",sep = "\t",row.names=FALSE,col.names = TRUE)
write.table(res05_1hvsctrl1_df,file = "RNAseq/RNAseq_H1C1_TNF_res.tsv",sep = "\t",row.names=FALSE,col.names = TRUE)
write.table(res05_2hvs1h_df,file = "RNAseq/RNAseq_H1H2_TNF_res.tsv",sep = "\t",row.names=FALSE,col.names = TRUE)
