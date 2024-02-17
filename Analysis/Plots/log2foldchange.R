#Log2 foldchange comparison of results of DEGs
#read results 
Result1<-read_tsv("result1.tsv")
Result2<-read_tsv("result2.tsv")

#adding gene name to the results
Result1$GENENAME = mapIds(org.Hs.eg.db,
                        keys = row.names(Result1),
                        column = "GENENAME",
                        keytype = "GENEID",
                        multiVals = "first")
Result2$GENENAME = mapIds(org.Hs.eg.db,
                          keys = row.names(Result2),
                          column = "GENENAME",
                          keytype = "GENEID",
                          multiVals = "first")
#correlation scatter
#compare_table should contain two results table
compare_table <- left_join(results1,results2,by ="Row.names") %>%
  filter(!is.na(padj.x)) %>%
  mutate(significant=case_when(
    padj.x<0.05 & padj.y<0.05 ~ "Both",
    padj.x<0.05 & padj.y>=0.05 ~ "HIT2vsHIT1",
    padj.x>=0.05 & padj.y<0.05 ~ "HIT2vsCTRL2",
    padj.x>=0.05 & padj.y>=0.05 ~ "None"
  ))
pearson_corr <-cor(compare_table$log2FoldChange.x,compare_table$log2FoldChange.y,method = "pearson")

plot<-ggplot(compare_table,aes(x=log2FoldChange.x,y=log2FoldChange.y))+
  geom_point(aes(col=significant))+
  scale_color_manual(values=c(Both="green",None="grey",HIT2vsHIT1="blue",HIT2vsCTRL2="red"))+
  ggrepel::geom_text_repel(data = filter(compare_table,significant=="Both"),aes(label=GENENAME.x))+
  geom_text(x = max(compare_table$log2FoldChange.x),y= min(compare_table$log2FoldChange.y),label=paste0("Pearson's correlation index: ", round(pearson_corr,2)),hjust=1,vjust=0,size=5)+
  ylab("LogFoldChange")+
  xlab("LogFoldChange")+
  ggtitle("Logfoldchange correlation Between conditions")
ggsave("Slog2foldchange_scatterplot.png",plot=plot,dpi = 300,height = 10,width = 14)
