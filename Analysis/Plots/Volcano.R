library(ggplot2)
library(readr)
library(org.Hs.eg.db)

#Reading the result table of DESeq2
RNA<-read_tsv("results_DESeq2.tsv")
View(RNA)

#Adding genename to the result table
RNA$GENENAME =  mapIds(org.Hs.eg.db,
                                    keys = row.names(RNA),
                                    column = "GENENAME",
                                    keytype = "GENEID",
                                    multiVals = "first")

#Adding the remaing column to seperate the significant genes
RNA$diffexpressed<-"NOT"
RNA$diffexpressed[RNA$log2FoldChange>0&RNA$padj<0.05]<-"UP"
RNA$diffexpressed[RNA$log2FoldChange<0&RNA$padj<0.05]<-"DOWN"
top15gene<-head(RNA[order(RNA$padj),"GENENAME"],20)
RNA$delabel<-ifelse(RNA$GENENAME%in%top15gene,RNA$GENENAME,NA)

#Creating volcano plot
#setting up the theme
theme_set(theme_bw(base_size = 20)+
            theme(
              axis.title.y = element_text(face = "bold",margin = margin(0,20,0,0),size = rel(1,1),colour = "black"),
              axis.title.x = element_text(hjust = 0.5,face = "bold",margin = margin(20,0,0,0),size = rel(1,1),colour = "black"),
              plot.title = element_text(hjust = 0.5)
            ))

#volcano plot
plot_volcano<-ggplot(data=RNA,aes(x=log2FoldChange,y= -log10(padj),col=diffexpressed,label=delabel))+
  geom_vline(xintercept = c(-0.5,0.5),col="black",linetype="dashed")+
  geom_hline(yintercept = c(0.5),col="black",linetype="dashed")+
  geom_point(size=3)+
  scale_color_manual(values = c("blueviolet","grey","darkorange"),
                     labels=c("DownRegulated","Not Significant","UpRegulated"))+
  coord_cartesian(xlim=c(-8,10))+
  scale_x_continuous(breaks = seq(-8,10,4))+
  labs(color="Conditions")+
  ggtitle("Volcano plot") +
  geom_text_repel(max.overlaps = Inf)

ggsave("volcano.png",plot = plot_s,dpi = 300,width = 12,height = 8)


#Enhanced volcano plot
EnhancedVolcano(RNA,
                lab = rownames(RNA),
                x="log2FoldChange",y ="pvalue",
                title = 'TNF-a 2HIT vs 1HIT',
                xlim = c(-4,6),
                ylim = c(0,20),
                labSize = 6.0,
                colAlpha = 1)

p2<-EnhancedVolcano(RNA,
                    lab = rownames(RNA),
                    x="log2FoldChange",y ="pvalue",
                    title = 'TNF-a 1HIT vs CTRL2',
                    xlim = c(-4,6),
                    ylim = c(0,20),
                    labSize = 6.0,
                    colAlpha = 1)
grid.arrange(p1,p2,
             ncol=2,
             top = textGrob('EnhancedVolcano',
                            just = c('center'),
                            gp = gpar(fontsize = 32)))


