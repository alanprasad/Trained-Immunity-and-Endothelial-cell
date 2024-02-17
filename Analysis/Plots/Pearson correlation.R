library(readr)

# Define the file names
file_names <- c("results_DESeq2_1.tsv", "results_DESeq2_2.tsv", "results_DESeq2_3.tsv")

# Initialize an empty list to store the data frames
result_list <- list()

# Read each file into the list
for (file in file_names) {
  result_list[[file]] <- read_tsv(file)
}

#defining conditions
conditions_TNF<-c("2hitvs1hit","2hitvsctrl2","1hitvsctrl1")

#creating empty matrix
corr_matrix<- matrix(0,nrow = length(conditions_TNF),ncol = length(conditions_TNF),dimnames = list(conditions_TNF,conditions_TNF))

for (i in 1:(length(conditions_TNF)-1)) {
  for (j in (i+1):length(conditions_TNF)){
    cat(conditions_TNF[i],"vs",conditions_TNF[j]," ")
    
    shared_genes <- intersect(rownames(result_list[[i]]),rownames(result_list[[j]]))
    
    df1_shared<- result_list[[i]][shared_genes, ]
    df2_shared<- result_list[[j]][shared_genes, ]
    
    #calculating correlation coefficients
    pearson_corr <-cor(df1_shared$log2FoldChange,df2_shared$log2FoldChange,method = "pearson")
    #storing correlation matrix
    corr_matrix[i,j]<- pearson_corr
    corr_matrix[j,i]<- pearson_corr
    
    df<- data.frame(EnsemblID = shared_genes,
                    method1_logfc = df1_shared$log2FoldChange,
                    method2_logfc = df2_shared$log2FoldChange)
    
    #plotting
    plot<-ggplot(data = df,aes(x=method1_logfc,y = method2_logfc))+
      geom_point(size=2)+
      labs(title=paste("comparison of log2FoldChange",conditions_TNF[i],"vs",conditions_TNF[j]),
           x = paste0(conditions_TNF[i],"log2FoldChange"), y = paste0(conditions_TNF[j],"log2FoldChange"))+
      theme(text = element_text(size = 14))+
      coord_cartesian(ylim = c(-5,10),xlim=c(-5,10))+
      geom_text(x = max(df$method2_logfc),y= min(df$method1_logfc),
                label=paste0("shared genes: ",length(shared_genes),"\nPearson's correlation index: ", round(pearson_corr,2)),hjust=1,vjust=0,size=5)
    ggsave(paste0("plot/",conditions_TNF[i],"vs",conditions_TNF[j],"_scatterplot.png"),plot=plot,dpi=300,width = 8,height = 8)                      
  }
  
}

#pearson heatmap
corr_matrix<-corr_matrix[order(rownames(corr_matrix)),order(rownames(corr_matrix))]
corr_df<-as.data.frame(corr_matrix)
corr_df<-replace(corr_df,corr_df==0,1)
corr_df$condition<-row.names(corr_df)

#melt dataframe into long format
melted_matrix<-melt(corr_df,id.vars = "condition",variable.name = "comparison",value.name = "correlation")

#plotting heat map
heatmap_plot<- ggplot(melted_matrix,aes(x=condition,y=comparison,fill=correlation))+
  geom_tile(color="black",linewidth=0.5)+
  geom_text(aes(label=sprintf("%.2f",correlation)),color="black",size=4)+
  scale_fill_gradient(low = "beige",high="darkorange")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size = 12),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(size = 16,face = "bold"))+
  ggtitle("Pearson's correlation")
ggsave("pearson heat map.png",plot = heatmap_plot,dpi = 300,height = 6,width = 6)
