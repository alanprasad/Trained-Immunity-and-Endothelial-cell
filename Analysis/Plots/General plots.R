#calculating mean for each sample
sample_mean<- rowMeans(count_data)
plot_data<- data.frame(mean = sample_mean)

#plotting density plot for data distribution
plot <- ggplot(plot_data,aes(x=mean))+
  geom_density(fill="red",alpha=0.5)+
  labs(x="Mean value of genes across all samples",y="Density")+
  scale_x_log10(labels=scales::comma_format())+
  ggtitle("Data distribution density plot")
ggsave("Density plot before norm.png",plot = plot,dpi = 300,width = 12)

#calculating the mean for each sample
gene_mean <- colMeans(count_data)
plot_scatter<- data.frame(sample= colnames(count_data),mean = gene_mean)

#combining pheno and count data
plot_scatter<- merge(x= plot_scatter,y=Metadata,by.x="sample",by.y = "sampleID")
#y_limits<- c(0,2000)
#plotting scatterplot
plot <-ggplot(plot_scatter,aes(x=sample,y=mean,fill=condition))+
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  labs(x="sample",y="Mean Value")+
  #scale_y_continuous(limits = y_limits)
  ggtitle('Scatterplot of mean sample')
ggsave("scatter plot.png",plot = plot,dpi = 300)

#performing PCA on filtered count data
pca_count <- t(count_data)
pca_count<- pca_count[,apply(pca_count,2,function(x)!all(x==0))]

#PCA
pcaResult<- prcomp(pca_count,scale. = TRUE)
pc1<- pcaResult$x[, 1]
pc2<- pcaResult$x[, 2]
pcaData <- data.frame(PC1 = pc1,PC2=pc2,sampleID=rownames(pca_count))
pcaData<- merge(pcaData,coldata, by = "sampleID")

#ploting PCA
plot1 <- ggplot(pcaData,aes(x= PC1,y=PC2,color=condition))+
  geom_point()+
  theme_bw()+
  labs(x="PC1",y="PC2",title = "PCA Plot")
ggsave("PCA.png",plot = plot1,dpi = 300)

