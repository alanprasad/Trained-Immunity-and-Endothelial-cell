library(readr)

# Define the file names
file_names <- c("results_DESeq2_1.tsv", "results_DESeq2_2.tsv", "results_DESeq2_3.tsv")

# Initialize an empty list to store the row names from each file
result_list <- list()

# Read each file and extract row names
for (file in file_names) {
  # Read the file
  result_df <- read_tsv(file)
  
  # Extract row names and store them in the list
  result_list[[file]] <- rownames(result_df)
}

# Viewing the list of row names
print(result_list)


#Venn diagram
venn.diagram(x=result_list,
             category.names = c("TNF-a 2HIT vs 1HIT","TNF-a 2HIT vs CTRL2","TNF-a 1HIT vs CTRL1"),
             filename = "Venn_HUVEV_TNF.png",
             output = TRUE,
             imagetype = "png",
             height = 680,
             width = 680,
             resolution = 500,
             compression = "lzw",
             lwd=1,
             col=c("#440154ff","#21908dff","#fde725ff"),
             fill=c(alpha("#440154ff",0.3),alpha("#21908dff",0.3),alpha("#fde725ff",0.3)),
             cex=0.5,
             fontfamily="sans",
             cat.cex=0.3,
             cat.default.pos="outer",
             cat.pos=c(-27,27,135),
             cat.dist=c(0.055,0.055,0.085),
             cat.fontfamily="sans",
             cat.col=c("#440154ff","#21908dff","#fde725ff"),
             rotation=1
)
