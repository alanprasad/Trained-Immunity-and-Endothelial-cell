install.packages("UpSetR")
library(UpSetR)

setwd("C:/Users/Alan")

#read the list of significant genes results
hih2t<-readLines("Result1.txt")
hic2t<-readLines("Result2.txt")
hic1t<-readLines("Result3.txt")
hih2i<-readLines("Result4.txt")
hic2i<-readLines("Result5.txt")
hic1i<-readLines("Result6.txt")

#creating Upset plot
upsetlist<-list(TNF_H2H1=hih2t,TNF_H2C2=hic2t,TNF_H1C1=hic1t,IFN_H2H1=hih2i,IFN_H2C2=hic2i,IFN_H1C1=hic1i)
colors<-c("red","blue","orange","violet","green","purple")
upset(fromList(upsetlist),nsets = 6,sets.x.label="Number of DEGs",main.bar.color = "brown",sets.bar.color = colors,
      set_size.show=TRUE,set_size.scale_max = 1000)
