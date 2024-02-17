#Reading the annotation file of ATAC-seq data
Annot_ATAC <- read.delim("AnnotationDB/ATACseq/HUVEC_IFN_consensus_peaks.mLb.clN.annotatePeaks.txt")

#Read the DESeq2 results of ATAC-seq
Result1<-read_tsv("deseq2_reslt.tsv")
Result2<-read_tsv("deseq2_reslt.tsv")
Result3<-read_tsv("deseq2_reslt.tsv")


Anno_h1c1 <-Annot_ATAC[Annot_ATAC$PeakID %in% rownames(Result1),]
Anno_h2c2 <-Annot_ATAC[Annot_ATAC$PeakID %in% rownames(Result2),]
Anno_h2h1 <-Annot_ATAC[Annot_ATAC$PeakID %in% rownames(Result3),]

#creating Grange object of ATAC annotation file
peaks_full<-makeGRangesFromDataFrame(Annot_ATAC,keep.extra.columns = TRUE)

#extracting the nearby gene of all peaks by setting the window size
peakAnno_full<-annotatePeak(peaks_full,tssRegion = c(-5000,5000),TxDb = txdb,annoDb = "org.Hs.eg.db" )
peak.bed_full<-as.data.frame(peakAnno_full)

#creating the Grange object for Deseq2 results
peaks<-makeGRangesFromDataFrame(Anno_h1c1,keep.extra.columns = TRUE)

#extracting the nearby gene of the significant peaks by setting the window size
peakAnno<-annotatePeak(peaks,tssRegion = c(-5000,5000),TxDb = txdb,annoDb = "org.Hs.eg.db" )
peak.bed<-as.data.frame(peakAnno)

#creating peak annotation plots of selected peaks
upsetplot(peakAnno,vennpie=TRUE)
plotDistToTSS(peakAnno)
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
