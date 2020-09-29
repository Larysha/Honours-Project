#### Data obtained from the ChIP-Atlas Peak Calling tool. Alignments with a significance threshold of 
#### 50 were set so that peaks with Q value < 1E-05 are shown on genome browser IGV. The BED file were 
#### extracted from the IGV genome browser and cleaned with the help of Shinya Oki from Kyoto University, Japan

##### Getting and Reading the Data

###### load the dplyr package 
###### set working directory to the TF folder that contains the BED files
library(dplyr)
wd <- setwd("~/R/Honours Project/Honours-Project/data BED/TF")

###### read the data into an R object
breast.bed <- read.table("breastTF.bed", sep = "\t", header = FALSE)
cardio.bed <- read.table("cardioTF.bed", sep = "\t", header = FALSE)
kidney.bed <- read.table("kidneyTF.bed", sep = "\t", header = FALSE)
prostate.bed <- read.table("prostateTF.bed", sep = "\t", header = FALSE)

###### quick look at the data in breast (repeat for cardio, kidney and prostate)
str(breast.bed)

###### assign variable names to the data objects
colnames(breast.bed) <- c("Chromosome", "Start", "End", "MACS2 Value", "ID", "TF", "Cell Type")
colnames(cardio.bed) <- c("Chromosome", "Start", "End", "MACS2 Value", "ID", "TF", "Cell Type")
colnames(kidney.bed) <- c("Chromosome", "Start", "End", "MACS2 Value", "ID", "TF", "Cell Type")
colnames(prostate.bed) <- c("Chromosome", "Start", "End", "MACS2 Value", "ID", "TF", "Cell Type")

###### Find out all the unique cell lines for the ChIP-Seq results in breast
unique(breast.bed$`Cell Type`)
###### contains 95 unique cell types, all cancerous to varying degrees 
unique(cardio.bed$`Cell Type`)
###### contains 16 unique cell types for cardio.bed
###### contains 30 unique cell types for kidney.bed
###### contains 22 unique cell types for prostate.bed

##### Peak Annotation and Visualisation 

##### loading ChIPseeker and associated packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
###### This exposes an annotation databases generated from UCSC by exposing these as TxDb objects
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
###### Limit this object to chromosome 2
seqlevels(txdb) <- "chr2"
library(clusterProfiler)

###### create a files object containing all 4 bed files/ tissue types
files = c('breastTF.bed', 'cardioTF.bed', 'kidneyTF.bed', 'prostateTF.bed')
###### make a list of these files 
samplefiles <- as.list(files)
###### name the values of the list
names(samplefiles) <- c("breast", "cardio", "kidney", "prostate")

###### annotate the peak data 1000 base pairs upstream and downstream of the TSS for each file in samplefiles
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb,tssRegion=c(-1000, 1000), verbose=FALSE)

###### plot feature distribution of peaks according to tissue type
plotAnnoPie(peakAnnoList[["breast"]], col = blues9)
plotAnnoPie(peakAnnoList[["cardio"]], col = blues9)
plotAnnoPie(peakAnnoList[["kidney"]], col = blues9)
plotAnnoPie(peakAnnoList[["prostate"]], col = blues9)

###### plot the location of each peak relative to the TSS
###### the x-axis gives the percentage of sites, while the color represents the distance from the TSS
plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci \n relative to TSS")

###### store the annotations for each tissue into a data frame
###### retrieve the information for the peakAnnoList object
breast_annot <- as.data.frame(peakAnnoList[["breast"]]@anno)
cardio_annot <- as.data.frame(peakAnnoList[["cardio"]]@anno)
kidney_annot <- as.data.frame(peakAnnoList[["kidney"]]@anno)
prostate_annot <- as.data.frame(peakAnnoList[["prostate"]]@anno)

###### write these data frames into files 
write.table(breast_annot, file="breast_annotation.txt", sep = "\t", row.names = FALSE)
write.table(cardio_annot, file="cardio_annotation.txt", sep = "\t", row.names = FALSE)
write.table(kidney_annot, file="kidney_annotation.txt", sep = "\t", row.names = FALSE)
write.table(prostate_annot, file="prostate_annotation.txt", sep = "\t", row.names = FALSE)

##### Functional Enrichment Analysis










---------------------------------------------------------------------------------
### load the peak data (from column 4) and store it into a GRanges object
peak <- readPeakFile(files[[4]])

### use the MACS score to find the peak locations in PXDN and visualize this
covplot(peak, weightCol="V4",chrs="chr2")



### prepare the TSS regions defined as the regions flanking the TSS
### alignment of the peaks that map to these regions used to create the tagMatrix
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)

### create a heatmap of the ChIP-Seq data binding to the TSS regions
peakHeatmap(files[[4]], TxDb=txdb, upstream=3000, downstream=3000, color="red")


