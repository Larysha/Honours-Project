#### Data obtained from the ChIP-Atlas Peak Calling tool. ChIP-seq read alignments with a significance threshold of 
#### 50 were set so that peaks with Q value < 1E-05 are shown on genome browser IGV. The MACS2 peak calling software 
#### generated the BED files, which were obtained with the help of Shinya Oki from Kyoto University, Japan.
#### These peaks can also be interactively viewed on the IGV genome browser using the Peak Calling ChIP-Atlas Tool.
-----------------------------------------------------------------------------------------------------------------

  ##### Getting and Reading the Data
---------------------------------------
  
###### load the dplyr package 
###### set working directory to the TF folder that contains the BED files
library(dplyr)
wd <- setwd("~/R/Honours Project/Data and Scripts/data BED/TF")

###### read the data into an R object
breast.bed <- read.table("breastTF.bed", sep = "\t", header = FALSE)
cardio.bed <- read.table("cardioTF.bed", sep = "\t", header = FALSE)
kidney.bed <- read.table("kidneyTF.bed", sep = "\t", header = FALSE)
prostate.bed <- read.table("prostateTF.bed", sep = "\t", header = FALSE)

###### quick look at the data in breast (repeat for cardio, kidney and prostate)
str(breast.bed)

###### assign variable names to the data objects
colnames(breast.bed) <- c("Chromosome", "Start", "End", "-log10(MACS2 Q Value)", "ID", "TF", "Cell Type")
colnames(cardio.bed) <- c("Chromosome", "Start", "End", "-log10(MACS2 Q Value)", "ID", "TF", "Cell Type")
colnames(kidney.bed) <- c("Chromosome", "Start", "End", "-log10(MACS2 Q Value)", "ID", "TF", "Cell Type")
colnames(prostate.bed) <- c("Chromosome", "Start", "End", "-log10(MACS2 Q Value)", "ID", "TF", "Cell Type")

###### Find out all the unique cell lines for the ChIP-Seq results in breast
unique(breast.bed$`Cell Type`)
###### contains 95 unique cell types, all cancerous to varying degrees 
unique(cardio.bed$`Cell Type`)
###### contains 16 unique cell types for cardio.bed
###### contains 30 unique cell types for kidney.bed
###### contains 22 unique cell types for prostate.bed

-------------------------------------------------------------------------------------------------------

#### Peak Annotation and Visualisation 
-----------------------------------------------------------
  
##### loading ChIPseeker and associated packages

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
###### This draws on an annotation databases generated from UCSC by exposing these as TxDb objects
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
###### Limit this object to chromosome 2
seqlevels(txdb) <- "chr2"
library(clusterProfiler)

###### create a files object containing all 4 bed files/ tissue types
files = c('breastTF.bed', 'cardioTF.bed', 'kidneyTF.bed', 'prostateTF.bed')
###### make a list of these files 
listfiles <- as.list(files)
###### name the values of the list
names(listfiles) <- c("breast", "cardio", "kidney", "prostate")

###### annotate the peak data 5000 base pairs upstream and 100 bps downstream (slightly into the 5'UTR) of the TSS for each file in samplefiles
peakAnnoList <- lapply(listfiles, annotatePeak, TxDb=txdb,tssRegion=c(-5000, 5000), level = "gene", verbose=FALSE)


###### plot feature distribution of peaks according to tissue type
plotAnnoPie(peakAnnoList[["breast"]], col = blues9)
plotAnnoPie(peakAnnoList[["cardio"]], col = blues9)
plotAnnoPie(peakAnnoList[["kidney"]], col = blues9)
plotAnnoPie(peakAnnoList[["prostate"]], col = blues9)

#### or plot altogether
plotAnnoBar(peakAnnoList, title = "Fetaure distribution according to tissue type")

###### plot the location of each peak relative to the TSS
###### the x-axis gives the percentage of sites, while the color represents the distance from the TSS
plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci relative to TSS and according to tissue type")

plot <- lapply(peakAnnoList, upsetplot, vennpie = TRUE)

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


###### read these files back into R objects of the same name for future use
breast_annot <- read.table("breast_annotation.txt", sep = "\t")
cardio_annot <- read.table("cardio_annotation.txt", sep = "\t")
kidney_annot <- read.table("kidney_annotation.txt", sep = "\t")
prostate_annot <- read.table("prostate_annotation.txt", sep = "\t")

###### assign names to the variables in the annotation files
colnames(breast_annot) <- c("chr", "start", "end", "width", "strand", "MACSqValue", "SRX_ID", "TF", "cell_type", "annotation", "geneCh", "geneS", "geneE", "geneLen", "geneSt", "geneID", "distanceToTSS")  
colnames(cardio_annot) <- c("chr", "start", "end", "width", "strand", "MACSqValue", "SRX_ID", "TF", "cell_type", "annotation", "geneCh", "geneS", "geneE", "geneLen", "geneSt", "geneID", "distanceToTSS") 
colnames(kidney_annot) <- c("chr", "start", "end", "width", "strand", "MACSqValue", "SRX_ID", "TF", "cell_type", "annotation", "geneCh", "geneS", "geneE", "geneLen", "geneSt", "geneID", "distanceToTSS")  
colnames(prostate_annot) <- c("chr", "start", "end", "width", "strand", "MACSqValue", "SRX_ID", "TF", "cell_type", "annotation", "geneCh", "geneS", "geneE", "geneLen", "geneSt", "geneID", "distanceToTSS") 

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###### Filtering Based on Annotation
--------------------------------------------

library(ChIPseeker)
library(dplyr)

###### filter the data from the annotation to include only the TFs binding to the Promoter region and 
###### Distal Intergenic regions upstream of PXDN
######. The aim of this is to narrow down the cis acting TFs 
breast_annot_filt <- breast_annot %>% filter(annotation == c("Promoter (4-5kb)", "Promoter (3-4kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", "Distal Intergenic"))
###### Change the class of the MACSqValue and check it
breast_annot_filt$MACSqValue <- as.numeric(breast_annot_filt$MACSqValue)
class(breast_annot_filt$MACSqValue)
###### filter according to the MACSqValue: only select experiments with a -log10(MACSq) greater than 150
breast_annot_filt <- breast_annot_filt %>% filter(MACSqValue > 200) %>% filter(TF != "Epitope tags")


###### repeat these steps for the other files:
cardio_annot_filt <- cardio_annot %>% filter(annotation == c("Promoter (4-5kb)", "Promoter (3-4kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", "Distal Intergenic"))
cardio_annot_filt$MACSqValue <- as.numeric(cardio_annot_filt$MACSqValue)
class(cardio_annot_filt$MACSqValue)
cardio_annot_filt <- cardio_annot_filt %>% filter(MACSqValue > 200) %>% filter(TF != "Epitope tags")


kidney_annot_filt <- kidney_annot %>% filter(annotation == c("Promoter (4-5kb)", "Promoter (3-4kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", "Distal Intergenic"))
kidney_annot_filt$MACSqValue <- as.numeric(kidney_annot_filt$MACSqValue)
kidney_annot_filt <- kidney_annot_filt %>% filter(MACSqValue > 200) %>% filter(TF != "Epitope tags")

prostate_annot_filt <- prostate_annot %>% filter(annotation == c("Promoter (4-5kb)", "Promoter (3-4kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", "Distal Intergenic"))
prostate_annot_filt$MACSqValue <- as.numeric(prostate_annot_filt$MACSqValue)
prostate_annot_filt <- prostate_annot_filt %>% filter(MACSqValue > 200) %>% filter(TF != "Epitope tags")


###### write these data frames into files 
write.table(breast_annot_filt, file="breast_annot_filt.txt", sep = "\t", row.names = FALSE)
write.table(cardio_annot_filt, file="cardio_annot_filt.txt", sep = "\t", row.names = FALSE)
write.table(kidney_annot_filt, file="kidney_annot_filt.txt", sep = "\t", row.names = FALSE)
write.table(prostate_annot_filt, file="prostate_annot_filt.txt", sep = "\t", row.names = FALSE)

###### read these files back into R objects of the same name for future use
breast_annot_filt <- read.table("breast_annot_filt.txt", sep = "\t", skip = 1)
cardio_annot_filt <- read.table("cardio_annot_filt.txt", sep = "\t", skip = 1)
kidney_annot_filt <- read.table("kidney_annot_filt.txt", sep = "\t", skip = 1)
prostate_annot_filt <- read.table("prostate_annot_filt.txt", sep = "\t", skip = 1)


#### rerun the colnames code to assign variables names
#### note that the gene ID 7837 does correspond to PXDN

colnames(breast_annot_filt) <- c("chr", "start", "end", "readL", "s", "MACSq", "SRX_ID", "TF", "cell_type", "annotation", "geneChr", "geneS", "geneE", "geneLen", "geneSt", "geneID", "distanceToTSS")  
colnames(cardio_annot_filt) <- c("chr", "start", "end", "readL", "s", "MACSq", "SRX_ID", "TF", "cell_type", "annotation", "geneChr", "geneS", "geneE", "geneLen", "geneSt", "geneID", "distanceToTSS") 
colnames(kidney_annot_filt) <- c("chr", "start", "end", "readL", "s", "MACSq", "SRX_ID", "TF", "cell_type", "annotation", "geneChr", "geneS", "geneE", "geneLen", "geneSt", "geneID", "distanceToTSS")  
colnames(prostate_annot_filt) <- c("chr", "start", "end", "readL", "s", "MACSq", "SRX_ID", "TF", "cell_type", "annotation", "geneChr", "geneS", "geneE", "geneL", "geneSt", "geneID", "distanceToTSS") 

####### make a list of the filtered files and name them
filt_annot_files = c('breast_annot_filt.txt', 'cardio_annot_filt.txt', 'kidney_annot_filt.txt', 'prostate_annot_filt.txt')
###### make a list of these files 
filt_annot_files_list <- as.list(filt_annot_files)
###### name the values of the list
names(filt_annot_files_list) <- c("breastFilt", "cardioFilt", "kidneyFilt", "prostateFilt")

-----------------------------------------------------------------------------------------------------------------------------------------------

### Plotting the filtered TFs according to significance and prevalence per tissue
--------------------------------------------------------------------------------------

##### Use ggplot to plot a bar graph of the TFs and their prevalence in each tissue type according to MACSq score. The most significant TFs can be analysed further.
  
  
library(ggplot2)

par(mfrow=c(4,1))

breast_annot_filt$n<-1
TFbreast <- breast_annot_filt %>% arrange(MACSq) %>% ggplot(aes(x = TF, n, fill = MACSq)) + geom_bar(stat = "identity") + scale_fill_continuous(low="blue", high="red")
TFbreast + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Bar graph of the prevalence of transcription factors \nregulating PXDN in breast tissue") + labs(y = "count", x = "Transcription Factors")

cardio_annot_filt$n<-1
TFcardio <- cardio_annot_filt %>% arrange(MACSq) %>% ggplot(aes(x = TF, n, fill = MACSq)) + geom_bar(stat = "identity") + scale_fill_continuous(low="blue", high="red")
TFcardio + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Bar graph of the prevalence of transcription factors \nregulating PXDN in cardiovascular tissue") + labs(y = "count", x = "Transcription Factors")


kidney_annot_filt$n<-1
TFkidney <- kidney_annot_filt %>% arrange(MACSq) %>% ggplot(aes(x = TF, n, fill = MACSq)) + geom_bar(stat = "identity") + scale_fill_continuous(low="blue", high="red")
TFkidney + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Bar graph of the prevalence of transcription factors \nregulating PXDN in kidney tissue") + labs(y = "count", x = "Transcription Factors")


prostate_annot_filt$n<-1
TFprostate <- prostate_annot_filt %>% arrange(MACSq) %>% ggplot(aes(x = TF, n, fill = MACSq)) + geom_bar(stat = "identity") + scale_fill_continuous(low="blue", high="red")
TFprostate + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Bar graph of the prevalence of transcription factors \nregulating PXDN in prostate tissue") + labs(y = "count", x = "Transcription Factors")

###### Use the bar graphs to narrrow down the final TFs per tissue

#### select the chosen TFs from the file 
FinalTFCardio = cardio_annot_filt[-c(14,48,49,50,52,56,57,58,80),]

FinalTFProstate <- prostate_annot_filt[-c(1,3,6,7,38,58,61,70,85,86,88,89,92,95,97,99,103,111,112 ),]

FinalTFKidney <- kidney_annot_filt[-c(1,2,4,5,10,12,14,18,19,21,22,23,25,27,30,31,32,33,40,41,42),]

FinalTFBreast <- breast_annot_filt[-c(20, 47, 53, 54, 105, 109, 114, 115, 144, 147,150, 152, 153, 155, 159, 161, 162, 164, 166, 168, 170,171, 172, 178, 181, 187, 188, 189, 196, 201, 203, 206,
                                      +                                       216, 219, 222, 231, 243, 244, 246, 249, 258, 286, 326, 328, 330, 331, 369, 375, 384 ), ]

#### Write into a file (excel and/or text)

library("writexl")
write.table(FinalTFCardio, file="FinalTFCardio.txt", sep = "\t", row.names = FALSE)
write_xlsx(FinalTFCardio, "FinalTfCardio.xlsx")

write.table(FinalTFProstate, file="FinalTFProstate.txt", sep = "\t", row.names = FALSE)
write_xlsx(FinalTFProstate, "FinalTfProstate.xlsx")

write.table(FinalTFKidney, file="FinalTFKidney.txt", sep = "\t", row.names = FALSE)
write_xlsx(FinalTFKidney, "FinalTfKidney.xlsx")

write.table(FinalTFBreast, file="FinalTFBreast.txt", sep = "\t", row.names = FALSE)
write_xlsx(FinalTFBreast, "FinalTfBreast.xlsx")

### Make list of files

FinalTFs = c('FinalTFBreast.txt', 'FinalTFCardio.txt', 'FinalTFKidney.txt', 'FinalTFProstate.txt')
FinalTFsList <- as.list(FinalTFs)


-----------------------------------------------------------------------------------------------------------------------------

####### Analyzing and Visualizing the Filtered Data in ChIPseeker
---------------------------------------------------------------------
library(GenomicRanges)
library(ChIPseeker)

###### Read peaks into GRanges R objects
breast_peak <- readPeakFile(FinalTFsList[[1]])
cardio_peak <- readPeakFile(FinalTFsList[[2]])
kidney_peak <- readPeakFile(FinalTFsList[[3]]) 
prostate_peak <- readPeakFile(FinalTFsList[[4]])

  
###### use the MACS score to find the peak locations in PXDN and visualize this for each tissue

covplot(breast_peak, weightCol="MACSq", title = "ChIP-Seq leak locations for breast tissue \naccording to binding significance")
covplot(cardio_peak, weightCol= "MACSq",title = "ChIP-Seq peak locations for cardiovascular tissue \naccording to binding significance")
covplot(kidney_peak, weightCol= "MACSq", title = "ChIP-Seq peak locations for kidney tissue \naccording to binding significance")
covplot(prostate_peak, weightCol= "MACSq", title = "ChIP-Seq peak locations for prostate tissue \naccording to binding significance")

------------------------------------------------------------------------------------------------------------------------------

###### Functional Enrichment Analysis for PXDN
-------------------------------------
  
library(ReactomePA)

##### functionally annotate the PXDN gene pathways. Only one tissue is necessary since they
##### all contain the same gene ID (PXDN). This does not provide any new information but acts as confirmation 
##### of annotation according to PXDN and is interesting to plot.
pathway1 <- enrichPathway(as.data.frame(breast_annot_filt)$geneID)
head(pathway1,2)
##### visualize this annotation 
dotplot(pathway1)

-------------------------------------------------------------------------------------------------------------------

##### Getting the FASTA sequences of the TFBS
-------------------------------------------------------------------------
##### TFs were chosen based on prevalence, MACSq significance values and biological function in the various tissue types
#### The TFs are filtered out separately to better be able to define the binding region overlaps from the 
#### different experiments- in this way each binding site sequence is retrieved only once from the ENSEMBL database per tissue
  

###  Breast Tissue TFs:
---------------------------------- 
ESR1Breast <- filter(FinalTFBreast, TF == "ESR1") %>% select(c(chr, start, end, TF, MACSq))
write.table(ESR1Breast, "ESR1Breast.BED" , sep = "\t", col.names = FALSE, row.names = FALSE)

CTCFBreast <- filter(FinalTFBreast, TF == "CTCF")%>% select(c(chr, start, end, TF, MACSq))
CDKN1Breast <- filter(FinalTFBreast, TF == "CDKN1B")%>% select(c(chr, start, end, TF, MACSq))
FOSBreast <- filter(FinalTFBreast, TF == "FOS")%>% select(c(chr, start, end, TF, MACSq))
FOXA1Breast <- filter(FinalTFBreast, TF == "FOXA1") %>% select(c(chr, start, end, TF, MACSq))
GATA3Breast <- filter(FinalTFBreast, TF == "GATA3")%>% select(c(chr, start, end, TF, MACSq))
GRHC2Breast <- filter(FinalTFBreast, TF == "GRHC2")%>% select(c(chr, start, end, TF, MACSq))
HIF1ABreast <- filter(FinalTFBreast, TF == "HIF1A") %>% select(c(chr, start, end, TF, MACSq))
HSF1Breast <- filter(FinalTFBreast, TF == "HSF1")%>% select(c(chr, start, end, TF, MACSq))
JUNBreast <- filter(FinalTFBreast, TF == "JUN")%>% select(c(chr, start, end, TF, MACSq))
NRF1Breast <- filter(FinalTFBreast, TF == "NRF1")%>% select(c(chr, start, end, TF, MACSq))
STAG2Breast <- filter(FinalTFBreast, TF == "STAG2")%>% select(c(chr, start, end, TF, MACSq))
SMC1ABreast <- filter(FinalTFBreast, TF == "SMC1A")%>% select(c(chr, start, end, TF, MACSq))
STAT3Breast <- filter(FinalTFBreast, TF == "STAT3")%>% select(c(chr, start, end, TF, MACSq))
TP53Breast <- filter(FinalTFBreast, TF == "TP53")%>% select(c(chr, start, end, TF, MACSq))
YAP1Breast <- filter(FinalTFBreast, TF == "YAP1")%>% select(c(chr, start, end, TF, MACSq))

###### All FASTA sequences obtained from the ENSEMBL database
--------------------------------------------------------------
library(httr)
library(jsonlite)
library(xml2)
server <- "https://rest.ensembl.org"

#### The genomic coordinates from the BED files can then be specified and the number of bases upstream and 
#### downstream of the binding site. In this case each read is set to 500bps in total since
#### this is the optimum length stated by MEME for de novo motif finding in downstream analysis of these sequences.

#### specify the region of interest 
ext <- "/sequence/region/human/2:1595689..1596091:-1?expand_3prime=44;expand_5prime=44"

#### Retrieve the genomic sequences in FASTA format
r <- GET(paste(server, ext, sep = ""), content_type("text/x-fasta"))

#### View the content/ sequence
print(content(r))

#### Each sequence is copied and saved into a file according to tissue type and TF. The list of these
#### sequences are used for MAST and MEME analysis.

##### Repeat this process for the other tissue types and their TFs:



##### This is a function used to determine how many bases upstream and downstream of the the binding site
##### need to be added for the read to be 500 bps in length (for MEME-ChIP)

funky <- function (start, end){
  diff <- 500 - (end - start)
  out <- diff/2
  print (out)
}

### Cardiovascular Tissue TFs:
---------------------------------
EP300Cardio <- filter(FinalTFCardio, TF == "EP300")%>% select(c(chr, start, end, TF, MACSq, s))

CTCFCardio <- filter(FinalTFCardio, TF == "CTCF")
FOSCardio <- filter(FinalTFCardio, TF == "FOS")
JUNDCardio <- filter(FinalTFCardio, TF == "JUND")
RAD21Cardio <- filter(FinalTFCardio, TF == "RAD21")
RELACardio <- filter(FinalTFCardio, TF == "RELA")
STAG2Cardio <- filter(FinalTFCardio, TF == "STAG2")
TAL1Cardio <- filter(FinalTFCardio, TF == "TAL1")
TCF21Cardio <- filter(FinalTFCardio, TF == "TCF21")



#### Kidney Tissue TFs
RELAkidney <- filter(FinalTFKidney, TF == "RELA")
RAD21kidney <- filter(FinalTFKidney, TF == "RAD21")
NR3C1kidney <- filter(FinalTFKidney, TF == "NR3C1")
CTCFkidney <- filter(FinalTFKidney, TF == "CTCF")
PAX8kidney <- filter(FinalTFKidney, TF == "PAX8")
FOXA1kidney <- filter(FinalTFKidney, TF == "FOXA1")
DPF2kidney <- filter(FinalTFKidney, TF == "DPF2")

#### Prostate Tissue TFs

ARprostate <- filter(FinalTFProstate, TF == "AR")
CTCFprostate <- filter(FinalTFProstate, TF == "CTCF")
FOXA1prostate <- filter(FinalTFProstate, TF == "FOXA1")
HOXB13prostate <- filter(FinalTFProstate, TF == "HOXB13")
TCF7L2prostate <- filter(FinalTFProstate, TF == "TCF7L2")

  
------------------------------------------------------------------------------

#### Binding matrices downloaded in MEME format from JASPER and FIMO (web-based) from
#### MEME-suite used to scan for binding sites in the FASTA sequences

  -------------------------------------------------------------------------

#### Retrieving the JASPER TF binding matrices
  
library(JASPAR2020)
library(TFBSTools)

##### retrieve JASPAR matrices 

CTCF <- getMatrixByName(JASPAR2020, name="CTCF")
FOS <- getMatrixByName(JASPAR2020, name="FOS")
FOXA1 <- getMatrixByName(JASPAR2020, name="FOXA1")
ESR1 <- getMatrixByName(JASPAR2020, name="ESR1")
HIF1A <- getMatrixByName(JASPAR2020, name="HIF1A")
HSF1 <- getMatrixByName(JASPAR2020, name="HSF1")
JUN <- getMatrixByName(JASPAR2020, name="JUN")
NRF1 <- getMatrixByName(JASPAR2020, name="NRF1")
STAT3 <- getMatrixByName(JASPAR2020, name="STAT3")
TP53 <- getMatrixByName(JASPAR2020, name="TP53")

JUND <- getMatrixByName(JASPAR2020, name="JUND")
RELA <- getMatrixByName(JASPAR2020, name="RELA")
TAL1 <- getMatrixByName(JASPAR2020, name="TAL1")
TCF21 <- getMatrixByName(JASPAR2020, name="TCF21")

NR3C1 <- getMatrixByName(JASPAR2020, name="NR3C1")
PAX8 <- getMatrixByName(JASPAR2020, name="PAX8")

AR <- getMatrixByName(JASPAR2020, name="AR")
HOXB13 <- getMatrixByName(JASPAR2020, name="HOXB13")
TCF7L2 <- getMatrixByName(JASPAR2020, name="TCF7L1")

------------------------------------------------------------------------------------------
  
  
  













    

