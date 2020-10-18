# Honours Project
### ChIP-Seq analysis workflow and code 

This project aims to identify potential transcription factors (TFs) involved in _PXDN_ regulation using publically available ChIP-Seq data. _PXDN_ is a gene that encodes the peroxidasin protein, which has been linked to several prominant diseases worldwide. These include breast and prostate cancers as well as kidney fibrosis and cardiovascular disease. Thus, these four tissue types are used in this analysis - breast, prostate, kidney and cardiovascular - to better understand the regulation of _PXDN_ in each. 

The raw data for this project was sourced from the ChIP-Atlas database, with the help of Prof. Shinya Oki - the co-author of these tools @inutano/chip-atlas. The ChIP-Seq reads were collected from major projects such as ENCODE as well as from smaller projects stored in the SRA. The reads were aligned to the human reference genome with Bowtie2 and MACS2 was used for peak calling of the significantly enriched aligned reads (link to ChIP-Atlas publication: https://www.embopress.org/doi/full/10.15252/embr.201846255).

These results were downloaded in BED file format for the region of _PXDN_ (hg19::chr2:1,577,805-1,806,141) for each tissue at a q-value greater than 1x10-5. Each of these files are available in the the TF folder. 

This script includes the steps of the downstream analysis of this data using a number of R packages from Bioconductor (ChIPseeker in particular). This includes the functional annotation, further filtering and visualisation of the data. The results generated here were used to find the known and _de novo_ binding sites for the final TFs using the MEME Suite command line tool (not included in this script).

This is my first time conducting bioinformatic analysis in R so all comments and suggestions are most welcome.
