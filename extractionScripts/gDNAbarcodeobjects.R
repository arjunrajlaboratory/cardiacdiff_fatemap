#Genomic DNA barcode extraction script for paper
#try different cutoffs
#July2021

library(tidyverse)

dataDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/extractedData/gDNAbarcodeobjects/"

#location of barcode datsets
gDNAbarDirectory <- "/Users/conniejiang/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/rawData/gDNAseq/" 

directories <- c('gDNA1_MOI_0.1/','gDNA1_MOI_0.5/','gDNA1_MOI_1.0/',
                 'gDNA2_MOI_0.1/','gDNA2_MOI_0.23/','gDNA2_MOI_1.0/')
datanames <- c('sampleTable_g1_MOI0.1','sampleTable_g1_MOI0.5','sampleTable_g1_MOI1.0',
              'sampleTable_g2_MOI0.1','sampleTable_g2_MOI0.23','sampleTable_g2_MOI1.0')
datasets <- list()

for (i in 1:length(datanames)){
  datasets[[i]] = as_tibble(read.table(paste0(gDNAbarDirectory,directories[i], 'stepThreeStarcodeShavedReads.txt'), stringsAsFactors=F, header = T)) %>%
    select(-cellID) %>%
    transmute(countsperUMI = UMI, OriginalBarcode, BC50StarcodeD8, BC40StarcodeD8, BC30StarcodeD8, BC30StarcodeD6, SampleNum) %>%
    filter(!str_detect(OriginalBarcode, "NNN"))
  saveRDS(datasets[[i]], file = paste0(dataDirectory,datanames[i],'.rds'))
}
