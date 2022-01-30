#10X Barcode extraction script for paper
#try different cutoffs
#July2021

library(Seurat)
library(tidyverse)
library(VennDiagram)

dataDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/extractedData/createbarcodeobj/"

#location of barcode reads
home2Directory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/rawData/10XData/" 
corebarcodes = as_tibble(read.table(paste0(home2Directory, 'stepThreeStarcodeShavedReads.txt'), stringsAsFactors=F, header = T)) 
corebarcodes = corebarcodes %>% na.omit() %>%
  filter(SampleNum %in% c("1","2"))#remove empty rows --> 6076685 rows
corebarcodes$SampleNum = sub("1", "sibA", corebarcodes$SampleNum)
corebarcodes$SampleNum = sub("2", "sibB", corebarcodes$SampleNum)

corebarcodes30d6 = corebarcodes %>%
  dplyr::select(cellID,UMI,BC30StarcodeD6,SampleNum) %>%
  unique() #remove PCR duplicates --> 781999

umiCut = 2 # Minimum UMI cutoff (per barcode) for reliable analysis. NOTE that unless this is taking place after unique() subset of barcodes it is not actually an umiCut but rather cutoff for total number of reads (without accounting for PCR duplicates)
fracumiCut = 0.3
#to a point, increasing fracumiCut decreases total cells recovered but increases proportion with 1 barcode

linCut = 1 # remove all the cases where one CellID is associated with X(one/two etc) barcodes.

coreupToLineageCounting30d6 = corebarcodes30d6 %>% 
  filter(SampleNum %in% c("sibA","sibB")) %>% #redundant
  group_by(cellID, BC30StarcodeD6, SampleNum) %>%
  summarise(nUMI = length(SampleNum)) %>%
  filter(nUMI >= umiCut) %>% 
  group_by(cellID,SampleNum) %>%
  mutate(fracUMI = nUMI/sum(nUMI)) %>%
  filter(fracUMI >= fracumiCut) %>%
  group_by(cellID, SampleNum) %>%
  mutate(nLineages = length(cellID)) 
#6952 for umiCut=2+fracumiCut=0.3

###taking only single cellid-->barcode mappings; removes doublets, ambient RNA (more than 1 barcode per cellID)
corelinCountTooverlaps30d6 = coreupToLineageCounting30d6 %>%
  ungroup() %>%
  filter(nLineages <= linCut) %>%
  unique() %>%
  mutate(cellIDSN = paste(SampleNum, cellID, sep="_"))
#5860 for umiCut=2+fracumiCut=0.3

#save RDS:
saveRDS(corelinCountTooverlaps30d6, file = paste0(dataDirectory,'corelinCountTooverlaps30d6.rds'))

corelinCountTooverlaps30d6 = readRDS(file = paste0(dataDirectory,'corelinCountTooverlaps30d6.rds'))
