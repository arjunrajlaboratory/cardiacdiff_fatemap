#Supplemental figure: gDNA barcode sequencing analysis/venn diagrams
#July2021

library(tidyverse)
library(VennDiagram)
survive <- function(datarow){ #reads in a row from sampleTable_g2_MOI1.0_30d8_iPS, designates binomial split by number of UMIs
  barcode_binom = rbinom(datarow[3],1,0.5)
}

#pull data from extract directory 
extractDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/extractedData/gDNAbarcodeobjects/"
#pull datasets not featured in Fig 3
sampleTable_g1_MOI0.1 <- readRDS(paste0(extractDirectory,'sampleTable_g1_MOI0.1.rds'))
sampleTable_g1_MOI1.0 <- readRDS(paste0(extractDirectory,'sampleTable_g1_MOI1.0.rds'))
sampleTable_g2_MOI0.1 <- readRDS(paste0(extractDirectory,'sampleTable_g2_MOI0.1.rds'))
sampleTable_g2_MOI0.23 <- readRDS(paste0(extractDirectory,'sampleTable_g2_MOI0.23.rds'))
sampleTable_g2_MOI1.0 <- readRDS(paste0(extractDirectory,'sampleTable_g2_MOI1.0.rds'))

#plotData directory
plotDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotData/Supplement/gDNA_VD/"

#g1_MOI0.1
#split barcodes into iPS and differentiated objects
sampleTable_g1_MOI0.1_30d8_diff = sampleTable_g1_MOI0.1 %>%
  filter(SampleNum %in% c("diff_01A_g1","diff_01B_g1","repdiff_01B_IC_g1")) %>% #consider altering this to collapse diff 0.1B normal and IC samples
  dplyr::select(countsperUMI,BC30StarcodeD8,SampleNum) %>% #33532 UMIs for A, 29416 for B using 30d8
  group_by(BC30StarcodeD8, SampleNum) %>%
  summarise(UMI = n(), counts = sum(countsperUMI))#collapse and count how many distinct UMIs per barcode
#32055

sampleTable_g1_MOI0.1_30d8_iPS = sampleTable_g1_MOI0.1 %>%
  filter(SampleNum %in% c("iPS_01_g1")) %>%
  dplyr::select(countsperUMI,BC30StarcodeD8,SampleNum) %>% #3901070 UMIs using 30d8
  group_by(BC30StarcodeD8, SampleNum) %>%
  summarise(UMI = n(), counts = sum(countsperUMI)) #collapse and count how many distinct UMIs per barcode
#59679 unique barcodes

#determine observed overlap
splitAbar_g1_MOI0.1 = sampleTable_g1_MOI0.1_30d8_diff %>%
  filter(SampleNum == "diff_01A_g1") %>%
  dplyr::select(BC30StarcodeD8) #15260 barcodes
splitBbar_g1_MOI0.1 = sampleTable_g1_MOI0.1_30d8_diff %>%
  filter(SampleNum %in% c("diff_01B_g1","repdiff_01B_IC_g1")) %>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #14119 barcodes; more than expected
split_initialbar_g1_MOI0.1 = sampleTable_g1_MOI0.1_30d8_iPS %>%
  filter(SampleNum == "iPS_01_g1") %>% #redundant
  dplyr::select(BC30StarcodeD8) #59679 barcodes 

venn.diagram(
  x = list(unlist(splitAbar_g1_MOI0.1), unlist(splitBbar_g1_MOI0.1)),
  category.names = c("splitA" , "splitB"),
  filename = paste0(plotDirectory,'g1MOI0.1/nocutoff/g1MOI0.1_VDdiffs.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.015, 0.015),
  cat.fontfamily = "sans"
) #more overlap than for MOI 0.1??
venn.diagram(
  x = list(unlist(split_initialbar_g1_MOI0.1), unlist(splitAbar_g1_MOI0.1)),
  category.names = c("iPS" , "splitA"),
  filename = paste0(plotDirectory,'g1MOI0.1/nocutoff/g1MOI0.1_VDsplitAiPS.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.015, 0.015),
  cat.fontfamily = "sans"
)
venn.diagram(
  x = list(unlist(split_initialbar_g1_MOI0.1), unlist(splitBbar_g1_MOI0.1)),
  category.names = c("iPS" , "splitB"),
  filename = paste0(plotDirectory,'g1MOI0.1/nocutoff/g1MOI0.1_VDsplitBiPS.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.015, 0.015),
  cat.fontfamily = "sans"
)

#simulated overlap
alliPS_g1_MOI0.1 = sampleTable_g1_MOI0.1_30d8_iPS %>%
  filter(SampleNum %in% c("iPS_01_g1")) %>%
  dplyr::select(BC30StarcodeD8) 
#for each barcode, split UMIs across 2 splits -- need only run this once
splitA_g1_MOI0.1 = vector()
splitB_g1_MOI0.1 = vector()
set.seed(2059)
binomsplit_g1_MOI0.1 = apply(sampleTable_g1_MOI0.1_30d8_iPS,1,survive)
for (i in 1:nrow(sampleTable_g1_MOI0.1_30d8_iPS)){
  barcode = sampleTable_g1_MOI0.1_30d8_iPS$BC30StarcodeD8[i]
  splitA_g1_MOI0.1 = c(splitA_g1_MOI0.1, rep(barcode,sum(binomsplit_g1_MOI0.1[[i]]==0)))
  splitB_g1_MOI0.1 = c(splitB_g1_MOI0.1, rep(barcode,sum(binomsplit_g1_MOI0.1[[i]]==1)))
} 
saveRDS(splitA_g1_MOI0.1,file=paste0(plotDirectory,'g1MOI0.1/simulation/splitA_g1_MOI0.1.rds'))
saveRDS(splitB_g1_MOI0.1,file=paste0(plotDirectory,'g1MOI0.1/simulation/splitB_g1_MOI0.1.rds'))
#skip to here to read in RDS files
splitA_g1_MOI0.1 = readRDS(file=paste0(plotDirectory,'g1MOI0.1/simulation/splitA_g1_MOI0.1.rds'))
splitB_g1_MOI0.1 = readRDS(file=paste0(plotDirectory,'g1MOI0.1/simulation/splitB_g1_MOI0.1.rds'))
set.seed(2059)
#factor in loss of some UMIs: for g1MOI0.1 survivability = 0.015
survivability_g1_MOI0.1 = 0.015
survivalA_g1_MOI0.1 <- sample(length(splitA_g1_MOI0.1),length(splitA_g1_MOI0.1)*survivability_g1_MOI0.1)
splitAsurvive_g1_MOI0.1 <- splitA_g1_MOI0.1[survivalA_g1_MOI0.1]

survivalB_g1_MOI0.1 <- sample(length(splitB_g1_MOI0.1),length(splitB_g1_MOI0.1)*survivability_g1_MOI0.1)
splitBsurvive_g1_MOI0.1 <- splitB_g1_MOI0.1[survivalB_g1_MOI0.1]

venn.diagram(
  x = list(unlist(alliPS_g1_MOI0.1), splitAsurvive_g1_MOI0.1),
  category.names = c("iPS" , "splitAsim"),
  filename = paste0(plotDirectory,'g1MOI0.1/simulation/splitAiPS_survivability',survivability_g1_MOI0.1,'.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)
venn.diagram(
  x = list(unlist(alliPS_g1_MOI0.1), splitBsurvive_g1_MOI0.1),
  category.names = c("iPS" , "splitBsim"),
  filename = paste0(plotDirectory,'g1MOI0.1/simulation/splitBiPS_survivability',survivability_g1_MOI0.1,'.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)

splitAbarcodes_g1_MOI0.1 = unique(splitAsurvive_g1_MOI0.1)
splitBbarcodes_g1_MOI0.1 = unique(splitBsurvive_g1_MOI0.1)

venn.diagram(
  x = list(splitAsurvive_g1_MOI0.1, splitBsurvive_g1_MOI0.1),
  category.names = c("splitAsim" , "splitBsim"),
  filename = paste0(plotDirectory,'g1MOI0.1/simulation/splitA_splitBoverlap',survivability_g1_MOI0.1,'.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-15, 15),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)

#g1_MOI1.0
#split barcodes into iPS and differentiated objects
sampleTable_g1_MOI1.0_30d8_diff = sampleTable_g1_MOI1.0 %>%
  filter(SampleNum %in% c("diff_10AGFPpos_g1","diff_10AGFPneg_g1","diff_10BGFPpos_g1","diff_10BGFPneg_g1")) %>%
  dplyr::select(countsperUMI,BC30StarcodeD8,SampleNum) %>% #2797919 UMIs for A GFP+, 2481109 for A GFP-, 2934332 for B GFP+, 3383230 for B GFP- using 30d8
  group_by(BC30StarcodeD8, SampleNum) %>%
  summarise(UMI = n(), counts = sum(countsperUMI)) #collapse and count how many distinct UMIs per barcode
#533297

sampleTable_g1_MOI1.0_30d8_iPS = sampleTable_g1_MOI1.0 %>%
  filter(SampleNum %in% c("iPS_10_g1")) %>%
  dplyr::select(countsperUMI,BC30StarcodeD8,SampleNum) %>% #12688701 UMIs using 30d8
  group_by(BC30StarcodeD8, SampleNum) %>%
  summarise(UMI = n(), counts = sum(countsperUMI)) #collapse and count how many distinct UMIs per barcode
#247405 barcodes

#determine observed overlap
splitAbar_g1_MOI1.0 = sampleTable_g1_MOI1.0_30d8_diff %>%
  filter(SampleNum %in% c("diff_10AGFPpos_g1","diff_10AGFPneg_g1")) %>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 168149 barcodes
splitBbar_g1_MOI1.0 = sampleTable_g1_MOI1.0_30d8_diff %>%
  filter(SampleNum %in% c("diff_10BGFPpos_g1","diff_10BGFPneg_g1")) %>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 173237 barcodes
split_initialbar_g1_MOI1.0 = sampleTable_g1_MOI1.0_30d8_iPS %>%
  filter(SampleNum == "iPS_10_g1")%>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 247405 barcodes

venn.diagram(
  x = list(unlist(splitAbar_g1_MOI1.0), unlist(splitBbar_g1_MOI1.0)),
  category.names = c("splitA" , "splitB"),
  filename = paste0(plotDirectory,'g1MOI1.0/nocutoff/g1MOI1.0_VDdiffs.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.015, 0.015),
  cat.fontfamily = "sans"
) 
venn.diagram(
  x = list(unlist(split_initialbar_g1_MOI1.0), unlist(splitAbar_g1_MOI1.0)),
  category.names = c("iPS" , "splitA"),
  filename = paste0(plotDirectory,'g1MOI1.0/nocutoff/g1MOI1.0_VDsplitAiPS.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.015, 0.015),
  cat.fontfamily = "sans"
)
venn.diagram(
  x = list(unlist(split_initialbar_g1_MOI1.0), unlist(splitBbar_g1_MOI1.0)),
  category.names = c("iPS" , "splitB"),
  filename = paste0(plotDirectory,'g1MOI1.0/nocutoff/g1MOI1.0_VDsplitBiPS.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.015, 0.015),
  cat.fontfamily = "sans"
)

#simulated overlap
alliPS_g1_MOI1.0 = sampleTable_g1_MOI1.0_30d8_iPS %>%
  filter(SampleNum %in% c("iPS_01_g1")) %>%
  dplyr::select(BC30StarcodeD8) 
#for each barcode, split UMIs across 2 splits -- need only run this once
splitA_g1_MOI1.0 = vector()
splitB_g1_MOI1.0 = vector()
set.seed(2059)
binomsplit_g1_MOI1.0 = apply(sampleTable_g1_MOI1.0_30d8_iPS,1,survive)
for (i in 1:nrow(sampleTable_g1_MOI1.0_30d8_iPS)){
  barcode = sampleTable_g1_MOI1.0_30d8_iPS$BC30StarcodeD8[i]
  splitA_g1_MOI1.0 = c(splitA_g1_MOI1.0, rep(barcode,sum(binomsplit_g1_MOI1.0[[i]]==0)))
  splitB_g1_MOI1.0 = c(splitB_g1_MOI1.0, rep(barcode,sum(binomsplit_g1_MOI1.0[[i]]==1)))
} 
saveRDS(splitA_g1_MOI1.0,file=paste0(plotDirectory,'g1MOI1.0/simulation/splitA_g1_MOI1.0.rds'))
saveRDS(splitB_g1_MOI1.0,file=paste0(plotDirectory,'g1MOI1.0/simulation/splitB_g1_MOI1.0.rds'))
#skip to here to read in RDS files
splitA_g1_MOI1.0 = readRDS(file=paste0(plotDirectory,'g1MOI1.0/simulation/splitA_g1_MOI1.0.rds'))
splitB_g1_MOI1.0 = readRDS(file=paste0(plotDirectory,'g1MOI1.0/simulation/splitB_g1_MOI1.0.rds'))
set.seed(2059)
#factor in loss of some UMIs: for g1MOI1.0 survivability = 0.18
survivability_g1_MOI1.0 = 0.18
survivalA_g1_MOI1.0 <- sample(length(splitA_g1_MOI1.0),length(splitA_g1_MOI1.0)*survivability_g1_MOI1.0)
splitAsurvive_g1_MOI1.0 <- splitA_g1_MOI1.0[survivalA_g1_MOI1.0]

survivalB_g1_MOI1.0 <- sample(length(splitB_g1_MOI1.0),length(splitB_g1_MOI1.0)*survivability_g1_MOI1.0)
splitBsurvive_g1_MOI1.0 <- splitB_g1_MOI1.0[survivalB_g1_MOI1.0]

venn.diagram(
  x = list(unlist(alliPS_g1_MOI1.0), splitAsurvive_g1_MOI1.0),
  category.names = c("iPS" , "splitAsim"),
  filename = paste0(plotDirectory,'g1MOI1.0/simulation/splitAiPS_survivability',survivability_g1_MOI1.0,'.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)
venn.diagram(
  x = list(unlist(alliPS_g1_MOI1.0), splitBsurvive_g1_MOI1.0),
  category.names = c("iPS" , "splitBsim"),
  filename = paste0(plotDirectory,'g1MOI1.0/simulation/splitBiPS_survivability',survivability_g1_MOI1.0,'.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)

splitAbarcodes_g1_MOI1.0 = unique(splitAsurvive_g1_MOI1.0)
splitBbarcodes_g1_MOI1.0 = unique(splitBsurvive_g1_MOI1.0)

venn.diagram(
  x = list(splitAsurvive_g1_MOI1.0, splitBsurvive_g1_MOI1.0),
  category.names = c("splitAsim" , "splitBsim"),
  filename = paste0(plotDirectory,'g1MOI1.0/simulation/splitA_splitBoverlap',survivability_g1_MOI1.0,'.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-15, 15),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)

#g2_MOI0.1
sampleTable_g2_MOI0.1_30d8_diff = sampleTable_g2_MOI0.1 %>%
  filter(SampleNum %in% c("diff_01A_g2","diff_01B_g2")) %>%
  dplyr::select(countsperUMI,BC30StarcodeD8,SampleNum) %>% #14448 UMIs for A, 13438 UMIs for B using 30d8
  group_by(BC30StarcodeD8, SampleNum) %>%
  summarise(UMI = n(), counts = sum(countsperUMI)) #collapse and count how many distinct UMIs per barcode
#16127

sampleTable_g2_MOI0.1_30d8_iPS = sampleTable_g2_MOI0.1 %>%
  filter(SampleNum %in% c("iPS_01_g2")) %>%
  dplyr::select(countsperUMI,BC30StarcodeD8,SampleNum) %>% #431323 UMIs using 30d8
  group_by(BC30StarcodeD8, SampleNum) %>%
  summarise(UMI = n(), counts = sum(countsperUMI)) #collapse and count how many distinct UMIs per barcode
#34195 

#determine observed overlap
splitAbar_g2_MOI0.1 = sampleTable_g2_MOI0.1_30d8_diff %>%
  filter(SampleNum == "diff_01A_g2") %>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 8181 barcodes
splitBbar_g2_MOI0.1 = sampleTable_g2_MOI0.1_30d8_diff %>%
  filter(SampleNum == "diff_01B_g2") %>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 7946 barcodes
split_initialbar_g2_MOI0.1 = sampleTable_g2_MOI0.1_30d8_iPS %>%
  filter(SampleNum == "iPS_01_g2")%>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 34195 barcodes

venn.diagram(
  x = list(unlist(splitAbar_g2_MOI0.1), unlist(splitBbar_g2_MOI0.1)),
  category.names = c("splitA" , "splitB"),
  filename = paste0(plotDirectory,'g2MOI0.1/nocutoff/g2MOI0.1_VDdiffs.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.015, 0.015),
  cat.fontfamily = "sans"
) #more overlap than for MOI 0.1??
venn.diagram(
  x = list(unlist(split_initialbar_g2_MOI0.1), unlist(splitAbar_g2_MOI0.1)),
  category.names = c("iPS" , "splitA"),
  filename = paste0(plotDirectory,'g2MOI0.1/nocutoff/g2MOI0.1_VDsplitAiPS.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.015, 0.015),
  cat.fontfamily = "sans"
)
venn.diagram(
  x = list(unlist(split_initialbar_g2_MOI0.1), unlist(splitBbar_g2_MOI0.1)),
  category.names = c("iPS" , "splitB"),
  filename = paste0(plotDirectory,'g2MOI0.1/nocutoff/g2MOI0.1_VDsplitBiPS.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.015, 0.015),
  cat.fontfamily = "sans"
)

#simulated overlap
alliPS_g2_MOI0.1 = sampleTable_g2_MOI0.1_30d8_iPS %>%
  filter(SampleNum %in% c("iPS_01_g2")) %>%
  dplyr::select(BC30StarcodeD8) 
#for each barcode, split UMIs across 2 splits -- need only run this once
splitA_g2_MOI0.1 = vector()
splitB_g2_MOI0.1 = vector()
set.seed(2059)
binomsplit_g2_MOI0.1 = apply(sampleTable_g2_MOI0.1_30d8_iPS,1,survive)
for (i in 1:nrow(sampleTable_g2_MOI0.1_30d8_iPS)){
  barcode = sampleTable_g2_MOI0.1_30d8_iPS$BC30StarcodeD8[i]
  splitA_g2_MOI0.1 = c(splitA_g2_MOI0.1, rep(barcode,sum(binomsplit_g2_MOI0.1[[i]]==0)))
  splitB_g2_MOI0.1 = c(splitB_g2_MOI0.1, rep(barcode,sum(binomsplit_g2_MOI0.1[[i]]==1)))
} 
saveRDS(splitA_g2_MOI0.1,file=paste0(plotDirectory,'g2MOI0.1/simulation/splitA_g2_MOI0.1.rds'))
saveRDS(splitB_g2_MOI0.1,file=paste0(plotDirectory,'g2MOI0.1/simulation/splitB_g2_MOI0.1.rds'))
#skip to here to read in RDS files
splitA_g2_MOI0.1 = readRDS(file=paste0(plotDirectory,'g2MOI0.1/simulation/splitA_g2_MOI0.1.rds'))
splitB_g2_MOI0.1 = readRDS(file=paste0(plotDirectory,'g2MOI0.1/simulation/splitB_g2_MOI0.1.rds'))
set.seed(2059)
#factor in loss of some UMIs: for g2MOI0.1 survivability = 0.05
survivability_g2_MOI0.1 = 0.05
survivalA_g2_MOI0.1 <- sample(length(splitA_g2_MOI0.1),length(splitA_g2_MOI0.1)*survivability_g2_MOI0.1)
splitAsurvive_g2_MOI0.1 <- splitA_g2_MOI0.1[survivalA_g2_MOI0.1]

survivalB_g2_MOI0.1 <- sample(length(splitB_g2_MOI0.1),length(splitB_g2_MOI0.1)*survivability_g2_MOI0.1)
splitBsurvive_g2_MOI0.1 <- splitB_g2_MOI0.1[survivalB_g2_MOI0.1]

venn.diagram(
  x = list(unlist(alliPS_g2_MOI0.1), splitAsurvive_g2_MOI0.1),
  category.names = c("iPS" , "splitAsim"),
  filename = paste0(plotDirectory,'g2MOI0.1/simulation/splitAiPS_survivability',survivability_g2_MOI0.1,'.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)
venn.diagram(
  x = list(unlist(alliPS_g2_MOI0.1), splitBsurvive_g2_MOI0.1),
  category.names = c("iPS" , "splitBsim"),
  filename = paste0(plotDirectory,'g2MOI0.1/simulation/splitBiPS_survivability',survivability_g2_MOI0.1,'.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)

splitAbarcodes_g2_MOI0.1 = unique(splitAsurvive_g2_MOI0.1)
splitBbarcodes_g2_MOI0.1 = unique(splitBsurvive_g2_MOI0.1)

venn.diagram(
  x = list(splitAsurvive_g2_MOI0.1, splitBsurvive_g2_MOI0.1),
  category.names = c("splitAsim" , "splitBsim"),
  filename = paste0(plotDirectory,'g2MOI0.1/simulation/splitA_splitBoverlap',survivability_g2_MOI0.1,'.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-15, 15),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)

#g2_MOI0.23
sampleTable_g2_MOI0.23_30d8_diff = sampleTable_g2_MOI0.23 %>%
  filter(SampleNum %in% c("diff_023A_g2","diff_023B_g2")) %>%
  dplyr::select(countsperUMI,BC30StarcodeD8,SampleNum) %>% #74884 UMIs for A, 44591 for B using 30d8 - why the discrepancy huh
  group_by(BC30StarcodeD8, SampleNum) %>%
  summarise(UMI = n(), counts = sum(countsperUMI)) #collapse and count how many distinct UMIs per barcode
#47540 

sampleTable_g2_MOI0.23_30d8_iPS = sampleTable_g2_MOI0.23 %>%
  filter(SampleNum %in% c("iPS_023_g2")) %>%
  dplyr::select(countsperUMI,BC30StarcodeD8,SampleNum) %>% #2816890 UMIs using 30d8
  group_by(BC30StarcodeD8, SampleNum) %>%
  summarise(UMI = n(), counts = sum(countsperUMI)) #collapse and count how many distinct UMIs per barcode
#70265 

#determine observed overlap
splitAbar_g2_MOI0.23 = sampleTable_g2_MOI0.23_30d8_diff %>%
  filter(SampleNum == "diff_023A_g2") %>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 26964 barcodes
splitBbar_g2_MOI0.23 = sampleTable_g2_MOI0.23_30d8_diff %>%
  filter(SampleNum == "diff_023B_g2") %>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 20576 barcodes
split_initialbar_g2_MOI0.23 = sampleTable_g2_MOI0.23_30d8_iPS %>%
  filter(SampleNum == "iPS_023_g2")%>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 70265 barcodes

venn.diagram(
  x = list(unlist(splitAbar_g2_MOI0.23), unlist(splitBbar_g2_MOI0.23)),
  category.names = c("splitA" , "splitB"),
  filename = paste0(plotDirectory,'g2MOI0.23/nocutoff/g2MOI0.23_VDdiffs.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.015, 0.015),
  cat.fontfamily = "sans"
) 
venn.diagram(
  x = list(unlist(split_initialbar_g2_MOI0.23), unlist(splitAbar_g2_MOI0.23)),
  category.names = c("iPS" , "splitA"),
  filename = paste0(plotDirectory,'g2MOI0.23/nocutoff/g2MOI0.23_VDsplitAiPS.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.015, 0.015),
  cat.fontfamily = "sans"
)
venn.diagram(
  x = list(unlist(split_initialbar_g2_MOI0.23), unlist(splitBbar_g2_MOI0.23)),
  category.names = c("iPS" , "splitB"),
  filename = paste0(plotDirectory,'g2MOI0.23/nocutoff/g2MOI0.23_VDsplitBiPS.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.015, 0.015),
  cat.fontfamily = "sans"
)

#simulated overlap
alliPS_g2_MOI0.23 = sampleTable_g2_MOI0.23_30d8_iPS %>%
  filter(SampleNum %in% c("iPS_01_g2")) %>%
  dplyr::select(BC30StarcodeD8) 
#for each barcode, split UMIs across 2 splits -- need only run this once
splitA_g2_MOI0.23 = vector()
splitB_g2_MOI0.23 = vector()
set.seed(2059)
binomsplit_g2_MOI0.23 = apply(sampleTable_g2_MOI0.23_30d8_iPS,1,survive)
for (i in 1:nrow(sampleTable_g2_MOI0.23_30d8_iPS)){
  barcode = sampleTable_g2_MOI0.23_30d8_iPS$BC30StarcodeD8[i]
  splitA_g2_MOI0.23 = c(splitA_g2_MOI0.23, rep(barcode,sum(binomsplit_g2_MOI0.23[[i]]==0)))
  splitB_g2_MOI0.23 = c(splitB_g2_MOI0.23, rep(barcode,sum(binomsplit_g2_MOI0.23[[i]]==1)))
} 
saveRDS(splitA_g2_MOI0.23,file=paste0(plotDirectory,'g2MOI0.23/simulation/splitA_g2_MOI0.23.rds'))
saveRDS(splitB_g2_MOI0.23,file=paste0(plotDirectory,'g2MOI0.23/simulation/splitB_g2_MOI0.23.rds'))
#skip to here to read in RDS files
splitA_g2_MOI0.23 = readRDS(file=paste0(plotDirectory,'g2MOI0.23/simulation/splitA_g2_MOI0.23.rds'))
splitB_g2_MOI0.23 = readRDS(file=paste0(plotDirectory,'g2MOI0.23/simulation/splitB_g2_MOI0.23.rds'))
set.seed(2059)
#factor in loss of some UMIs: for g2MOI0.23 survivability = 0.03
survivability_g2_MOI0.23 = 0.03
survivalA_g2_MOI0.23 <- sample(length(splitA_g2_MOI0.23),length(splitA_g2_MOI0.23)*survivability_g2_MOI0.23)
splitAsurvive_g2_MOI0.23 <- splitA_g2_MOI0.23[survivalA_g2_MOI0.23]

survivalB_g2_MOI0.23 <- sample(length(splitB_g2_MOI0.23),length(splitB_g2_MOI0.23)*survivability_g2_MOI0.23)
splitBsurvive_g2_MOI0.23 <- splitB_g2_MOI0.23[survivalB_g2_MOI0.23]

venn.diagram(
  x = list(unlist(alliPS_g2_MOI0.23), splitAsurvive_g2_MOI0.23),
  category.names = c("iPS" , "splitAsim"),
  filename = paste0(plotDirectory,'g2MOI0.23/simulation/splitAiPS_survivability',survivability_g2_MOI0.23,'.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)
venn.diagram(
  x = list(unlist(alliPS_g2_MOI0.23), splitBsurvive_g2_MOI0.23),
  category.names = c("iPS" , "splitBsim"),
  filename = paste0(plotDirectory,'g2MOI0.23/simulation/splitBiPS_survivability',survivability_g2_MOI0.23,'.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)

splitAbarcodes_g2_MOI0.23 = unique(splitAsurvive_g2_MOI0.23)
splitBbarcodes_g2_MOI0.23 = unique(splitBsurvive_g2_MOI0.23)

venn.diagram(
  x = list(splitAsurvive_g2_MOI0.23, splitBsurvive_g2_MOI0.23),
  category.names = c("splitAsim" , "splitBsim"),
  filename = paste0(plotDirectory,'g2MOI0.23/simulation/splitA_splitBoverlap',survivability_g2_MOI0.23,'.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-15, 15),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)

#g2_MOI1.0
#split barcodes into iPS and differentiated objects
sampleTable_g2_MOI1.0_30d8_diff = sampleTable_g2_MOI1.0 %>%
  filter(SampleNum %in% c("diff_10A_g2","diff_10B_g2")) %>%
  dplyr::select(countsperUMI,BC30StarcodeD8,SampleNum) %>% #404941 UMIs for A, 357501 UMIs for B using 30d8
  group_by(BC30StarcodeD8, SampleNum) %>%
  summarise(UMI = n(), counts = sum(countsperUMI)) #collapse and count how many distinct UMIs per barcode
#175112

sampleTable_g2_MOI1.0_30d8_iPS = sampleTable_g2_MOI1.0 %>%
  filter(SampleNum %in% c("iPS_10_g2")) %>%
  dplyr::select(countsperUMI,BC30StarcodeD8,SampleNum) %>% #14774677 UMIs using 30d8
  group_by(BC30StarcodeD8, SampleNum) %>%
  summarise(UMI = n(), counts = sum(countsperUMI)) #collapse and count how many distinct UMIs per barcode
#213210

#determine observed overlap
splitAbar_g2_MOI1.0 = sampleTable_g2_MOI1.0_30d8_diff %>%
  filter(SampleNum == "diff_10A_g2") %>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 90324 barcodes
splitBbar_g2_MOI1.0 = sampleTable_g2_MOI1.0_30d8_diff %>%
  filter(SampleNum == "diff_10B_g2") %>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 84788 barcodes
split_initialbar_g2_MOI1.0 = sampleTable_g2_MOI1.0_30d8_iPS %>%
  filter(SampleNum == "iPS_10_g2")%>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 213210 barcodes

venn.diagram(
  x = list(unlist(splitAbar_g2_MOI1.0), unlist(splitBbar_g2_MOI1.0)),
  category.names = c("splitA" , "splitB"),
  filename = paste0(plotDirectory,'g2MOI1.0/nocutoff/g2MOI1.0_VDdiffs.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.015, 0.015),
  cat.fontfamily = "sans"
) 
venn.diagram(
  x = list(unlist(split_initialbar_g2_MOI1.0), unlist(splitAbar_g2_MOI1.0)),
  category.names = c("iPS" , "splitA"),
  filename = paste0(plotDirectory,'g2MOI1.0/nocutoff/g2MOI1.0_VDsplitAiPS.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.015, 0.015),
  cat.fontfamily = "sans"
)
venn.diagram(
  x = list(unlist(split_initialbar_g2_MOI1.0), unlist(splitBbar_g2_MOI1.0)),
  category.names = c("iPS" , "splitB"),
  filename = paste0(plotDirectory,'g2MOI1.0/nocutoff/g2MOI1.0_VDsplitBiPS.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.015, 0.015),
  cat.fontfamily = "sans"
)

#simulated overlap
alliPS_g2_MOI1.0 = sampleTable_g2_MOI1.0_30d8_iPS %>%
  filter(SampleNum %in% c("iPS_01_g2")) %>%
  dplyr::select(BC30StarcodeD8) 
#for each barcode, split UMIs across 2 splits -- need only run this once
splitA_g2_MOI1.0 = vector()
splitB_g2_MOI1.0 = vector()
set.seed(2059)
binomsplit_g2_MOI1.0 = apply(sampleTable_g2_MOI1.0_30d8_iPS,1,survive)
for (i in 1:nrow(sampleTable_g2_MOI1.0_30d8_iPS)){
  barcode = sampleTable_g2_MOI1.0_30d8_iPS$BC30StarcodeD8[i]
  splitA_g2_MOI1.0 = c(splitA_g2_MOI1.0, rep(barcode,sum(binomsplit_g2_MOI1.0[[i]]==0)))
  splitB_g2_MOI1.0 = c(splitB_g2_MOI1.0, rep(barcode,sum(binomsplit_g2_MOI1.0[[i]]==1)))
} 
saveRDS(splitA_g2_MOI1.0,file=paste0(plotDirectory,'g2MOI1.0/simulation/splitA_g2_MOI1.0.rds'))
saveRDS(splitB_g2_MOI1.0,file=paste0(plotDirectory,'g2MOI1.0/simulation/splitB_g2_MOI1.0.rds'))
#skip to here to read in RDS files
splitA_g2_MOI1.0 = readRDS(file=paste0(plotDirectory,'g2MOI1.0/simulation/splitA_g2_MOI1.0.rds'))
splitB_g2_MOI1.0 = readRDS(file=paste0(plotDirectory,'g2MOI1.0/simulation/splitB_g2_MOI1.0.rds'))
set.seed(2059)
#factor in loss of some UMIs: for g2MOI1.0 survivability = 0.03
survivability_g2_MOI1.0 = 0.03
survivalA_g2_MOI1.0 <- sample(length(splitA_g2_MOI1.0),length(splitA_g2_MOI1.0)*survivability_g2_MOI1.0)
splitAsurvive_g2_MOI1.0 <- splitA_g2_MOI1.0[survivalA_g2_MOI1.0]

survivalB_g2_MOI1.0 <- sample(length(splitB_g2_MOI1.0),length(splitB_g2_MOI1.0)*survivability_g2_MOI1.0)
splitBsurvive_g2_MOI1.0 <- splitB_g2_MOI1.0[survivalB_g2_MOI1.0]

venn.diagram(
  x = list(unlist(alliPS_g2_MOI1.0), splitAsurvive_g2_MOI1.0),
  category.names = c("iPS" , "splitAsim"),
  filename = paste0(plotDirectory,'g2MOI1.0/simulation/splitAiPS_survivability',survivability_g2_MOI1.0,'.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)
venn.diagram(
  x = list(unlist(alliPS_g2_MOI1.0), splitBsurvive_g2_MOI1.0),
  category.names = c("iPS" , "splitBsim"),
  filename = paste0(plotDirectory,'g2MOI1.0/simulation/splitBiPS_survivability',survivability_g2_MOI1.0,'.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)

splitAbarcodes_g2_MOI1.0 = unique(splitAsurvive_g2_MOI1.0)
splitBbarcodes_g2_MOI1.0 = unique(splitBsurvive_g2_MOI1.0)

venn.diagram(
  x = list(splitAsurvive_g2_MOI1.0, splitBsurvive_g2_MOI1.0),
  category.names = c("splitAsim" , "splitBsim"),
  filename = paste0(plotDirectory,'g2MOI1.0/simulation/splitA_splitBoverlap',survivability_g2_MOI1.0,'.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-15, 15),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)


##############silencing venn diagrams from g1 MOI 1.0
#overlap in GFP+ vs GFP- conditions
splitAposbar = sampleTable_g1_MOI1.0_30d8_diff %>%
  filter(SampleNum %in% c("diff_10AGFPpos_g1")) %>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 141210 barcodes - more than expected
splitAnegbar = sampleTable_g1_MOI1.0_30d8_diff %>%
  filter(SampleNum %in% c("diff_10AGFPneg_g1")) %>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 118304 barcodes - more than expected
splitBposbar = sampleTable_g1_MOI1.0_30d8_diff %>%
  filter(SampleNum %in% c("diff_10BGFPpos_g1")) %>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 126046 barcodes - more than expected
splitBnegbar = sampleTable_g1_MOI1.0_30d8_diff %>%
  filter(SampleNum %in% c("diff_10BGFPneg_g1")) %>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 147737 barcodes - more than expected

venn.diagram(
  x = list(unlist(splitAposbar), unlist(splitAnegbar),unlist(splitBposbar), unlist(splitBnegbar)),
  category.names = c("splitApos" , "splitAneg","splitBpos" , "splitBneg"),
  filename = paste0(plotDirectory,'g1MOI1.0/nocutoff/g1MOI1.0_VDdiffsorts.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  #cat.pos = c(-27, 27),
  #cat.dist = c(0.105, 0.105),
  cat.fontfamily = "sans"
)
#pairwise
venn.diagram(
  x = list(unlist(splitAposbar), unlist(splitAnegbar)),
  category.names = c("splitApos" , "splitAneg"),
  filename = paste0(plotDirectory,'g1MOI1.0/nocutoff/g1MOI1.0_VDsplitAsilencing.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.105, 0.105),
  cat.fontfamily = "sans"
)

venn.diagram(
  x = list(unlist(splitBposbar), unlist(splitBnegbar)),
  category.names = c("splitBpos" , "splitBneg"),
  filename = paste0(plotDirectory,'g1MOI1.0/nocutoff/g1MOI1.0_VDsplitBsilencing.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.105, 0.105),
  cat.fontfamily = "sans"
)

venn.diagram(
  x = list(unlist(splitAposbar), unlist(splitBposbar)),
  category.names = c("splitApos" , "splitBpos"),
  filename = paste0(plotDirectory,'g1MOI1.0/nocutoff/g1MOI1.0_VDGFPposdiffs.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.105, 0.105),
  cat.fontfamily = "sans"
)

venn.diagram(
  x = list(unlist(splitAnegbar), unlist(splitBnegbar)),
  category.names = c("splitAneg" , "splitBneg"),
  filename = paste0(plotDirectory,'g1MOI1.0/nocutoff/g1MOI1.0_VDGFPnegdiffs.pdf'),
  output=TRUE,
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.105, 0.105),
  cat.fontfamily = "sans"
)



