#try to determine whether survival overlap is significant
#July2021

extractDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotData/Figure3/"
plotDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/plotData/Supplement/gDNA_VD/"
barcodeDirectory <- "~/Dropbox (RajLab)/cardiac_diff_scRNAseq/paper/extractedData/gDNAbarcodeobjects/"

#pull simulated splitA and splitB from g1 MOI0.5 (Figure 3)
#NEED TO RUN FIGURE 3 SCRIPT BEFORE RUNNING THIS SCRIPT
splitA = readRDS(file=paste0(extractDirectory,'VD_gDNA/g1MOI0.5/simulation/splitA.rds'))
splitB = readRDS(file=paste0(extractDirectory,'VD_gDNA/g1MOI0.5/simulation/splitB.rds'))

#pull observed data
sampleTable_g1_MOI0.5 <- readRDS(paste0(barcodeDirectory,'sampleTable_g1_MOI0.5.rds'))
####split barcodes into iPS and differentiated objects
sampleTable_g1_MOI0.5_30d8_diff = sampleTable_g1_MOI0.5 %>%
  filter(SampleNum %in% c("diff_05A_g1","diff_05B_g1")) %>%
  dplyr::select(countsperUMI,BC30StarcodeD8,SampleNum) %>% #293974 UMIs for A, 269517 for B, using 30d8
  group_by(BC30StarcodeD8, SampleNum) %>%
  summarise(UMI = n(), counts = sum(countsperUMI)) #collapse and count how many distinct UMIs per barcode
#142854 barcodes bc not collapsing between samples
####determine actual observed overlap
splitAbar = sampleTable_g1_MOI0.5_30d8_diff %>%
  filter(SampleNum == "diff_05A_g1") %>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 72782 barcodes -- a bit more than expected w/o selection
splitBbar = sampleTable_g1_MOI0.5_30d8_diff %>%
  filter(SampleNum == "diff_05B_g1") %>%
  dplyr::select(BC30StarcodeD8) %>%
  unique() #sometimes same barcode will have multiple UMIs: 70072 barcodes

#resample overlap 1000 times with survivability = 0.11
n = 1000
seeds = 1:n
fractionoverlapA = c()
fractionoverlapB = c()
survive = 0.11
for (j in seeds){
  set.seed(j)
  survivalA <- sample(length(splitA),length(splitA)*survive)
  splitAsurvive <- splitA[survivalA]
  survivalB <- sample(length(splitB),length(splitB)*survive)
  splitBsurvive <- splitB[survivalB]
  splitAbarcodes = unique(splitAsurvive)
  splitBbarcodes = unique(splitBsurvive)
  fractionoverlapA[j] = length(intersect(splitAsurvive,splitBsurvive))/length(splitAbarcodes)
  fractionoverlapB[j] = length(intersect(splitAsurvive,splitBsurvive))/length(splitBbarcodes)
}
obsoverlapA = nrow(intersect(splitAbar,splitBbar))/nrow(splitAbar)
obsoverlapB = nrow(intersect(splitAbar,splitBbar))/nrow(splitBbar)

meanoverlapA = mean(fractionoverlapA) #average of average random distance matrix
sdoverlapA = sd(fractionoverlapA)
pA = pnorm(obsoverlapA,mean=meanoverlapA,sd=sdoverlapA,lower.tail=FALSE)
if(round(pA,3)==0){
  pA3 <- "< 0.001"
}else{
  pA3 <- paste("=",round(pA,3))
}
overlapAplot = ggplot() +
  geom_histogram(aes(x=fractionoverlapA),binwidth = 0.001,fill='gray80') +
  geom_vline(xintercept=obsoverlapA,color='#009292') +
  ggtitle(paste0('split A overlap','\np ',pA3)) +
  xlab("fraction of split A barcodes found in split B") +
  coord_cartesian(xlim = c(0.75, 0.8),ylim = c(0,300))+ #without coord_cartesian throws an error about rows being removed
  theme_bw()
pdf(paste(plotDirectory,"fractionoverlapsplitA.pdf", sep=""), width=5, height=5)
print(overlapAplot)
dev.off()

meanoverlapB = mean(fractionoverlapB) #average of average random distance matrix
sdoverlapB = sd(fractionoverlapB)
pB = pnorm(obsoverlapB,mean=meanoverlapB,sd=sdoverlapB,lower.tail=FALSE)
if(round(pB,3)==0){
  pB3 <- "< 0.001"
}else{
  pB3 <- paste("=",round(pB,3))
}
overlapBplot = ggplot() +
  geom_histogram(aes(x=fractionoverlapB),binwidth = 0.001,fill='gray80') +
  geom_vline(xintercept=obsoverlapB,color='#009292') +
  ggtitle(paste0('split B overlap','\np ',pB3)) +
  xlab("fraction of split B barcodes found in split A") +
  coord_cartesian(xlim = c(0.75, 0.8), ylim = c(0,300))+ #without coord_cartesian throws an error about rows being removed
  theme_bw()
pdf(paste(plotDirectory,"fractionoverlapsplitB.pdf", sep=""), width=5, height=5)
print(overlapBplot)
dev.off()
