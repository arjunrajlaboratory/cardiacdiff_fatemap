README — Figures_12.18.21_SuppFig5.txt

A. Venn diagrams for observed and simulated overlap between differentiated splits from additional infection conditions (MOI 0.1 x 2 replicates, MOI 0.23, MOI 1.0 x 2 replicates) generated using the script plotScripts/survivalVD.R. By condition, venn diagrams can be found at: 
	MOI 0.1, replicate 1: 
		simulated: plotData/Supplement/gDNA_VD/g1MOI0.1/simulation/splitA_splitBoverlap0.015.pdf; 
		observed: plotData/Supplement/gDNA_VD/g1MOI0.1/nocutoff/g1MOI0.1_VDdiffs.pdf
	MOI 0.1, replicate 2:
		simulated: plotData/Supplement/gDNA_VD/g2MOI0.1/simulation/splitA_splitBoverlap0.05.pdf; 
		observed: plotData/Supplement/gDNA_VD/g2MOI0.1/nocutoff/g2MOI0.1_VDdiffs.pdf
	MOI 0.23:
		simulated: plotData/Supplement/gDNA_VD/g2MOI0.23/simulation/splitA_splitBoverlap0.03.pdf; 
		observed: plotData/Supplement/gDNA_VD/g2MOI0.23/nocutoff/g2MOI0.23_VDdiffs.pdf
	MOI 1.0, replicate 1:
		simulated: plotData/Supplement/gDNA_VD/g1MOI1.0/simulation/splitA_splitBoverlap0.18.pdf; 
		observed: plotData/Supplement/gDNA_VD/g1MOI1.0/nocutoff/g1MOI1.0_VDdiffsorts.pdf
	MOI 1.0, replicate 2:
		simulated: plotData/Supplement/gDNA_VD/g2MOI1.0/simulation/splitA_splitBoverlap0.03.pdf; 
		observed: plotData/Supplement/gDNA_VD/g2MOI1.0/nocutoff/g2MOI1.0_VDdiffs.pdf

B. Histograms (grey) demonstrating the range associated with 10000 simulated fractions of barcodes that overlap between splits A and B as compared to observed (teal line) for MOI 0.5 condition (corresponding venn diagram in Figure 3B) were generated with the script plotScripts/survivaloverlapsignificance.R and can be found at plotData/Supplement/gDNA_VD/fractionoverlapsplitA.pdf and plotData/Supplement/gDNA_VD/fractionoverlapsplitB.pdf.

C. All venn diagrams generated using the annotated section of the script plotScripts/survivalVD.R. Venn diagrams demonstrating observed overlap in barcodes found in GFP positive and GFP negative sorted MOI 1.0 differentiated cells for each split can be found at plotData/Supplement/gDNA_VD/g1MOI0.1/nocutoff/g1MOI1.0_VDsplitAsilencing for split A and plotData/Supplement/gDNA_VD/g1MOI0.1/nocutoff/g1MOI1.0_VDsplitBsilencing for split B. Venn diagram demonstrating observed overlap in barcodes found in GFP positive sorted differentiated cells across splits can be found at plotData/Supplement/gDNA_VD/g1MOI0.1/nocutoff/g1MOI1.0_VDGFPposdiffs.pdf.

D. Table of the top 20 barcodes (by total number of cells labeled) generated in Adobe Illustrator CC 2020 based on "print(top100,n=20)" command in plotScripts/Figure3.R.