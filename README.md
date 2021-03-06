# cardiacdiff_fatemap

This repo contains scripts needed to produce all graphs and images contained in figures for the paper:

Jiang CL, et al. Cell type determination for cardiac differentiation occurs soon after seeding of human induced pluripotent stem cells. (2022)

For any questions, please contact CLJ (cjiang93 at gmail) and the corresponding authors for this paper, AR (arjunrajlab at gmail) and RJ (jainr at pennmedicine.upenn.edu)

In order to reproduce all graphs and images in this paper:

Download and install dependencies for this code:
R v4.0.5
R packages VennDiagram_1.6.20, reshape2_1.4.4, DescTools_0.99.41, gridExtra_2.3, Seurat_4.0.2, dplyr_1.0.6, tidyverse_1.3.1, and their associated dependencies

Download all extractedData (35.7GB of intermediate-processed experiment data) files in the structure provided for this repository from Dropbox: https://www.dropbox.com/sh/pcihymkterrvvz3/AACTEQWG2KKQw3bi0NqShsWYa?dl=0.
Your top-level project directory should have subdirectories finalFigures, rawData, plotScripts, extractedData, contractileMovies, additionalgDNAbarcodeScripts, 10XbarcodepipelineScripts, extractionScripts, and plotData, as well as 3 supplementary table Excel files. 

Edit directory variables (i.e. "SeuratDirectory", "barcodeDirectory", "dataDirectory", "extractDirectory") variables to reflect your local file paths for this project repository and the associated downloaded files. Set "plotDirectory" to reflect your desired location for finished plots.

Additional annotation/information available in individual R scripts.

Scripts to produce extractedData are also included in extractionScripts - these use some but not all of the files included in the rawData subdirectory of the Dropbox folder referenced above.

Scripts to produce figures from extractedData can be found in plotScripts.

If running scripts using downloaded GEO processed single cell RNA sequencing files (of the structure <samplename>.matrix.mtx.gz, <samplename>.features.tsv.gz, <samplename>.barcodes.tsv.gz), please first make separate directories for each <samplename> containing the matrix, features, and barcodes files and then remove <samplename> from those; i.e., <samplename>/matrix.mtx.gz and other files).
