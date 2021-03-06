# VNC_scRNAseq
Analysis of scRNA-seq of the Drosophila melanogaster adult ventral nerve cord from females and males.

These data and analyses were published in "A single-cell transcriptomic atlas of the adult Drosophila ventral nerve cord" (Allen, Neville, et al., 2020), avaible [here](https://elifesciences.org/articles/54074).

Fastq sequencing data were processed using Drop-seq_tools v2.1.0, and was implemented with the '10x_parallelprocessing.sh', '10x_preprocessing.sh', and '10x_postprocessing.sh'. Resulting digital expression matrices were analysed with Seurat v2.3.4 using the 'R_preprocessing.R' script. The 'all_analysis.R' script contains information regarding the subsequent analysis for each figure of the manuscript. 'myDotPlot.R' is a modified version of Seurat v2.3.4's 'DotPlot' function to accept 3 colours when plotting the scaled expression. 'R_preprocessing_brain_data.R' script includes the analysis of the mid-brain drop-seq data (Croset et al 2018) and the brain 10x data (Davie et al 2018) to compare to the VNC data. 'sessionInfo.txt' contains the package version used in this analysis.  
