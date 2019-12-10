#!/bin/bash


PathToDropSeqTools=/home/goodwintracking/Genomics/Software/Drop-seq_tools-2.1.0
PICARD=/home/goodwintracking/Genomics/Software/Drop-seq-master/lib/picard-2.18.14.jar
WorkingDirectory=/mnt/LocalData/Genomics/VNC/20190709_analysis
DataLocation1=/mnt/LocalData/Genomics/VNC/20190709_analysis/RawData/190109_K00181_0103_AHYKH7BBXX/FASTQ
DataLocation2=/mnt/LocalData/Genomics/VNC/20190709_analysis/RawData/190206_K00198_0396_BH3CN2BBXY/FASTQ
DataLocation3=/mnt/LocalData/Genomics/VNC/20190709_analysis/RawData/190214_K00198_0398_BH3C3KBBXY/FASTQ
PathToGenome=/home/goodwintracking/Genomics/Genome/dm613_TG_clean

CurrDir=$(pwd)
ExpName=${CurrDir##*/}


# Stage 3: sort aligned reads (STAR does not necessarily emit reads in the same order as the input)
echo Sorting the aligned reads
if [ -f 08_${ExpName}_star_Aligned.out.sam ]; then
	java -Dsamjdk.buffer_size=131072 \
	 -XX:GCTimeLimit=50 \
	 -XX:GCHeapFreeLimit=10 \
	 -Xmx4000m -jar \
	 ${PICARD} SortSam \
		INPUT=08_${ExpName}_star_Aligned.out.sam \
		OUTPUT=09_${ExpName}_aligned_sorted.bam \
		SORT_ORDER=queryname
fi

# Stage 4: merge and tag aligned reads
echo Merging the tags with the aligned reads
if [ -f 09_${ExpName}_aligned_sorted.bam ]; then
	java -Xmx4000m -jar \
	 ${PICARD} MergeBamAlignment \
		REFERENCE_SEQUENCE=${PathToGenome}/dm613_TG_clean.fa \
		UNMAPPED_BAM=06_${ExpName}_unaligned_mc_tagged_polyA_filtered.bam \
		ALIGNED_BAM=09_${ExpName}_aligned_sorted.bam \
		INCLUDE_SECONDARY_ALIGNMENTS=false \
		PAIRED_RUN=false \
		OUTPUT=10_${ExpName}_merged.bam
fi

echo Tag read with gene expression
if [ -f 10_${ExpName}_merged.bam ]; then
	${PathToDropSeqTools}/./TagReadWithGeneFunction \
		O=11_${ExpName}_merged_gene_exon_tagged.bam \
		ANNOTATIONS_FILE=${PathToGenome}/dm613_TG_clean.refFlat \
		CREATE_INDEX=true \
		INPUT=10_${ExpName}_merged.bam
fi

echo Make digital expression matrix
if [ -f 11_${ExpName}_merged_gene_exon_tagged.bam ]; then
	${PathToDropSeqTools}/./DigitalExpression \
		I=11_${ExpName}_merged_gene_exon_tagged.bam \
		O=12_${ExpName}_merged_gene_exon_tagged_top30000cells.dge.txt.gz \
		SUMMARY=${ExpName}_summaries/12_out_gene_exon_tagged_top30000cells_dge_summary.txt \
		NUM_CORE_BARCODES=30000
fi



echo some sort of histogram ...
if [ -f 12_${ExpName}_merged_gene_exon_tagged_top30000cells.dge.txt.gz ]; then
	${PathToDropSeqTools}/./BamTagHistogram \
		I=11_${ExpName}_merged_gene_exon_tagged.bam \
		O=13_${ExpName}_merged_gene_exon_tagged_ReadsPerCellBarcode.dge.txt \
		TAG=XC
fi


exit
