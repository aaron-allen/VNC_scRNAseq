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
mkdir ${ExpName}_summaries

echo Runing Fastq To Sam for ${ExpName}
if [[ -f ${ExpName}_1.fastq.gz && -f ${ExpName}_2.fastq.gz ]]; then
	java -jar ${PICARD} FastqToSam \
		F1=${ExpName}_1.fastq.gz \
		F2=${ExpName}_2.fastq.gz \
		O=01_${ExpName}_unmapped.bam \
		SM=${ExpName} 
fi



# Stage 1: pre-alignment tag and trim

# cellular tag
echo Extracting Cellular Tag
if [ -f 01_${ExpName}_unmapped.bam ]; then
	${PathToDropSeqTools}/./TagBamWithReadSequenceExtended \
		SUMMARY=${ExpName}_summaries/01to02_unaligned_tagged_Cellular.bam_summary.txt \
		BASE_RANGE=1-16 \
		BASE_QUALITY=10 \
		BARCODED_READ=1 \
		DISCARD_READ=false \
		TAG_NAME=XC \
		NUM_BASES_BELOW_QUALITY=1 \
		INPUT=01_${ExpName}_unmapped.bam \
		OUTPUT=02_${ExpName}_unaligned_tagged_Cell.bam
fi

# molecular tag
echo Extracting Molecular Tag
if [ -f 02_${ExpName}_unaligned_tagged_Cell.bam ]; then
	${PathToDropSeqTools}/./TagBamWithReadSequenceExtended \
		SUMMARY=${ExpName}_summaries/02to03_unaligned_tagged_Molecular.bam_summary.txt \
		BASE_RANGE=17-26 \
		BASE_QUALITY=10 \
		BARCODED_READ=1 \
		DISCARD_READ=true \
		TAG_NAME=XM \
		NUM_BASES_BELOW_QUALITY=1 \
		INPUT=02_${ExpName}_unaligned_tagged_Cell.bam \
		OUTPUT=03_${ExpName}_unaligned_tagged_CellMolecular.bam
fi

# quality filter
echo Quality Filter
if [ -f 03_${ExpName}_unaligned_tagged_CellMolecular.bam ]; then
	${PathToDropSeqTools}/./FilterBam \
		TAG_REJECT=XQ \
		INPUT=03_${ExpName}_unaligned_tagged_CellMolecular.bam \
		OUTPUT=04_${ExpName}_unaligned_tagged_filtered.bam
fi

# read trimming
echo Trimming polyA tail
if [ -f 04_${ExpName}_unaligned_tagged_filtered.bam ]; then
	${PathToDropSeqTools}/./PolyATrimmer \
		OUTPUT=06_${ExpName}_unaligned_mc_tagged_polyA_filtered.bam \
		OUTPUT_SUMMARY=${ExpName}_summaries/04to06_polyA_trimming_report.txt \
		MISMATCHES=0 \
		NUM_BASES=12 \
		INPUT=04_${ExpName}_unaligned_tagged_filtered.bam
fi

# Stage 2: alignment
echo Converting Sam to Fastq
if [ -f 06_${ExpName}_unaligned_mc_tagged_polyA_filtered.bam ]; then
	java -Xmx500m -jar \
	 ${PICARD} SamToFastq \
		INPUT=06_${ExpName}_unaligned_mc_tagged_polyA_filtered.bam \
		FASTQ=07_${ExpName}_unaligned_mc_tagged_polyA_filtered.fastq
fi


exit