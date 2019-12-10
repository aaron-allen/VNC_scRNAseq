#!/bin/bash

printf "\n"
echo $(date)
echo "GOOD MORNING/AFTERNOON/EVENING (DELETE WHERE APPROPRIATE)"
printf "\n"
echo TIME TO SET UP THE FILE PATHS!
printf "\n"



PathToDropSeqTools=/home/goodwintracking/Genomics/Software/Drop-seq_tools-2.1.0
PICARD=/home/goodwintracking/Genomics/Software/Drop-seq-master/lib/picard-2.18.14.jar
WorkingDirectory=/mnt/LocalData/Genomics/VNC/20190709_analysis
DataLocation1=/mnt/LocalData/Genomics/VNC/20190709_analysis/RawData/190109_K00181_0103_AHYKH7BBXX/FASTQ
DataLocation2=/mnt/LocalData/Genomics/VNC/20190709_analysis/RawData/190206_K00198_0396_BH3CN2BBXY/FASTQ
DataLocation3=/mnt/LocalData/Genomics/VNC/20190709_analysis/RawData/190214_K00198_0398_BH3C3KBBXY/FASTQ
PathToGenome=/home/goodwintracking/Genomics/Genome/dm613_TG_clean

MasterDirectory=$(pwd)



########
######## To generate reference genome, start with .fa file, add transgenes, then run:
STAR --runMode genomeGenerate --genomeDir ${PathToGenome} --runThreadN 12 --genomeFastaFiles ${PathToGenome}/dm613_TG_clean.fa
java -jar ${PICARD} CreateSequenceDictionary R=${PathToGenome}/dm613_TG_clean.fa O=${PathToGenome}/dm613_TG_clean.dict





# # ########
# # # CHANGE NAME OF FASTQs:

echo MAKING DIRECTORIES AND CHANGING NAMES OF FASTQ FILES

cd ${WorkingDirectory}
mkdir -p Preprocessing
cd ${WorkingDirectory}/Preprocessing


echo MaleVNCRep1
mkdir MaleVNCRep1
cat ${DataLocation1}/628685_22/628685_22_S6_L003_R1_001.fastq.gz \
	${DataLocation2}/640576_22/640576_22_S66_L008_R1_001.fastq.gz \
	${DataLocation3}/643569_22/643569_22_S26_L002_R1_001.fastq.gz \
	${DataLocation3}/643570_22/643570_22_S10_L003_R1_001.fastq.gz \
	${DataLocation3}/643571_22/643571_22_S42_L004_R1_001.fastq.gz \
	${DataLocation3}/643572_22/643572_22_S50_L005_R1_001.fastq.gz \
	${DataLocation3}/643573_22/643573_22_S2_L006_R1_001.fastq.gz \
	${DataLocation3}/643574_22/643574_22_S18_L007_R1_001.fastq.gz \
	> ${WorkingDirectory}/Preprocessing/MaleVNCRep1/MaleVNCRep1_1.fastq.gz

cat ${DataLocation1}/628685_22/628685_22_S6_L003_R2_001.fastq.gz \
	${DataLocation2}/640576_22/640576_22_S66_L008_R2_001.fastq.gz \
	${DataLocation3}/643569_22/643569_22_S26_L002_R2_001.fastq.gz \
	${DataLocation3}/643570_22/643570_22_S10_L003_R2_001.fastq.gz \
	${DataLocation3}/643571_22/643571_22_S42_L004_R2_001.fastq.gz \
	${DataLocation3}/643572_22/643572_22_S50_L005_R2_001.fastq.gz \
	${DataLocation3}/643573_22/643573_22_S2_L006_R2_001.fastq.gz \
	${DataLocation3}/643574_22/643574_22_S18_L007_R2_001.fastq.gz \
	> ${WorkingDirectory}/Preprocessing/MaleVNCRep1/MaleVNCRep1_2.fastq.gz



echo MaleVNCRep2
mkdir MaleVNCRep2
cat ${DataLocation1}/628685_24/628685_24_S8_L003_R1_001.fastq.gz \
	${DataLocation2}/640576_24/640576_24_S68_L008_R1_001.fastq.gz \
	${DataLocation3}/643569_24/643569_24_S28_L002_R1_001.fastq.gz \
	${DataLocation3}/643570_24/643570_24_S12_L003_R1_001.fastq.gz \
	${DataLocation3}/643571_24/643571_24_S44_L004_R1_001.fastq.gz \
	${DataLocation3}/643572_24/643572_24_S52_L005_R1_001.fastq.gz \
	${DataLocation3}/643573_24/643573_24_S4_L006_R1_001.fastq.gz \
	${DataLocation3}/643574_24/643574_24_S20_L007_R1_001.fastq.gz \
	> ${WorkingDirectory}/Preprocessing/MaleVNCRep2/MaleVNCRep2_1.fastq.gz

cat ${DataLocation1}/628685_24/628685_24_S8_L003_R2_001.fastq.gz \
	${DataLocation2}/640576_24/640576_24_S68_L008_R2_001.fastq.gz \
	${DataLocation3}/643569_24/643569_24_S28_L002_R2_001.fastq.gz \
	${DataLocation3}/643570_24/643570_24_S12_L003_R2_001.fastq.gz \
	${DataLocation3}/643571_24/643571_24_S44_L004_R2_001.fastq.gz \
	${DataLocation3}/643572_24/643572_24_S52_L005_R2_001.fastq.gz \
	${DataLocation3}/643573_24/643573_24_S4_L006_R2_001.fastq.gz \
	${DataLocation3}/643574_24/643574_24_S20_L007_R2_001.fastq.gz \
	> ${WorkingDirectory}/Preprocessing/MaleVNCRep2/MaleVNCRep2_2.fastq.gz



echo FemaleVNCRep1
mkdir FemaleVNCRep1
cat ${DataLocation1}/628685_26/628685_26_S10_L003_R1_001.fastq.gz \
	${DataLocation2}/640576_26/640576_26_S70_L008_R1_001.fastq.gz \
	${DataLocation3}/643569_26/643569_26_S30_L002_R1_001.fastq.gz \
	${DataLocation3}/643570_26/643570_26_S14_L003_R1_001.fastq.gz \
	${DataLocation3}/643571_26/643571_26_S46_L004_R1_001.fastq.gz \
	${DataLocation3}/643572_26/643572_26_S54_L005_R1_001.fastq.gz \
	${DataLocation3}/643573_26/643573_26_S6_L006_R1_001.fastq.gz \
	${DataLocation3}/643574_26/643574_26_S22_L007_R1_001.fastq.gz \
	> ${WorkingDirectory}/Preprocessing/FemaleVNCRep1/FemaleVNCRep1_1.fastq.gz

cat ${DataLocation1}/628685_26/628685_26_S10_L003_R2_001.fastq.gz \
	${DataLocation2}/640576_26/640576_26_S70_L008_R2_001.fastq.gz \
	${DataLocation3}/643569_26/643569_26_S30_L002_R2_001.fastq.gz \
	${DataLocation3}/643570_26/643570_26_S14_L003_R2_001.fastq.gz \
	${DataLocation3}/643571_26/643571_26_S46_L004_R2_001.fastq.gz \
	${DataLocation3}/643572_26/643572_26_S54_L005_R2_001.fastq.gz \
	${DataLocation3}/643573_26/643573_26_S6_L006_R2_001.fastq.gz \
	${DataLocation3}/643574_26/643574_26_S22_L007_R2_001.fastq.gz \
	> ${WorkingDirectory}/Preprocessing/FemaleVNCRep1/FemaleVNCRep1_2.fastq.gz



echo FemaleVNCRep2
mkdir FemaleVNCRep2
cat ${DataLocation1}/628685_28/628685_28_S12_L003_R1_001.fastq.gz \
	${DataLocation2}/640576_28/640576_28_S72_L008_R1_001.fastq.gz \
	${DataLocation3}/643569_28/643569_28_S32_L002_R1_001.fastq.gz \
	${DataLocation3}/643570_28/643570_28_S16_L003_R1_001.fastq.gz \
	${DataLocation3}/643571_28/643571_28_S48_L004_R1_001.fastq.gz \
	${DataLocation3}/643572_28/643572_28_S56_L005_R1_001.fastq.gz \
	${DataLocation3}/643573_28/643573_28_S8_L006_R1_001.fastq.gz \
	${DataLocation3}/643574_28/643574_28_S24_L007_R1_001.fastq.gz \
	> ${WorkingDirectory}/Preprocessing/FemaleVNCRep2/FemaleVNCRep2_1.fastq.gz

cat ${DataLocation1}/628685_28/628685_28_S12_L003_R2_001.fastq.gz \
	${DataLocation2}/640576_28/640576_28_S72_L008_R2_001.fastq.gz \
	${DataLocation3}/643569_28/643569_28_S32_L002_R2_001.fastq.gz \
	${DataLocation3}/643570_28/643570_28_S16_L003_R2_001.fastq.gz \
	${DataLocation3}/643571_28/643571_28_S48_L004_R2_001.fastq.gz \
	${DataLocation3}/643572_28/643572_28_S56_L005_R2_001.fastq.gz \
	${DataLocation3}/643573_28/643573_28_S8_L006_R2_001.fastq.gz \
	${DataLocation3}/643574_28/643574_28_S24_L007_R2_001.fastq.gz \
	> ${WorkingDirectory}/Preprocessing/FemaleVNCRep2/FemaleVNCRep2_2.fastq.gz









cd ${WorkingDirectory}/Preprocessing
printf "\n\n"
echo STARTING PREPROCESSING
for A in */
do
	echo CHANGING DIRECTORY INTO $A
	cd $A
	cp -r $MasterDirectory/10x_preprocessing.sh 10x_preprocessing.sh
	bash ./10x_preprocessing.sh &
	sleep 5s    # 5 second lag to allow JAVA to open
	# Check if JAVA is running
	while [ $(pgrep -c "java") -gt 3 ]
	do
		sleep 5s
	done
	# echo Changing directory to ${WorkingDirectory}
	cd ${WorkingDirectory}/Preprocessing
done

# Wait for all preprocessing to finish before going on to next section
echo JUST WAITING FOR PREPROCESSING TO FINISH
while pgrep -x "java" > /dev/null
do
	sleep 5s
done



cd ${WorkingDirectory}/Preprocessing
printf "\n\n"
echo STARTING ALIGNMENT
for A in */
do
	echo CHANGING DIRECTORY INTO $A
	cd $A
	echo Aligning with Star
	if [ -f 07_${A%/}_unaligned_mc_tagged_polyA_filtered.fastq ]; then
		STAR --genomeDir ${PathToGenome}/ \
		 --runThreadN 12 \
		 --outFileNamePrefix 08_${A%/}_star_ \
		 --alignIntronMax 200000 \
		 --sjdbGTFfile ${PathToGenome}/dmel-all-r6.13.gtf \
		 --readFilesIn ./07_${A%/}_unaligned_mc_tagged_polyA_filtered.fastq
	fi
	mv 08_${A%/}_star_Log* ./${A%/}_summaries/
	# echo CHANGING DIRECTORY TO ${WorkingDirectory}
	cd ${WorkingDirectory}/Preprocessing
done




cd ${WorkingDirectory}
mkdir -p R/
mkdir -p R/InputFiles

cd ${WorkingDirectory}/Preprocessing


printf "\n\n"
echo STARTING POSTPROCESSING
for A in */
do
	echo Changing Directory into $A
	cd $A
	cp -r $MasterDirectory/10x_postprocessing.sh 10x_postprocessing.sh
	bash ./10x_postprocessing.sh &
	sleep 5s    # 5 second lag to allow JAVA to open
	# Check if JAVA is running
	while [ $(pgrep -c "java") -gt 3 ]
	do
		sleep 5s
	done
	# echo Changing directory to ${WorkingDirectory}
	cd ${WorkingDirectory}/Preprocessing
done
# Wait for all preprocessing to finish before going on to next section
echo JUST WAITING FOR POSTPROCESSING TO FINISH
while pgrep -x "java" > /dev/null
do
	sleep 5s
done



printf "\n\n"
echo All Done.
echo $(date)

