#!/bin/bash

<<Usage

service=$(basename $(dirname $(pwd)))
cat samples_id.txt | while read -r line
do echo "qsub -V -b y -j y -cwd -l h_vmem=10G -pe openmp 8 -q all.q -N ASSM_${line} bash 01_Trimmer_assembler_validater_BU_ISCIII.sh $line $service"
done > _01_commands_assembler_$service.sh

bash _01_commands_assembler_$service.sh

Usage


sample=$1
service=$2
serviceDir="/processing_Data/bioinformatics/services_and_colaborations/CNM/bacteriologia"
sourcedir=$serviceDir/$service/ANALYSIS/00-reads
trimmeddir=$serviceDir/$service/ANALYSIS/02-preprocessing
contigdir=$serviceDir/$service/ANALYSIS/05-assembly
fastQCdir=$serviceDir/$service/ANALYSIS/01-fastQC
preprocQCdir=$serviceDir/$service/ANALYSIS/03-preprocQC


#echo find /processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/RAW/Prueba_ensamblado -name  "*_1*" -type f | awk '/'"${sample}"'/ && /fastq.gz/'


mkdir -p $contigdir/$sample
mkdir -p $trimmeddir/$sample


##Check quality with Fastqc of raw samples

echo "################################## Checking $sample ###################################################################"


mkdir -p $fastQCdir/$sample/
fastqc $(find -L $sourcedir -name "*fastq*" -type f | awk '/'"${sample}"'/') -o $fastQCdir/$sample/

#fastqc $sourcedir/$sample/$sample_*


echo "################################## DONE Checking $sample ###################################################################"


##Trim with Trimmomatic in order to get better quality sequences

echo "################################ Trimming $sample ###########################################################"
R1=$(find -L $sourcedir -name "*fastq*" -type f | awk '/'"${sample}"'/ && /1.fastq/')
R2=$(find -L $sourcedir -name "*fastq*" -type f | awk '/'"${sample}"'/ && /2.fastq/')


echo "Sourcedir = " $sourcedir
echo "R1 = " $R1
echo "R2 = " $R2

java -jar -Xmx10G /opt/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 8 \
$R1 \
$R2 \
$trimmeddir/$sample/$sample"_R1_filtered.fastq.gz" \
$trimmeddir/$sample/$sample"_R1_unpaired.fastq.gz" \
$trimmeddir/$sample/$sample"_R2_filtered.fastq.gz" \
$trimmeddir/$sample/$sample"_R2_unpaired.fastq.gz" \
ILLUMINACLIP:/opt/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:40

echo "################################ DONE Trimming $sample ###########################################################"



##Check quality with Fastqc of trimmed samples

echo "################################## Checking trimmed $sample ###################################################################"

mkdir -p $preprocQCdir/$sample/
fastqc $(find $trimmeddir/$sample -name "*fastq.gz" -type f) -o $preprocQCdir/$sample/


#fastqc $trimmeddir/$sample/$sample_*

echo "################################## DONE Checking trimmed $sample ###################################################################"


##Asemble with spades

R1P=$(find $trimmeddir/$sample -name  "*1_filtered.fastq.gz" -type f)
R2P=$(find $trimmeddir/$sample -name  "*2_filtered.fastq.gz" -type f)
R1U=$(find $trimmeddir/$sample -name  "*1_unpaired.fastq.gz" -type f)
R2U=$(find $trimmeddir/$sample -name  "*2_unpaired.fastq.gz" -type f)

echo "################################ Assembling $sample with SPAdes ###########################################################"

spades.py --careful -t 8 -k 61,81,101,121 \
--pe1-1 $R1P \
--pe1-2 $R2P \
--s1 $R1U \
--s2 $R2U \
-o $contigdir/$sample

echo "################################ DONE Assembling $sample ###########################################################"


#echo "################################ Assembling $sample with Velvet ###########################################################"

#velveth $contigdir 21,121,10 -fastq.gz -shortPaired -separate $R1P $R2P

#for f in  $contigdir_*
#do 
#  velvetg $f -cov_cutoff auto -unused_reads yes
#done



echo "################################ ALL DONE with $sample ###########################################################"



