#!/bin/bash

<<Usage
#Copy this command to generate a file (commands_placnet_by_sample.txt) with the order to execute this script for each sample


group=references_simulation
samplesid=/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/DOC/SAMPLES/samplesid_$group.txt


cat $samplesid | while read -r line
do echo "qsub -V -b y -j y -l h_vmem=20G -pe openmp 10 -cwd -q all.q -N MAP_${line} bash 1_filter_reads_to_plasmid_ddbb.sh $line $group"
done > 1_commands_plasmid_mapping_$group.txt

bash 1_commands_plasmid_mapping_$group.txt
Usage

###############1 Map reads on plasmid database & filter positive matches

#-o/--offrate If reporting many alignments per read, try reducing bowtie2-build offrate
#--threads Increasing the number of threads will speed up the index building considerably in most cases: Use qsub flag -pe openmp 10 and --threads 10 (same number of threads)
#-N <int> Sets the number of mismatches to allowed in a seed alignment during multiseed alignment. 

#samtools index uses BAI index by default

###############2 Find plasmid with more coverage
###############2 de novo assembly with matched reads
###############3 Blast contigs with all plasmid DDBB
###############4 Find the most represented plasmid using Megablast output
###############5 Cluster hit plasmids (clust after in order to get more species specific plasmids)
###############6 Map reads against those
###############7 Check which plasmids have better coverage


###############1 Map reads on plasmid database & filter positive matches
#IN: illumina reads
#OUT: sorted BAM and fastq of each pair matching plasmid

#plasmidDdbb="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/plasmid.all.genomic.fna"
#plasmidDdbb="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/DOC/DOC_AB/Db_mcr1_plasmids.bin"
plasmidDdbb="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/plasmid.all.genomic.feb212017.fasta"
sample=$1
group=$2
trimmedDir="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/ANALYSIS/TRIMMED"
mappedDir="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/ANALYSIS/MAPPING/PLASMIDS"


mkdir -p $mappedDir/$group/$sample

#<<How_to_index_plasmid_ddbb

bowtie2-build --offrate 3 $plasmidDdbb $plasmidDdbb

#How_to_index_plasmid_ddbb


echo "##### Mapping $sample of $group Group "

bowtie2 -a -x $plasmidDdbb \
-1 $trimmedDir/$group/$sample/$sample"_1_paired.fastq.gz" \
-2 $trimmedDir/$group/$sample/$sample"_2_paired.fastq.gz" \
-S $mappedDir/$group/$sample/$sample.sam

echo "##### DONE Mapping $sample of $group Group"



echo "###### Converting SAM to sorted indexed BAM in $sample of $group Group"

samtools view -Sb $mappedDir/$group/$sample/$sample.sam -o $mappedDir/$group/$sample/$sample.bam

samtools sort -T $mappedDir/$group/$sample/$sample".sorted.bam" -o $mappedDir/$group/$sample/$sample".sorted.bam" $mappedDir/$group/$sample/$sample.bam

samtools index $mappedDir/$group/$sample/$sample".sorted.bam"

rm $mappedDir/$group/$sample/$sample.sam
rm $mappedDir/$group/$sample/$sample.bam

echo "###### DONE Converting SAM to sorted indexed BAM in $sample of $group Group"

#echo "##### Filtering matching reads in $sample of $group Group"
#-f 67 = read paired + read mapped in proper pair + first in pair (66 also work)
#samtools view -f 67 $mappedDir/$group/$sample/$sample".sorted.bam" | awk '{if($3 != "*") print "@" $1 "\n" $10 "\n" "+" $1 "\n" $11}' > $mappedDir/$group/$sample/$sample"_R1.fastq"
#-f 131 = read paired + read mapped in proper pair + second in pair (230 also work)
#samtools view -f 131 $mappedDir/$group/$sample/$sample".sorted.bam" | awk '{if($3 != "*") print "@" $1 "\n" $10 "\n" "+" $1 "\n" $11}' > $mappedDir/$group/$sample/$sample"_R2.fastq"

#echo "##### DONE Filtering matching reads in $sample of $group Group"