#!/bin/bash


<<Usage

#serviceDir="/processing_Data/bioinformatics/services_and_colaborations/CNM/bacteriologia"
#if executed from SERVICEFOLDER/ANALYSIS
service=$(basename $(dirname $(pwd)))
qsub -V -b y -j y -cwd -l h_vmem=10G -q all.q -N QST_$service bash 05_Quast_report_BU_ISCIII.sh $service

Usage


service=$1
serviceDir="/processing_Data/bioinformatics/services_and_colaborations/CNM/bacteriologia"
referenceDir="/processing_Data/bioinformatics/services_and_colaborations/CNM/bacteriologia/$service/REFERENCES"


contigdir=$serviceDir/$service/ANALYSIS/05-assembly
nametofind="scaffolds.fasta"


reference=$(find $referenceDir -name "*genomic.fna.gz" | awk '!/cds/&&!/rna/')
gff=$(find $referenceDir -name "*genomic.gff")

echo $reference
echo $gff

echo "Evaluating contig assembly in $service"

quast.py $(find $contigdir -maxdepth 2 -name $nametofind | sort) \
-R $reference \
-G $gff \
-L -o $contigdir"/QUAST_REPORT"

echo "DONE Eevaluating contig assembly in $service"