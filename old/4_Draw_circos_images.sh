#!/bin/bash

<<Usage
#Copy this command to generate a file (commands_placnet_by_sample.txt) with the order to execute this script for each sample


group=references_simulation
samplesid=/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/DOC/SAMPLES/samplesid_$group.txt

cat $samplesid | while read -r line
do echo "qsub -V -b y -j y -l h_vmem=20G -pe openmp 10 -cwd -q all.q -N CRC_${line}_PLSMD bash 4_Draw_circos_images.sh $line $group"
done > 4_commands_Draw_circos_images_$group.txt

bash 4_commands_Draw_circos_images_$group.txt
Usage

###############1 Map reads on plasmid database & filter positive matches
###############2 Find plasmid with more coverage
###############3 Extract plasmid with more coverage (>90%) and tracks



#IN: 


sample=$1
group=$2
plasmidDdbb="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/plasmid.all.genomic.feb212017.fasta"
abReference="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/ARGannot.fasta"
trimmedDir="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/ANALYSIS/TRIMMED" 
mappedDir="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/ANALYSIS/MAPPING/PLASMIDS"
imageDir="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/ANALYSIS/IMAGES/PLASMIDS"
contigDir="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/ANALYSIS/ASSEMBLY"
#plasmidLenght="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/NCBI_PLASMID_LENGTH.tsv"
plasmidLenght="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/plasmid.all.genomic.feb212017.length"


circosConfFileAll="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/DOC/CONFIG_FILES/circos_all_plasmids.conf"
circosConfFileInd="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/DOC/CONFIG_FILES/circos_individual_plasmid.conf"
circosDir=$imageDir/$group/$sample"/CIRCOS"
imageName=$sample".png"

plasmidMapped=$mappedDir/$group/$sample/$sample".ac.covered.txt"
plasmidMappedFasta=$mappedDir/$group/$sample/$sample".ac.covered.fasta"
sampleCoverage=$mappedDir/$group/$sample/$sample".coverage.txt"

karyotypeFile_80=$imageDir/$group/$sample/$sample".karyotype.txt_80"
karyotypeFile=$imageDir/$group/$sample/$sample".karyotype.txt"
abrFile=$imageDir/$group/$sample/$sample".abr.bed.coords"
replisomeFile=$imageDir/$group/$sample/$sample".pfinder.bed.coords"
isFile=$imageDir/$group/$sample/$sample".isfinder.bed.coords"
coverageFile=$imageDir/$group/$sample/$sample".plasmid.bedgraph_80"
cdsContigFile=$imageDir/$group/$sample/$sample".bed.coords"
cdsDdbbFile="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/plasmid.all.genomic.feb212017.bed"
contigFile=$imageDir/$group/$sample/$sample".plasmids.blast.contigs"
linksFile=$imageDir/$group/$sample/$sample".plasmids.blast.links.sorted"


#PLASMID_KARYOTYPE
#karyotype = /processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/ANALYSIS/IMAGES/PLASMIDS/blaOXA48/A14_03/A14_03_karyotype.txt
#PLASMID_ANTIBIOTIC_RESISTANCE
#file = /processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/TMP/A14_03_CIRCOS/A14_03.abr.bed.coords
#PLASMID_COVERAGE_GRAPH
#file = /processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/TMP/A14_03_plasmid.bedgraph_80
#PLASMID_CDS_CONTIG
#file = /processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/TMP/A14_03_CIRCOS/A14_03.bed.coords
#PLASMID_CDS_DDBB
#file = /processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/plasmid.all.genomic.ac.id.bed
#PLASMID_CONTIGS
#file = processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/TMP/A14_03_CIRCOS/A14_03.blast.contigs
#PLASMID_LINKS
#file = /processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/TMP/A14_03_CIRCOS/A14_03.blast.links.sorted

mkdir -p $circosDir

echo "Creating config file for circos in SAMPLE $sample FILE $circosDir/$sample.circos.conf"

awk '{gsub("PLASMID_KARYOTYPE","'$karyotypeFile'"); \
gsub("PLASMID_ANTIBIOTIC_RESISTANCE","'$abrFile'"); \
gsub("PLASMID_REPLISOME_PLASMIDFINDER","'$replisomeFile'"); \
gsub("PLASMID_IS_ISFINDER","'$isFile'"); \
gsub("PLASMID_COVERAGE_GRAPH","'$coverageFile'"); \
gsub("PLASMID_CDS_CONTIG","'$cdsContigFile'"); \
gsub("PLASMID_CDS_DDBB","'$cdsDdbbFile'"); \
gsub("PLASMID_CONTIGS","'$contigFile'"); \
gsub("PLASMID_LINKS","'$linksFile'"); \
gsub("OUTPUTDIR","'$circosDir'"); \
gsub("IMAGENAME","'$imageName'"); \
print $0}' $circosConfFileAll > $circosDir/$sample".circos.conf"

echo "DONE Creating config file for circos in SAMPLE $sample"


echo "Creating config file for circos in SAMPLE $sample FILE $circosDir/$sample.circos.conf"

awk '{gsub("PLASMID_KARYOTYPE","'$karyotypeFile_80'"); \
gsub("PLASMID_ANTIBIOTIC_RESISTANCE","'$abrFile'"); \
gsub("PLASMID_REPLISOME_PLASMIDFINDER","'$replisomeFile'"); \
gsub("PLASMID_IS_ISFINDER","'$isFile'"); \
gsub("PLASMID_COVERAGE_GRAPH","'$coverageFile'"); \
gsub("PLASMID_CDS_CONTIG","'$cdsContigFile'"); \
gsub("PLASMID_CDS_DDBB","'$cdsDdbbFile'"); \
gsub("PLASMID_CONTIGS","'$contigFile'"); \
gsub("PLASMID_LINKS","'$linksFile'"); \
gsub("OUTPUTDIR","'$circosDir'"); \
print $0}' $circosConfFileInd > $circosDir/$sample"_individual.circos.conf"

echo "DONE Creating config file for circos in SAMPLE $sample"

echo "Executing circos in SAMPLE $sample FILE $circosDir/$sample.circos.conf"

circos -conf $circosDir/$sample".circos.conf"

#SAMPLE_SHOWN
#IMAGENAME_SAMPLE

for i in $(cat $plasmidMapped"_80")
do
	awk '{gsub("SAMPLE_SHOWN","'$i'"); \
	gsub("IMAGENAME_SAMPLE_PLASMID","'$sample'_'$i'.png"); \
	print $0}' $circosDir/$sample"_individual.circos.conf" > $circosDir/$sample"_"$i"_individual.circos.conf"
	circos -conf $circosDir/$sample"_"$i"_individual.circos.conf"
done 



echo "ALL DONE"
