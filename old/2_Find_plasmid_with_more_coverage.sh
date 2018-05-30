#!/bin/bash

<<Usage
#Copy this command to generate a file (commands_placnet_by_sample.txt) with the order to execute this script for each sample


group=eferences_simulation
samplesid=/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/DOC/SAMPLES/samplesid_$group.txt

cat $samplesid | while read -r line
do echo "qsub -V -b y -j y -l h_vmem=20G -pe openmp 10 -cwd -q all.q -hold_jid MAP_${line} -N COV_${line}_PLSMD bash 2_Find_plasmid_with_more_coverage.sh $line $group"
done > 2_commands_plasmid_coverage_$group.txt

bash 2_commands_plasmid_coverage_$group.txt
Usage

###############1 Map reads on plasmid database & filter positive matches

###############2 Find plasmid with more coverage
#bedtools genomecov [OPTIONS] -i <bed/gff/vcf> -g <genome>
#-ibam           The input file is in BAM format.
#bedtools output
#1.chromosome (or entire genome)
#2 depth of coverage from features in input file
#3.number of bases on chromosome (or genome) with depth equal to column 2.
#4.size of chromosome (or entire genome) in base pairs
#5.fraction of bases on chromosome (or entire genome) with depth equal to column 2.
#$ bedtools genomecov -i A.bed -g my.genome
#ref_name depth numberOfReads totalLength percentageWithDepth
#chr1   0  980  1000  0.98
#chr1   1  20   1000  0.02
#chr2   1  500  500   1
#genome 0  980  1500  0.653333
#genome 1  520  1500  0.346667

#CD-HIT
#-i input filename in fasta format, required
#-o output filename, required
#-c sequence identity threshold, default 0.9 
#-M max available memory (Mbyte), default 400
#-d length of description in .clstr file, default 20. If set to 0, it takes the fasta defline and stops at first space
#-s length difference cutoff, default 0.0 if set to 0.9, the shorter sequences need to be at least 90% length of the representative of the cluster  
#-B 1 or 0, default 0, by default, sequences are stored in RAM if set to 1, sequence are stored on hard drive it is recommended to use -B 1 for huge databases

###############2 de novo assembly with matched reads
###############3 Blast contigs with all plasmid DDBB
###############4 Find the most represented plasmid using Megablast output
###############5 Cluster hit plasmids (clust after in order to get more species specific plasmids)
###############6 Map reads against those
###############7 Check which plasmids have better coverage



###############2 Find plasmid with more coverage
#IN: sorted BAM
#OUT: list of ac numbers covered more than 90%
#OUT: bedgraph of those sequences

sample=$1
group=$2
#plasmidReference="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/plasmid.all.genomic.fna"
#plasmidDdbb="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/DOC/DOC_AB/Db_mcr1_plasmids.bin"
plasmidDdbb="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/plasmid.all.genomic.feb212017.fasta"
trimmedDir="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/ANALYSIS/TRIMMED"
mappedDir="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/ANALYSIS/MAPPING/PLASMIDS"
#plasmidLenght="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/NCBI_PLASMID_LENGTH.tsv"
plasmidLenght="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/plasmid.all.genomic.feb212017.length"
plasmidMapped=$mappedDir/$group/$sample/$sample".ac.covered.txt"
plasmidMappedFasta=$mappedDir/$group/$sample/$sample".ac.covered.fasta"
sampleCoverage=$mappedDir/$group/$sample/$sample".coverage.txt"



echo "###### Calculating coverage for every plamid"

bedtools genomecov -ibam $mappedDir/$group/$sample/$sample".sorted.bam" -g $plasmidLenght > $mappedDir/$group/$sample/$sample".coverage.txt"

#NZ_CP011540.1   0       165286  174066  0.949559
#NZ_CP011540.1   1       3450    174066  0.0198201
#NZ_CP011540.1   2       3445    174066  0.0197913
#NZ_CP011540.1   3       558     174066  0.00320568
#NZ_CP011540.1   4       449     174066  0.00257948
#NZ_CP011540.1   5       88      174066  0.000505555
#NZ_CP011540.1   6       27      174066  0.000155114
#NZ_CP011540.1   8       2       174066  1.14899e-05
#NZ_CP011540.1   9       2       174066  1.14899e-05
#NZ_CP011540.1   10      7       174066  4.02146e-05

echo "##### DONE Calculating coverage for every plamid"

echo "##### Extracting list of mapped plamid"

#Filter coverage file, taking all AC numbers whose unmaped reads (column 2 = 0) are all (column 5 = 1)
# sorting is optional 
##IN THIS STEP SEQUENCES STARTING WITH NG_ (Sequence fragments) ARE REMOVED

#OLD DDBB awk '{split($1, ac,"|")};  {if ($2 == 0 && $5 < 0.5 && $1 != "genome" && $1 !~ /NG_/) {print ac[4]}}' $mappedDir/$group/$sample/$sample".coverage.txt" > $plasmidMapped

awk '{if ($2 == 0 && $5 < 0.5 && $1 != "genome" && $1 !~ /NG_/) {print $1}}' $mappedDir/$group/$sample/$sample".coverage.txt" > $plasmidMapped

#NZ_CP016764.1                                                                                                                        
#NC_012555.1                                                                                                                          
#NC_012556.1
#NC_025004.1
#NZ_CP008825.1
#NZ_CP008899.1
#NZ_CP008906.1
#NZ_CP012170.1
#NZ_CP013028.1
#NZ_CP015021.1

echo "##### DONE Extracting list of mapped plamid"

echo "##### Obtaining FASTA from mapped plasmids"

for i in $(cat $plasmidMapped)
do
	awk 'BEGIN {RS=">"} /'"$i"'/ {print ">"$0}' $plasmidDdbb
done > $plasmidMappedFasta

#>NZ_CP016764.1 Citrobacter freundii strain B38 plasmid pOZ181, complete sequence
#TATCTTAAATAAAAAAGGGCGAGTTCTCTCGCCCTTCGTCACAAGAACACGTTTAGTTCAGCCGGTTACC
#TTTTTCATCAACTACCTTCTCACCGTCCTCTTTCGTGAACGCGCTTTTCTGCGCATCGGGAAGAATATCC
#AGCACGACTTCAGAAGGGCGGCATAGCCGGGTACCCAGCGGCGTTACCACGATAGGGCGATTGATCAGAA
#CAGGGTGCTGCAACATAAAATCGATTAACTGCTCGTCAGTAAAGCGGTCTTCCGCCAGCCCTAACGCTTC
#AAAAGGTTCGACATTCTTACGCAGCAGCGCTCGTACCGTGATCCCCATATCCGCAATGAGTTTTACCAGC
#TCAGCTCGTGATGGTGGGGTTTCAAGATAATGAATAACGGTCGGCTCTGTACCGCTGTTGCGGATCATTT
#CAAGGGTATTGCGTGACGTGCCGCAGGCCGGGTTGTGATAAATGGTAATGTTGCTCATATCAGTATCTCA
#TTACAAAGTGAAAGACAGACGAAGCGCCAGTGCTGCAAGCGTGACAAACAGCACGGGGATTGTCATGACA
#ATGCCCACCCGGAAGTAATATCCCCAGGTAATTTTGATATTTTTCTGCGACAGAACGTGCAGCCACAGTA

echo "##### DONE Obtaining FASTA from mapped plasmids"



echo "##### Clustering plasmids with mapped reads"

cd $mappedDir/$group/$sample/

cd-hit -i $plasmidMappedFasta -o $plasmidMappedFasta"_80" -c 0.70 -n 4 -d 0 -s 0.8 -B 1 -M 100000 -T 0

#cd-hit-est -i $plasmidMappedFasta -o $plasmidMappedFasta"_80" -c 0.80 -n 4 -d 0 -B 1 -M 100000 -T 0

#cd-hit-est -i $plasmidMappedFasta -o $plasmidMappedFasta"_85" -c 0.85 -n 10 -d 0 -s 0.9 -B 1 -M 100000 -T 0

echo "##### DONE Clustering plasmids with mapped reads"


echo "##### Obtaining coverage of clustered plasmids"

for i in $(cat $plasmidMappedFasta"_80"| grep ">" | awk '{gsub(">","");print $1}');do awk '/'"$i"'/' $sampleCoverage; done > $sampleCoverage"_80"

#NZ_CP016764.1                                                                                                                        
#NC_012555.1                                                                                                                          
#NC_012556.1
#NC_025004.1
#NZ_CP008825.1
#NZ_CP008899.1
#NZ_CP008906.1
#NZ_CP012170.1
#NZ_CP013028.1
#NZ_CP015021.1

#Other way to extract AC using Regex: grep ">" | awk '{gsub(">","");print substr($0,index($0,NC_), (index($0,".1 ")+1))}'

#OLD DDBB for i in $(cat $plasmidMappedFasta"_85"| grep ">" | cut -d "|" -f 4);do awk '/'"$i"'/' $sampleCoverage; done > $sampleCoverage"_85"

echo "##### DONE Obtaining coverage of clustered plasmids"


echo "##### Retrieving ac number of sequences mapped more than 80% and its FASTA"


awk '{if ($2 == 0 && $5 < 0.2 && $1 != "genome") {print $1}}' $sampleCoverage"_80" > $plasmidMapped"_80"

awk '{if ($2 == 0 && $5 < 0.2 && $1 != "genome") {print $1, ((1 - $5)*100)}}' $sampleCoverage"_80" > $plasmidMapped"_percentage_80"

#NZ_CP011540.1   0       165286  174066  0.949559
#NZ_CP011540.1   1       3450    174066  0.0198201
#NZ_CP011540.1   2       3445    174066  0.0197913
#NZ_CP011540.1   3       558     174066  0.00320568
#NZ_CP011540.1   4       449     174066  0.00257948
#NZ_CP011540.1   5       88      174066  0.000505555
#NZ_CP011540.1   6       27      174066  0.000155114
#NZ_CP011540.1   8       2       174066  1.14899e-05
#NZ_CP011540.1   9       2       174066  1.14899e-05
#NZ_CP011540.1   10      7       174066  4.02146e-05


#awk '{split($1, ac,"|")};  {if ($2 == 0 && $5 < 0.2 && $1 != "genome") {print ac[4]}}' $sampleCoverage"_80" > $plasmidMapped"_80"

echo "##### DONE Retrieving ac number of sequences mapped more than 80% and its FASTA"

#in ordet to know the complete fasta header
#for i in $(awk '{split($1, ac,"|")};  {if ($2 == 0 && $5 < 0.1 && $1 != "genome") {print ac[4]}}' A14_03_coverage.txt_80); do awk '/'"$i"'/' A14_03_ac_covered.fasta_80; done

echo "##### Obtaining FASTA from plasmid mapped more than 80%"
 
for i in $(cat $plasmidMapped"_80")
do
	awk 'BEGIN {RS=">"} /'"$i"'/ {print ">"$0}' $plasmidMappedFasta"_80"
done > $plasmidMappedFasta"_80_80"

#>NZ_CP016764.1 Citrobacter freundii strain B38 plasmid pOZ181, complete sequence
#TATCTTAAATAAAAAAGGGCGAGTTCTCTCGCCCTTCGTCACAAGAACACGTTTAGTTCAGCCGGTTACC
#TTTTTCATCAACTACCTTCTCACCGTCCTCTTTCGTGAACGCGCTTTTCTGCGCATCGGGAAGAATATCC
#AGCACGACTTCAGAAGGGCGGCATAGCCGGGTACCCAGCGGCGTTACCACGATAGGGCGATTGATCAGAA
#CAGGGTGCTGCAACATAAAATCGATTAACTGCTCGTCAGTAAAGCGGTCTTCCGCCAGCCCTAACGCTTC
#AAAAGGTTCGACATTCTTACGCAGCAGCGCTCGTACCGTGATCCCCATATCCGCAATGAGTTTTACCAGC
#TCAGCTCGTGATGGTGGGGTTTCAAGATAATGAATAACGGTCGGCTCTGTACCGCTGTTGCGGATCATTT
#CAAGGGTATTGCGTGACGTGCCGCAGGCCGGGTTGTGATAAATGGTAATGTTGCTCATATCAGTATCTCA
#TTACAAAGTGAAAGACAGACGAAGCGCCAGTGCTGCAAGCGTGACAAACAGCACGGGGATTGTCATGACA
#ATGCCCACCCGGAAGTAATATCCCCAGGTAATTTTGATATTTTTCTGCGACAGAACGTGCAGCCACAGTA

echo "##### DONE Obtaining FASTA from plasmid mapped more than 80%"

#echo "##### Obtaining FASTA from plasmid mapped more than 80%"
 
#for i in $(cat $plasmidMapped"_80")
#do
#	awk 'BEGIN {RS=">"} /'"$i"'/ {print ">"$0}' $plasmidMappedFasta"_80"
#done > $plasmidMappedFasta"_80_80"

#echo "##### DONE Obtaining FASTA from plasmid mapped more than 80%"

