#!/bin/bash

<<Usage
#Copy this command to generate a file (commands_$group_by_sample.txt) with the order to execute this script for each sample


group=ACIN_first
samplesid=/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/DOC/SAMPLES/samplesid_$group.txt

cat $samplesid | while read -r line
do echo "qsub -V -b y -j y -cwd -l h_vmem=40G -pe openmp 5 -q all.q -hold_jid COV_${line}_PLSMD -N DRW_${line}_PLSMD bash 3_Draw_plasmid_with_more_coverage.sh $line $group"
done > 3_commands_plasmid_drawing_$group.txt

bash 3_commands_plasmid_drawing_$group.txt
Usage

#-l h_vmem=20G -pe openmp 5 

###############1 Map reads on plasmid database & filter positive matches
###############2 Find plasmid with more coverage
###############3 Draw plasmid with more coverage (>90%)

<<outfmt6
#1	 	Query label.(qseqid)
#2	 	Target or subject(database sequence or cluster centroid) label. (sseqid)
#3	 	Percent identity. (pident)
#4	 	Alignment length. (length)
#5	 	Number of mismatches. (mismatch)
#6	 	Number of gap opens. (gapopen)
#7	 	Start position in query. Query coordinates start with 1 at the first base in the sequence as it appears in the input file. For translated searches (nucleotide queries, protein targets), query start<end for +ve frame and start>end for -ve frame. (qstart)
#8	 	End position in query. (qend)
#9	 	Start position in target. Target coordinates start with 1 at the first base in sequence as it appears in the database. For untranslated nucleotide searches, target start<end for plus strand, start>end for a reverse-complement alignment. (sstart)
#10	 	End position in target. (send)
#11	 	E-value calculated using Karlin-Altschul statistics. (evalue)
#12	 	Bit score calculated using Karlin-Altschul statistics. (bitscore)
#13		length of query (qlen)
#14		Length of target (slen)

outfmt6

#IN: List of ac mapped more than 90% in plasmid DDBB clustered with CCD-HIT (90%)

#set -e

#set -u

#set -o pipefail 

sample=$1
group=$2
clusterDir="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T"

#module load bedtools2/bedtools2-2.25.0
#module load prokka/prokka-1.12
#module load ncbi-blast/ncbi_blast-2.2.30+

abReference="/processing_Data/bioinformatics/references/resistance/ARGANNOT/20170213/ARGannot.r1.fasta"
plasmidFinderReference="/processing_Data/bioinformatics/references/plasmids/plasmidfinder/plasmidFinder_01_26_2018.fsa"
isFinderReference="$clusterDir/REFERENCES/PLASMIDS/ISFINDER/ISFinder_2.fasta"
#PRIMERS_COLE1
#isFinderReference="clusterDir/REFERENCES/PLASMIDS/ISFINDER/primers_ColE1s.fasta"
<<PREVDIR
mappedDir="clusterDir/ANALYSIS/MAPPING/PLASMIDS"
imageDir="clusterDir/ANALYSIS/IMAGES/PLASMIDS"
contigDir="clusterDir/ANALYSIS/ASSEMBLY"
PREVDIR

clusterDir="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T"
plasmidIdDir=$clusterDir"/ANALYSIS/PLASMIDID"
mappedDir=$plasmidIdDir/$group/$sample/mapping
imageDir=$plasmidIdDir/$group/$sample/data
referenceDir=$plasmidIdDir/references
contigDir=$clusterDir"/ANALYSIS/ASSEMBLY"


plasmidMapped=$mappedDir/$sample".ac.covered.txt"
plasmidMappedFasta=$mappedDir/$sample".ac.covered.fasta"
sampleCoverage=$mappedDir/$sample".coverage.txt"


#For 80% use 0.2
coverageCutoff=0.3
coverageCutoffSummary=0.2
coverageCutoffPercentage=$(echo "(1 - $coverageCutoff)*100" | bc -l)
coverageCutoffSummaryPercentage=$(echo "(1 - $coverageCutoffSummary)*100" | bc -l)
#For 80% use 0.8
thresholdCdHit=0.8
thresholdCdHitPercentage=$(echo "$thresholdCdHit*100" | bc -l)

blastIdentityCutoff=90
blastLengthCutoff=500
#For 20% use 0.2
minimumLengthComplete=0.4
minimumLengthCompletePercentage=$(echo "$minimumLengthComplete*100" | bc -l)


mkdir -p $referenceDir
mkdir -p $imageDir


echo "##### Obtain list of cromosomes (idiogram) for CIRCOS karyotype file"
echo "Generating summary karyotype file with plasmids that mapped more than" $coverageCutoffSummaryPercentage"%"

awk '{if ($2 == 0 && $5 < '"${coverageCutoffSummary}"' ) {print "chr -", $1, $1, "0", $4, "id="$1}}' $sampleCoverage"_80" > $imageDir/$sample".karyotype.txt"


echo "Generating individual karyotype file with plasmids that mapped more than" $coverageCutoffPercentage"%"

awk '{if ($2 == 0 && $5 < '"${coverageCutoff}"' ) {print "chr -", $1, $1, "0", $4, "id="$1}}' $sampleCoverage"_80" > $imageDir/$sample".karyotype.txt_80"

#chr - NZ_CP015021.1 NZ_CP015021.1 0 81401  id=NZ_CP015021.1                                                           
#chr - NZ_CP015022.1 NZ_CP015022.1 0 95170  id=NZ_CP015022.1
#chr - NZ_CP015856.1 NZ_CP015856.1 0 115432  id=NZ_CP015856.1
#chr - NZ_CP017252.1 NZ_CP017252.1 0 92691  id=NZ_CP017252.1
#chr - NZ_CP017670.1 NZ_CP017670.1 0 92755  id=NZ_CP017670.1

echo "##### DONE Obtain list of cromosomes (idiogram) for CIRCOS karyotype file"

if [ -f $imageDir/$sample".plasmid.bedgraph" ];then \
	echo "Found a bedgraph file for sample" $sample;
	echo "Omitting bedgraph step"

elif [ ! -f $mappedDir/$sample".sorted.bam" ]; then \

	echo "ERROR: BAM file was not found for sample" $sample;
	echo "Exit"
	exit

else

	echo "Obtaining coverage coordinates from sequences"

	bedtools genomecov -ibam $mappedDir/$sample".sorted.bam" -bga -max 500 > $imageDir/$sample".plasmid.bedgraph"

	echo "DONE obtaining coverage coordinates from sequences"
fi


#NZ_CP011540.1   0       44983   0
#NZ_CP011540.1   44983   44985   1
#NZ_CP011540.1   44985   44989   3
#NZ_CP011540.1   44989   44990   5
#NZ_CP011540.1   44990   44991   6



echo "##### Filtering coordinate file with plasmid sequences that mapped more than" $coverageCutoffPercentage"%"
#obtain a list with matched plasmids to filter bedgraph and make drawing faster

listAcRegexp=$(printf "%s|" $(cat $plasmidMapped"_80") | sed 's/|$//g')

awk '{if ($1 ~ /'"${listAcRegexp}"'/) {print $1, $2, $3, $4}}' $imageDir/$sample".plasmid.bedgraph" > $imageDir/$sample".plasmid.bedgraph_80"


echo "##### DONE Filtering coordinate file with matched sequences"


#contigFile=$(find $contigDir/ -name "contigs.fasta" -type f | awk '/'"${sample}"'/ ')
contigFile=$(find -L $contigDir/ -name "scaffolds.fasta" -type f 2> /dev/null| awk '/'"${sample}"'/' | awk 'NR==1')
#2> /dev/null avoid permission denied error pront for non contig directories
#contigFile=$nadalesDir/A_257_index.final.scaffolds.fasta
if [ -f $imageDir/$sample".fna" -a -f $imageDir/$sample".gff"  ];then \

	echo "Found an annotation file for sample" $sample;
	echo "Omitting annotation with prokka"

elif [ ! -f $contigFile ]; then \

	echo "ERROR: File with contigs was not found for sample" $sample;
	echo "Exit"
	exit

else

	echo -e "contig file found at:" $contigFile
	echo -e "##### Anotating contigs with Prokka \n"

	prokka --force --outdir $imageDir \
	--prefix $sample \
	--addgenes \
	--kingdom Bacteria \
	--usegenus \
	--centre CNM \
	--locustag $sample \
	--compliant \
	--cpus 5 \
	$contigFile

	echo "##### DONE Anotating contigs with Prokka"

fi

#--genus Klebsiella \
#--species pneumoniae \

#From now on fasta sequences from contigs are those reformatted in prokka (program doesn't acept too large sequences names outputed by spades)


echo "##### Generatin contigs for  CIRCOS image"

if [ ! -s $plasmidMappedFasta"_80_80" ];then \
	echo "ERROR:" $(basename $plasmidMappedFasta"_80_80") "not found, exiting";
	exit
else

	echo "Building blast ddbb for" $plasmidMappedFasta"_80_80"
	makeblastdb -in $plasmidMappedFasta"_80_80" -out $imageDir/$sample".blast.tmp" -dbtype nucl
	echo "BLASTing" $sample "contigs against" $(basename $plasmidMappedFasta"_80_80")
	blastn -query $imageDir/$sample".fna" \
	-db $imageDir/$sample".blast.tmp" \
	-out $imageDir/$sample".plasmids.blast" -evalue 0.0001 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

fi






#Resulting file from blasi in outputfmt 6:
#gnl|CNM|A14_03_2	NZ_CP012884.1	100.00	122	0	0	1	122	133535	133414	2e-56	  226	484815	194742
#gnl|CNM|A14_03_2	NZ_CP012884.1	80.60	232	41	4	283550	283779	159259	159488	2e-41	  176	484815	194742
#gnl|CNM|A14_03_2	NZ_CP012884.1	87.60	121	15	0	1	121	46556	46436	6e-31	  141	484815	194742
#gnl|CNM|A14_03_2	NZ_CP018365.1	100.00	122	0	0	1	122	83332	83453	2e-56	  226	484815	260772
#gnl|CNM|A14_03_2	NZ_CP018365.1	100.00	121	0	0	1	121	10288	10408	6e-56	  224	484815	260772

awk 'BEGIN {OFS="\t"} {split($1, contigname, "_"); {if ($3 > '"${blastIdentityCutoff}"' && $4 > '"${blastLengthCutoff}"') print $2, $9, $10, contigname[length(contigname)], "id=contig_"contigname[length(contigname)]}}' \
$imageDir/$sample".plasmids.blast" > $imageDir/$sample".plasmids.blast.contigs"

#NG_048025.1     1178    1476    contig_1    id=contig_1
#NG_048025.1     1361    1127    contig_1    id=contig_1
#NZ_CP010574.1   66599   66479   contig_2    id=contig_2
#NZ_CP008930.1   122307  122187  contig_2    id=contig_2
#NZ_CP011577.1   112628  112748  contig_2    id=contig_2

echo "##### DONE Generatin contigs for CIRCOS image"

echo "##### Generatin full length contigs for CIRCOS image with cmtigs that aligned more than" $minimumLengthCompletePercentage"%"

cat $imageDir/$sample".plasmids.blast" | awk 'BEGIN{OFS="\t"}{split($1,contigname, "_")}(($3 > '"${blastIdentityCutoff}"') && (($4/$13)>'"${minimumLengthComplete}"') && (!x[$1$2]++)){{isInverted=($10-$9);ext2=($13-$8)}; \
{if (isInverted < 0) {pos1=$10; pos2=$9;} else {pos1 =$9;pos2=$10}; \
{if ((isInverted < 0) && (($14 - pos2) > $7)) {coordChr2=(pos2 + $7);} else if ((isInverted < 0) && (($14 - pos2) <= $7)) {coordChr2=$14}; \
{if ((isInverted < 0) && (ext2 <= pos1)) {coordChr1= pos1 - ext2;} else if ((isInverted < 0) && (ext2 > pos1)) {coordChr1= 1}; \
{if ((isInverted > 0) && (pos1 > $7)) {coordChr1=(pos1 - $7);} else if ((isInverted > 0) && (pos1 <= $7)) {coordChr1=1}; \
{if ((isInverted > 0) && (ext2 > ($14-pos2))) {coordChr2= $14;} else if ((isInverted > 0) && (ext2 <= ($14-pos2))) {coordChr2= (pos2 + ext2)}; \
{print $2, coordChr1, coordChr2, contigname[length(contigname)], "id="$13} }}}}}}' \
>$imageDir/$sample".plasmids.complete.coords"


cat $imageDir/$sample".plasmids.blast" | awk 'BEGIN{OFS="\t"}{split($1,contigname, "_")}(($3 > '"${blastIdentityCutoff}"') && (($4/$13)>'"${minimumLengthComplete}"') && (!x[$1$2]++)){{isInverted=($10-$9);ext2=($13-$8)}; \
{if (isInverted < 0) {pos1=$10; pos2=$9;} else {pos1 =$9;pos2=$10}; \
{if ((isInverted < 0) && (($14 - pos2) < $7)) {coordChr1=1; coordChr2=($7-($14-pos2)); {print $2, coordChr1, coordChr2, contigname[length(contigname)], "id="$13}}; \
{if ((isInverted < 0) && (ext2 > pos1)) {coordChr1=($14-(ext2-pos1)); coordChr2=$14; {print $2, coordChr1, coordChr2, contigname[length(contigname)], "id="$13}}; \
{if ((isInverted > 0) && (pos1 < $7)) {coordChr1=($14-($7-pos1)); coordChr2=$14; {print $2, coordChr1, coordChr2, contigname[length(contigname)], "id="$13}}; \
{if ((isInverted > 0) && (ext2 > ($14-pos2))) {coordChr1=1; coordChr2=(ext2-($14-pos2)); {print $2, coordChr1, coordChr2, contigname[length(contigname)], "id="$13}}}}}}}}' \
>>$imageDir/$sample".plasmids.complete.coords"

#NC_019986.1	1	3510	1	id=contig_1
#NZ_CP018444.1	1	5234	1	id=contig_1
#NZ_CP014300.1	1	4744	1	id=contig_1
#JQ319775.1		1	5225	1	id=contig_1
#NZ_CP018357.1	1	3851	1	id=contig_1

echo "##### DONE Generatin full length contigs for CIRCOS image"

echo "##### Generating link track for CIRCOS image"

#Parse blast output and obtain links between both datasets

awk '(($4/$13)>'"${minimumLengthComplete}"') && !contigPlasmid[$1$2]++ {print $1$2}' $imageDir/$sample".plasmids.blast" > $imageDir/$sample".dictionary_covered.txt"

awk 'NR==FNR{contigPlasmid[$1]=$1;next} \
{split($1, contigname, "_"); header=$1$2}{if ((header in contigPlasmid) && ($3>90) && (($4/$13)>0.05)) print contigname[length(contigname)], $7,$8,$2,$9,$10, "id=contig_"contigname[length(contigname)]}' \
$imageDir/$sample".dictionary_covered.txt" $imageDir/$sample".plasmids.blast" > $imageDir/$sample".plasmids.blast.links"

#sample.blast.links (file with matching coordinates between coordinates in contigs and plasmids)
#contig_1    915668  915969  NG_048025.1     1178    1476    id=contig_1
#contig_1    123340  123574  NG_048025.1     1361    1127    id=contig_1
#contig_2    1       121     NZ_CP010574.1   66599   66479   id=contig_2
#contig_2    1       121     NZ_CP008930.1   122307  122187  id=contig_2
#contig_2    1       121     NZ_CP011577.1   112628  112748  id=contig_2

##Change coordinates from contig --> plasmid to plasmid-->plasmid in order to represent them in CIRCOSS

awk 'BEGIN{OFS="\t";}{if($1 != savedNode){savedNode= $1; delete chr} else{for(i in chr){print $4" "$5" "$6" "chr[i]" id="savedNode}}chr[$4$5$6] = $4" "$5" "$6}' \
$imageDir/$sample".plasmids.blast.links" > $imageDir/$sample".plasmids.blast.links.sorted"


#NG_048025.1 1361 1127 NG_048025.1 1178 1476 id=contig_1
#NZ_CP008930.1 122307 122187 NZ_CP010574.1 66599 66479 id=contig_2
#NZ_CP011577.1 112628 112748 NZ_CP010574.1 66599 66479 id=contig_2
#NZ_CP011577.1 112628 112748 NZ_CP008930.1 122307 122187 id=contig_2
#NZ_CP006927.1 76749 76629 NZ_CP011577.1 112628 112748 id=contig_2


echo "##### DONE Generating link track for CIRCOS image"


echo "##### Creating annotated track from contigs assembled (annotated with Prokka)"


#Filter Gff(3) file from prokka and create a bed file with coordinates with annotated genes (WITH NAMES)
awk 'BEGIN{OFS="\t";}{split($1,query,"_"); \
split($9,description,"Name=");split(description[2],name,";"); \
split($9,locustag,"locus_tag=");split(locustag[2],locus,";");split(locus[1],locusname,"_");\
split(name[1],nameLowBar,"_")};{if ($3 == "gene" && $1 != "#" && $9 ~ /Name/) {print query[length(query)],$4,$5, name[1]}; \
{if ($3 == "gene" && $1 != "#" && $9 !~ /Name/) {print query[length(query)],$4,$5, "CDS_"locusname[length(locusname)]}}}' \
$imageDir/$sample".gff" > $imageDir/$sample".gff.bed"


#Create a file with matches between file .bed and file .blast.link based on contig name
#$sample.bed (file with annotation from prokka on contigs assembled)
#contig_1    1005    2300    kgtP
#contig_1    2669    4024    pssA
#contig_1    7599    8024    trxC
#contig_1    9371    10060   ung
#contig_1    10373   10756   grcA


echo "##### Adapting coords in annotated contigs"


#create a file with first hit of each blast combination

awk 'NR==FNR{a[NR]=$1;b[NR]=$0;next}{for(i = 1; i <= NR; ++i){if (a[i] == $1) print b[i],"\t", $0}}' \
$imageDir/$sample".gff.bed" $imageDir/$sample".plasmids.blast.links" > $imageDir/$sample".bed.coords.tmp"

#The result file is a combination of annotation and coordinates:
#contig_1    1005    2300    kgtP_1   contig_1   915668  915969  NG_048025.1     1178    1476    id=contig_1
#contig_1    2669    4024    pssA     contig_1   915668  915969  NG_048025.1     1178    1476    id=contig_1
#contig_1    7599    8024    trxC     contig_1   915668  915969  NG_048025.1     1178    1476    id=contig_1
#contig_1    9371    10060   ung      contig_1   915668  915969  NG_048025.1     1178    1476    id=contig_1
#contig_1    10373   10756   grcA     contig_1   915668  915969  NG_048025.1     1178    1476    id=contig_1


awk '(($2 >= $6 - 5000 && $2 <= $7) || ($3 >= $6 && $3 <= $7 + 5000)) {{isInverted=($10-$9); \
genelength=($3-$2)};{if (isInverted < 0) {coordChr1=(($7-$3)+$10);} else {coordChr1=(($2-$6)+$9)}}; \
coordChr2=(coordChr1+genelength); {print $8, coordChr1, coordChr2, $4}}' $imageDir/$sample".bed.coords.tmp" > $imageDir/$sample".bed.all.coords"

#resulting in a bed file with coordinated of plasmid bur refering to contig annotation:
#NZ_CP010574.1 34820 33528 arsB_1
#NZ_CP008930.1 90527 89235 arsB_1
#NZ_CP006927.1 44969 43677 arsB_1
#NZ_CP010574.1 81021 82508 ltrA_1
#NZ_CP008930.1 144220 145707 ltrA_1


#Remove duplicate of several matches 

awk '($2 > 0) && ($3 > 0) && !uniq[$1$4]++' /$imageDir/$sample".bed.all.coords" \
> $imageDir/$sample".bed.name.coords"

awk '{split($4, namelowbar, "_")} {$4=($4 !~ /CDS/) ? namelowbar[1] : $4}1' $imageDir/$sample".bed.name.coords" \
> $imageDir/$sample".bed.coords"

#NZ_CP010574.1 34820 33528 arsB
#NZ_CP008930.1 90527 89235 arsB
#NZ_CP006927.1 44969 43677 arsB
#NZ_CP010574.1 81021 82508 ltrA
#NZ_CP008930.1 144220 145707 ltrA


echo "##### DONE Adapting coords in annotated contigs"
echo "##### DONE Creating annotated track from contigs assembled (annotated with Prokka)"





#blast contigs on ARGannot.fasta DDBB

if [ -f $referenceDir/ARGannot.nhr -a -f $referenceDir/ARGannot.nin -a -f $referenceDir/ARGannot.nsq ];then \
	echo "Found a Blast ddbb for ARGannot";
	echo "Omitting makeblast step"
else

	echo "Making ARGannot Blast DDBB";
	makeblastdb -in $abReference -out $referenceDir/ARGannot -dbtype nucl

fi

echo "##### Creating antibiotic resistance track from contigs assembled"

blastn -query $imageDir/$sample".fna" \
-db $referenceDir/ARGannot \
-out $imageDir/$sample".abr.blast" -evalue 0.0001 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

#/opt/srst2-0.1.8/data/ARGannot.fasta

#gnl|CNM|ECLCNM_contig_1     0__OqxBgb_Flq__OqxBgb__49       73.53   306     70      9       915668  915969  1078    1376    1e-20     106
#gnl|CNM|ECLCNM_contig_1     0__OqxBgb_Flq__OqxBgb__49       73.75   240     53      9       123340  123574  1261    1027    2e-14   86.1
#gnl|CNM|ECLCNM_contig_7     86__AmpH_Bla__AmpH__634 99.22   1161    9       0       62169   63329   1161    1       0.0      2095
#gnl|CNM|ECLCNM_contig_14     371__CatBx_Phe__CatB4__602      100.00  108     0       0       96486   96593   442     549     6e-50     200
#gnl|CNM|ECLCNM_contig_17     0__OqxBgb_Flq__OqxBgb__49       98.86   3153    36      0       31982   35134   1       3153    0.0      5624


#Blast optput format 6 to BED coordinates obtained with prokka (Contig name = gnl|center(CNM)|locustag|contig_{%d}6)
#and matching ARGannot DDBB
#Output = contig name - coordinates on contig - name of ABR gene
#Output is unique since ARGannot DDBB has several serotypes of AR genes

cat $imageDir/$sample".abr.blast" |sort -k 12 -nr |\
awk '{OFS="\t";split($1, contigname, "_"); split($2, ABR, "__"); split(ABR[2], ABRType, "_")};($3 > 95)&&(($4/$14)>0.8)&&(!x[ABR[2]contigname[length(contigname)]]++){print contigname[length(contigname)], $7, $8, ABRType[2]"_"ABR[3]}' \
> $imageDir/$sample".abr.bed"

#contig_17 31982 35134 OqxBgb_Flq
#contig_55 809 1993 TetD_Tet
#contig_17 30783 31958 OqxA_Flq
#contig_7 62169 63329 AmpH_Bla
#contig_36 1073 1948 CTX-M-1_Bla


#Transfors contigs coord to plasmid coords 

#Create a intermediate file with all matches in contig name
awk 'NR==FNR{a[NR]=$1;b[NR]=$0;next}{for(i = 1; i <= NR; ++i){if (a[i] == $1) print b[i],"\t", $0}}' \
$imageDir/$sample".abr.bed" \
$imageDir/$sample".plasmids.blast.links" \
> $imageDir/$sample".abr.bed.coords.tmp"

#contig_7	168907	170064	AMPH_Ecoli_Bla 	 contig_7	74669	75873	NC_019114.1	21471	20266	id=contig_7
#contig_7	168907	170064	AMPH_Ecoli_Bla 	 contig_7	74676	75870	NZ_CP008906.1		139589	140784	id=contig_7
#contig_7	168907	170064	AMPH_Ecoli_Bla 	 contig_7	172420	172671	CP019443.1		26728	26477	id=contig_7
#contig_7	168907	170064	AMPH_Ecoli_Bla 	 contig_7	172420	172671	NZ_CP011062.1	234981	235232	id=contig_7
#contig_8	45602	47227	Mcr1_Colistin 	 contig_8	101753	157404	NZ_CP015833.1	13262	68917	id=contig_8



#If annotation is within blasted fragment of the contig, coordinates are translated from contig to plasmid in order to represent them on CIRCOS


awk '(($2 >= $6 - 5000 && $2 <= $7) || ($3 >= $6 && $3 <= $7 + 5000)) {{isInverted=($10-$9); \
genelength=($3-$2)};{if (isInverted < 0) {coordChr1=(($7-$3)+$10);} else {coordChr1=(($2-$6)+$9)}}; \
coordChr2=(coordChr1+genelength); {print $8, coordChr1, coordChr2, $4}}' $imageDir/$sample".abr.bed.coords.tmp" > $imageDir/$sample".abr.bed.all.coords"


#NG_048025.1 1178 1479 OqxBgb_Flq
#NG_048025.1 1361 1127 OqxBgb_Flq
#NZ_CP010574.1 143152 143045 CatBx_Phe
#NZ_CP008930.1 37212 37105 CatBx_Phe
#NZ_CP011577.1 21393 21286 CatBx_Phe

#Remove duplicate of several matches 

awk '!uniq[$1$4]++'  /$imageDir/$sample".abr.bed.all.coords" \
> $imageDir/$sample".abr.bed.coords"


echo "##### DONE Creating antibiotic resistance track from contigs assembled"

if [ -f $referenceDir/plasmidFinder.nhr -a -f $referenceDir/plasmidFinder.nin -a -f $referenceDir/plasmidFinder.nsq ];then \
	echo "Found a Blast ddbb for PlasmidFinder";
	echo "Omitting makeblast step"
else

	echo "Making PlasmidFinder Blast DDBB";
	makeblastdb -in $plasmidFinderReference -out $referenceDir/plasmidFinder -dbtype nucl

fi

echo "##### Creating replisome track from contigs assembled (PlasmidFinder DDBB)"

#blast contigs on PlasmidFinder DDBB


blastn -query $imageDir/$sample".fna" \
-db $referenceDir/plasmidFinder \
-out $imageDir/$sample".pfinder.blast" -evalue 0.0001 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"


#Obtain bed coords in contigs
#Filter by identity greater than 95%
#Get the first ocurrence of each contig since only one replisome is expected per contig
cat $imageDir/$sample".pfinder.blast" | sort -k 12 -rn \
| awk '{OFS="\t"; split($1,contigname,"_"); split($2,repname,"_")};($3 > 95)&&(($4/$14)>0.8)&&(!contig[repname[1]"_"repname[2]contigname[length(contigname)]]++){print contigname[length(contigname)], $7, $8, repname[1]}' \
> $imageDir/$sample".pfinder.bed"

#contig_107      4433    5114    IncFIB(AP001918)_1
#contig_14       124847  125476  IncHI2A_1
#contig_57       11563   12205   IncFIB(S)_1
#contig_52       46975   47301   IncHI2_1
#contig_116      122     389     IncFIA_1

#Transfors contigs coord to plasmid coords 

#Create a intermediate file with all matches in contig name
awk 'NR==FNR{a[NR]=$1;b[NR]=$0;next}{for(i = 1; i <= NR; ++i){if (a[i] == $1) print b[i],"\t", $0}}' \
$imageDir/$sample".pfinder.bed" \
$imageDir/$sample".plasmids.blast.links" \
> $imageDir/$sample".pfinder.bed.coords.tmp"


#contig_25	59506	59653	IncFII(K)_1 	 contig_25	12675	62112	NZ_CP015132.1	92402	42983	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	6277	11583	NZ_CP015132.1	97699	92393	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	2947	6212	NZ_CP015132.1	101088	97819	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	2	2826	NZ_CP015132.1	103977	101153	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	16775	17556	NZ_CP015132.1	98354	97573	id=contig_25


#If annotation is within blasted fragment of the contig, coordinates are translated from contig to plasmid in order to represent them on CIRCOS


awk '(($2 >= $6 - 2000 && $2 <= $7) || ($3 >= $6 && $3 <= $7 + 2000)) {{isInverted=($10-$9); \
genelength=($3-$2)};{if (isInverted < 0) {coordChr1=(($7-$3)+$10);} else {coordChr1=(($2-$6)+$9)}}; \
coordChr2=(coordChr1+genelength); {print $8, coordChr1, coordChr2, $4}}' $imageDir/$sample".pfinder.bed.coords.tmp" \
> $imageDir/$sample".pfinder.bed.coords"


#NZ_CP015132.1 139233 139086 IncFII(K)_1
#NZ_CP012884.1 125294 125441 IncFII(K)_1
#NZ_CP018365.1 207768 207915 IncFII(K)_1
#NZ_CP018430.1 165980 166127 IncFII(K)_1
#NC_009649.1 5612 5759 IncFII(K)_1


echo "##### DONE Creating replisome track from contigs assembled (PlasmidFinder DDBB)"

echo "Done with data generation"


#blast contigs on ISFinder DDBB
<<ISFINDER

echo "##### Creating IS track from contigs assembled (ISFinder DDBB)"
makeblastdb -in $isFinderReference -out $clusterDir/REFERENCES/PLASMIDS/ISFINDER -dbtype nucl

blastn -word_size 10 -query $imageDir/$sample".fna" \
-db $clusterDir/REFERENCES/PLASMIDS/ISFINDER \
-out $imageDir/$sample".isfinder.blast" -evalue 0.0001 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"


#Obtain bed coords in contigs
#Filter by identity greater than 95%
#Get the first ocurrence of each contig since only one replisome is expected per contig
cat $imageDir/$sample".isfinder.blast" | sort -k 12 -rn \
| awk '{OFS="\t"; split($1,contigname,"_"); split($2,isname,"_")};($3 > 90)&&(($4/$14)>0.8)&&(!contig[isname[1]"_"isname[2]contigname[length(contigname)]]++){print contigname[length(contigname)], $7, $8, $2}' \
> $imageDir/$sample".isfinder.bed"

#| awk '{OFS="\t"; split($1,contigname,"_"); split($2,repname,"_")};($3 > 95)&&(($4)>500)&&(!contig[isname[1]"_"isname[2]contigname[length(contigname)]]++){print "contig_"contigname[length(contigname)], $7, $8, $2}' \


#contig_107      4433    5114    IncFIB(AP001918)_1
#contig_14       124847  125476  IncHI2A_1
#contig_57       11563   12205   IncFIB(S)_1
#contig_52       46975   47301   IncHI2_1
#contig_116      122     389     IncFIA_1

#Transfors contigs coord to plasmid coords 

#Create a intermediate file with all matches in contig name
awk 'NR==FNR{a[NR]=$1;b[NR]=$0;next}{for(i = 1; i <= NR; ++i){if (a[i] == $1) print b[i],"\t", $0}}' \
$imageDir/$sample".isfinder.bed" \
$imageDir/$sample".plasmids.blast.links" \
> $imageDir/$sample".isfinder.bed.coords.tmp"


#contig_25	59506	59653	IncFII(K)_1 	 contig_25	12675	62112	NZ_CP015132.1	92402	42983	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	6277	11583	NZ_CP015132.1	97699	92393	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	2947	6212	NZ_CP015132.1	101088	97819	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	2	2826	NZ_CP015132.1	103977	101153	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	16775	17556	NZ_CP015132.1	98354	97573	id=contig_25

#If annotation is within blasted fragment of the contig, coordinates are translated from contig to plasmid in order to represent them on CIRCOS


awk '(($2 >= $6 - 5000 && $2 <= $7) || ($3 >= $6 && $3 <= $7 + 5000)) {{isInverted=($10-$9); \
genelength=($3-$2)};{if (isInverted < 0) {coordChr1=(($7-$3)+$10);} else {coordChr1=(($2-$6)+$9)}}; \
coordChr2=(coordChr1+genelength); {print $8, coordChr1, coordChr2, $4}}' $imageDir/$sample".isfinder.bed.coords.tmp" \
> $imageDir/$sample".isfinder.bed.coords"

rm $imageDir/$sample".isfinder.bed.coords.tmp"

#NZ_CP015132.1 139233 139086 IncFII(K)_1
#NZ_CP012884.1 125294 125441 IncFII(K)_1
#NZ_CP018365.1 207768 207915 IncFII(K)_1
#NZ_CP018430.1 165980 166127 IncFII(K)_1
#NC_009649.1 5612 5759 IncFII(K)_1

echo "##### DONE Creating IS track from contigs assembled (ISFinder DDBB)"

ISFINDER

<<C
rm $imageDir/$sample".plasmids.blast.links"
rm $imageDir/$sample".bed.coords.tmp"
rm $imageDir/$sample".bed.all.coords"
rm $imageDir/$sample".abr.bed.coords.tmp"
rm $imageDir/$sample".abr.bed.all.coords"
rm $imageDir/$sample".pfinder.bed.coords.tmp"
rm $imageDir/$sample".plasmid.bedgraph"
C