#!/bin/bash

<<Usage
#Copy this command to generate a file (commands_placnet_by_sample.txt) with the order to execute this script for each sample


group=NLAB
samplesid=/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/DOC/SAMPLES/samplesid_$group.txt

cat $samplesid | while read -r line
do echo "qsub -V -b y -j y -cwd -q all.q -N STS_${line}_PLSMD bash 5_Stats_coverage.sh $line $group"
done > 5_commands_stats_coverage_$group.txt

bash 5_commands_stats_coverage_$group.txt
Usage

###############1 Map reads on plasmid database & filter positive matches
###############2 Find plasmid with more coverage
###############3 Draw plasmid with more coverage (>90%)


#IN: List of ac mapped more than 90% in plasmid DDBB clustered with CCD-HIT (80%)


#sample=$1
group=$1
samplesid=/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/DOC/SAMPLES/samplesid_$group.txt
plasmidDdbb="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/plasmid.all.genomic.feb212017.fasta"
trimmedDir="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/ANALYSIS/TRIMMED"
mappedDir="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/ANALYSIS/MAPPING/PLASMIDS"
imageDir="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/ANALYSIS/IMAGES/PLASMIDS"
contigDir="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/ANALYSIS/ASSEMBLY"
#plasmidLenght="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/NCBI_PLASMID_LENGTH.tsv"
plasmidLenght="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/plasmid.all.genomic.feb212017.length"

plasmidMapped=$mappedDir/$group/$sample/$sample".ac.covered.txt"
plasmidMappedFasta=$mappedDir/$group/$sample/$sample".ac.covered.fasta"


echo "obtaining list and description of all plasmids matching requisites in $group"

#check full description of each plasmid matching the requeriments
#Obtain full description of plasmids present in all samples in one file
for sampleid in $(cat $samplesid);do \
  for ac in $(cat $mappedDir/$group/$sampleid/$sampleid".ac.covered.txt_80");do \
    awk 'split(FILENAME,sample,"/") {if($0 ~ /'$ac'/) {print sample[length(sample)-1],"\t", $0}}' $mappedDir/$group/$sampleid/$sampleid".ac.covered.fasta_80_80" 
  done
done > $imageDir/$group/sample_description_list_$group.txt

#A14_03    >NZ_CP007733.1 Klebsiella pneumoniae subsp. pneumoniae KPNIH27 plasmid pKPN-068, complete sequence
#A14_03   >NC_009649.1 Klebsiella pneumoniae subsp. pneumoniae MGH 78578 plasmid pKPN3, complete sequence
#A14_03   >NZ_CP011990.1 Klebsiella pneumoniae UHKPC33 plasmid pUHKPC33-162.533kb, complete sequence
#A14_03   >NZ_CP014698.2 Klebsiella quasipneumoniae strain ATCC 700603 plasmid pKQPS2, complete sequence
#A14_03   >NZ_CP011640.1 Serratia marcescens strain CAV1492 plasmid pCAV1492-73, complete sequence
#C01_03   >NZ_CM004622.1 Escherichia coli strain Ec47VL plasmid pEC47a, whole genome shotgun sequence
#C01_03   >NZ_CP011334.1 Escherichia coli O104:H4 str. C227-11 plasmid, complete sequence
#C01_03   >NZ_CP010390.1 Klebsiella pneumoniae strain 6234 plasmid p6234-198.371kb, complete sequence
#C01_03   >NZ_CP018443.1 Klebsiella pneumoniae strain Kp_Goe_822917 plasmid pKp_Goe_917-2, complete sequence
#C01_03   >NZ_CP018694.1 Klebsiella pneumoniae strain Kp_Goe_821588 plasmid pKp_Goe_588-2, complete sequence 


#Obtain a list of all common plasmids present in all samples within the group
#split name and retrieve AC number (NX_######)
#create an array with no overlap
#print array

echo "obtaining list of common plamsids in all samples $group"

awk '{gsub(">","")} {uniq[$2]++} END {for (i in uniq) print i}' $imageDir/$group/sample_description_list_$group.txt | sort > $imageDir/$group/sample_ac_common_list_$group.txt

#NZ_CP006801.1                                                                                                                                                        
#NC_019089.1                                                                                                                                                          
#NC_023027.1                                                                                                                                                          
#NC_019251.1                                                                                                                                                          
#NZ_CP009857.1                                                                                                                                                        
#NC_020552.1                                                                                                                                                          
#NZ_CP011990.1                                                                                                                                                        
#NZ_CP011977.1                                                                                                                                                        
#NC_025184.1                                                                                                                                                          
#NC_007414.1                                                                                                                                                          
#NC_004833.1

echo "obtaining percentages of common plamsids in all samples $group"

for i in $(find $mappedDir/$group -type f -name "*_percentage*"); do \
awk 'BEGIN{OFS="\t"} split(FILENAME,sample,"/") {print sample[length(sample)-1],$1,$2}' $i ;
done > $imageDir/$group/sample_ac_percentage_$group.txt


#MPVECO8516_S58  NC_017630.1     84.0323                                                                                                                                                               
#MPVECO8516_S58  NC_011749.1     87.3656                                                                                                                                                               
#MPVECO8516_S58  NC_021363.1     91.4132                                                                                                                                                               
#MPVECO8516_S58  NC_021364.1     96.1917                                                                                                                                                               
#MPVECO8516_S58  NC_016840.1     92.2687                                                                                                                                                               
#MPVECO8516_S58  NC_010421.1     91.3282                                                                                                                                                               
#MPVECO8516_S58  NZ_CP017726.1   81.3283                                                                                                                                                               
#MPVECO8516_S58  NC_024986.1     80.4921                                                                                                                                                               
#MPVECO8517_S44  NC_020280.1     95.1215                                                                                                                                                               
#MPVECO8517_S44  NZ_CP008906.1   81.6608

echo "obtaining length of common plamsids in all samples $group"

if [ -f $imageDir/$group/sample_ac_common_lenght_$group.txt ];then \
	rm $imageDir/$group/sample_ac_common_lenght_$group.txt;
	rm $imageDir/$group/sample_ac_common_join_$group.txt;
fi

for ac in $(cat $imageDir/$group/sample_ac_common_list_$group.txt); do \
awk '/'$ac'/ {print $0}' $plasmidLenght >> $imageDir/$group/sample_ac_common_lenght_$group.txt ;
awk '/'$ac'/ {print "'$ac'", $0}' $imageDir/$group/sample_description_list_$group.txt >> $imageDir/$group/sample_ac_common_join_$group.txt;
done

join $imageDir/$group/sample_ac_common_lenght_$group.txt $imageDir/$group/sample_ac_common_join_$group.txt > $imageDir/$group/ac_length_sample_description_$group.txt


#NC_005246.1 60145 A14_03 >NC_005246.1 Erwinia amylovora LebB66 plasmid pEL60, complete sequence
#NC_005246.1 60145 G00_07 >NC_005246.1 Erwinia amylovora LebB66 plasmid pEL60, complete sequence
#NC_005246.1 60145 K3087 >NC_005246.1 Erwinia amylovora LebB66 plasmid pEL60, complete sequence
#NC_005246.1 60145 K3167 >NC_005246.1 Erwinia amylovora LebB66 plasmid pEL60, complete sequence
#NC_009649.1 175879 A14_03 >NC_009649.1 Klebsiella pneumoniae subsp. pneumoniae MGH 78578 plasmid pKPN3, complete sequence
#NC_009649.1 175879 K3167 >NC_009649.1 Klebsiella pneumoniae subsp. pneumoniae MGH 78578 plasmid pKPN3, complete sequence
#NC_011410.1 12656 C01_03 >NC_011410.1 Salmonella enterica subsp. enterica serovar Typhimurium plasmid pFPTB1 delta tnpA, ORF294, tetA, tetR, blaTEM, tnpR and tnp genes
#NC_011410.1 12656 G00_07 >NC_011410.1 Salmonella enterica subsp. enterica serovar Typhimurium plasmid pFPTB1 delta tnpA, ORF294, tetA, tetR, blaTEM, tnpR and tnp genes
#NC_011410.1 12656 K3087 >NC_011410.1 Salmonella enterica subsp. enterica serovar Typhimurium plasmid pFPTB1 delta tnpA, ORF294, tetA, tetR, blaTEM, tnpR and tnp genes
#NC_011640.1 4211 K3087 >NC_011640.1 Klebsiella pneumoniae plasmid pKpn114, complete sequence
#awk '{print substr($0, index($0,$3))}'


echo "preparing table with plasmid description $group"

awk '{gsub(/complete sequence/,"") gsub(",","")} {print $1, "\t", $2, "\t", $3, "\t", $5, $6, "\t",substr($0, index($0,$7))}' \
$imageDir/$group/ac_length_sample_description_$group.txt | sort > $imageDir/$group/table_ac_length_sample_description_$group.tsv

awk '{gsub(/complete sequence/,"") gsub(",","")} {print $1, "\t", $2, "\t", $5, $6, "\t", substr($0, index($0,$7))}' \
$imageDir/$group/ac_length_sample_description_$group.txt | sort -u | cut -d "," -f 1 > $imageDir/$group/table_ac_length_NO_sample_description_$group.tsv

#NZ_CP011990.1    162533          3060_S7         Klebsiella pneumoniae   UHKPC33 plasmid pUHKPC33-162.533kb, NZ_CP011990.1 
#NZ_CP011990.1    162533          3594_S5         Klebsiella pneumoniae   UHKPC33 plasmid pUHKPC33-162.533kb, NZ_CP011990.1 
#NZ_CP011990.1    162533          3734_S4         Klebsiella pneumoniae   UHKPC33 plasmid pUHKPC33-162.533kb, NZ_CP011990.1 
#NZ_CP011990.1    162533          3950_S3         Klebsiella pneumoniae   UHKPC33 plasmid pUHKPC33-162.533kb, NZ_CP011990.1 
#NZ_CP011990.1    162533          K5496_S4        Klebsiella pneumoniae   UHKPC33 plasmid pUHKPC33-162.533kb, NZ_CP011990.1 
#NZ_CP011990.1    162533          K5761_S8        Klebsiella pneumoniae   UHKPC33 plasmid pUHKPC33-162.533kb, NZ_CP011990.1 
#NZ_CP012803.1    92739   2013-16_S18     Escherichia coli        O157:H7 strain WS4202 NZ_CP012803.1 plasmid pO157-WS4202
#NZ_CP014006.1    46161   K5761_S8        Klebsiella pneumoniae   subsp. pneumoniae strain NZ_CP014006.1 NUHL24835 plasmid unnamed2
#NZ_CP015021.1    81401   2013-16_S18     Escherichia coli        strain 28RC1 plasmid NZ_CP015021.1 p28RC1, 
#NZ_CP015022.1    95170   2013-16_S18     Escherichia coli        strain SRCC 1675 NZ_CP015022.1
#NZ_CP015075.1    63543   K5050_S6        Escherichia coli        strain Ecol_745 plasmid NZ_CP015075.1 pEC745_OXA48
#NZ_CP015133.1    26450   3061_S6         Klebsiella pneumoniae   strain Kpn555 plasmid NZ_CP015133.1 pKPN-d6b 
#NZ_CP015133.1    26450   9517_7#47       Klebsiella pneumoniae   strain Kpn555 plasmid NZ_CP015133.1 pKPN-d6b
#NZ_CP015133.1    26450   9517_7#48       Klebsiella pneumoniae   strain Kpn555 plasmid NZ_CP015133.1 pKPN-d6b
#NZ_CP015133.1    26450   K5410_S7        Klebsiella pneumoniae   strain Kpn555 plasmid NZ_CP015133.1 pKPN-d6b


echo "obtaining table with presence/ausence $group"

#create a file with all combinations 

for sampleid in $(cat /processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/DOC/SAMPLES/samplesid_$group.txt);do \
  for ac in $(cat $imageDir/$group/sample_ac_common_list_$group.txt); do \
    echo -e $sampleid "\t" $ac
  done
done > $imageDir/$group/sample_ac_combination_list_$group.txt

#2013-16_S18    NC_023027.1
#2013-16_S18    NC_023314.1
#2013-16_S18    NC_024992.1
#2013-16_S18    NC_025007.1
#2013-16_S18    NC_025009.1
#2013-16_S18    NC_025105.1
#2013-16_S18    NC_025130.1
#2013-16_S18    NC_025184.1
#2013-16_S18    NZ_CP006800.1
#2013-16_S18    NZ_CP006801.1
#2013-16_S18    NZ_CP006924.1
#2013-16_S18    NZ_CP006927.1
#2013-16_S18    NZ_CP008721.1



#awk 'NR==FNR{comb[$1" "$2];next}{table[$3" "$1]}END{for (i in comb) {for (j in table) {if (i ~ j) print i, 1; else print i, 0}' #previous syntax same result

awk 'NR==FNR{comb[$1" "$2];next}{table[$3" "$1]}END{for (i in comb) {for (j in table) {if (i ~ j) {print i, 1;} else {print i, 0}}}}' \
$imageDir/$group/sample_ac_combination_list_$group.txt $imageDir/$group/table_ac_length_sample_description_$group.tsv | sort -u \
> $imageDir/$group/table_ac_sample_presence_$group.tsv

#2013-16_S18 NZ_LN824135.1 0
#2013-16_S18 NZ_LN868944.1 0
#2013-16_S18 NZ_LN868944.1 1
#2013-16_S18 NZ_LN868945.1 0
#2013-16_S18 NZ_LN868945.1 1
#2013-16_S18 NZ_LN868946.1 0
#2013-16_S18 NZ_LN868946.1 1
#2719_S10 NC_002145.1 0
#2719_S10 NC_004833.1 0
#2719_S10 NC_004989.1 0
#2719_S10 NC_006855.1 0



awk 'NR==FNR{comb[$1" "$2];next}{table[$1" "$2]}{percentage[$1" "$2]=$3}END{for (i in comb) {for (j in table) {if (i ~ j) {print i, percentage[i];} else {print i, 0}}}}' \
$imageDir/$group/sample_ac_combination_list_$group.txt $imageDir/$group/sample_ac_percentage_$group.txt | sort -u \
> $imageDir/$group/table_ac_sample_percentage_$group.tsv

#MPVECO111_S45 NC_025179.1 0
#MPVECO111_S45 NC_025179.1 80.9484
#MPVECO111_S45 NC_025180.1 0
#MPVECO111_S45 NC_025186.1 0
#MPVECO111_S45 NC_025198.1 0
#MPVECO111_S45 NZ_AFET01000005.1 0
#MPVECO111_S45 NZ_AFET01000005.1 92.0264
#MPVECO111_S45 NZ_AFVX01000096.1 0
#MPVECO111_S45 NZ_AFVX01000096.1 82.4395
#MPVECO111_S45 NZ_AGTD01000002.1 0

#Remove false negatives by merging them with true positives and format output as table



echo "sorting presence/ausence table $group"


columnNumber=$(cat $samplesid | wc -l)
rowNumber=$(cat $imageDir/$group/sample_ac_common_list_$group.txt | wc -l )

awk '{a[$1" "$2]++};END{for (i in a) print i, a[i]}' $imageDir/$group/table_ac_sample_presence_$group.tsv | sort > $imageDir/$group/table_ac_sample_presence_fixed_$group.tsv


awk -v nr=$rowNumber '{gsub("1","0") gsub("2","1")} {a[NR]=$3} END{for (i=1;i<=nr;i++) {for (j=i;j<=NR;j+=nr) printf "%s\t",a[j]; print""}} END{gsub("1","0"); gsub("2","1")}' \
$imageDir/$group/table_ac_sample_presence_fixed_$group.tsv > $imageDir/$group/table_presence.plasmids_$group.tsv

#awk -v n=$columnNumber 'BEGIN { row = row $3 "\t"; if (NR % n == 0) { print row; row = "" } }' $imageDir/$group/table_ac_sample_presence_fixed_$group.tsv > $imageDir/$group/table_presence.plasmids_$group.tsv

#1 1 2 2 2 1 1 2 1 1 1 1 1 1 1 1 1 1 1 2 1 1                                                                                                                                                                                                  
#1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2                                                                                                                                                                                                  
#1 1 1 1 1 1 1 1 1 2 1 2 2 1 1 1 2 2 1 2 2 2                                                                                                                                                                                                  
#1 1 1 1 1 2 1 1 1 1 1 1 1 1 2 2 1 1 1 1 2 1                                                                                                                                                                                                  
#1 1 1 1 2 1 1 2 2 1 2 1 2 1 1 2 2 2 1 2 1 1                                                                                                                                                                                                  
#2 2 1 1 1 1 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1                                                                                                                                                                                                  
#1 1 1 1 1 2 1 1 1 1 1 1 1 1 2 2 1 1 1 1 2 1                                                                                                                                                                                                  
#1 1 1 1 1 1 1 1 2 1 2 1 1 1 1 2 2 1 1 2 1 1                                                                                                                                                                                                  
#2 2 1 1 2 1 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1                                                                                                                                                                                                  
#1 1 1 1 1 2 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1



echo "sorting percentage table $group"

sort -k 3 -nr $imageDir/$group/table_ac_sample_percentage_$group.tsv | awk '(!first[$1$2]++)' | sort > $imageDir/$group/table_ac_sample_percentage_fixed_$group.tsv

awk -v nr=$rowNumber '{a[NR]=$3} END{for (i=1;i<=nr;i++) {for (j=i;j<=NR;j+=nr) printf "%s\t",a[j]; print""}}' \
$imageDir/$group/table_ac_sample_percentage_fixed_$group.tsv > $imageDir/$group/table_percentage.plasmids_$group.tsv


#echo "merging all data in a final report table $group"

#paste <(awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' $imageDir/$group/table_presence.plasmids_$group.tsv) <(cat $imageDir/$group/table_presence.plasmids_$group.tsv)  > $imageDir/$group/table_presence.plasmids_$group.tsv

#printf "%s\t" AC_Number Length Species Description N $(cat $samplesid | sort) > $imageDir/$group/header_column_$group.tsv

#printf "\n" >> $imageDir/$group/header_column_$group.tsv

#paste $imageDir/$group/table_ac_length_NO_sample_description_$group.tsv $imageDir/$group/table_presence.plasmids_$group.tsv > $imageDir/$group/FINAL_REPORT_NO_samples_$group.tsv

#cat $imageDir/$group/header_column_$group.tsv $imageDir/$group/FINAL_REPORT_NO_samples_$group.tsv > $imageDir/$group/FINAL_REPORT_$group.tsv



echo "merging all data in a final report table $group"

paste <(awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' $imageDir/$group/table_presence.plasmids_$group.tsv) <(cat $imageDir/$group/table_percentage.plasmids_$group.tsv)  > $imageDir/$group/table_report.plasmids_$group.tsv

printf "%s\t" AC_Number Length Species Description N $(cat $samplesid | sort) > $imageDir/$group/header_column_$group.tsv

printf "\n" >> $imageDir/$group/header_column_$group.tsv

paste $imageDir/$group/table_ac_length_NO_sample_description_$group.tsv $imageDir/$group/table_report.plasmids_$group.tsv > $imageDir/$group/FINAL_REPORT_NO_samples_$group.tsv

cat $imageDir/$group/header_column_$group.tsv $imageDir/$group/FINAL_REPORT_NO_samples_$group.tsv > $imageDir/$group/FINAL_TMP_$group.tsv



#PROCESS CLUST FILE

echo "processing cluster file to append in a final report table $group"

cat $(find $mappedDir/$group -name "*.clstr") | \
awk 'BEGIN{RS=">Cluster "} {gsub("\n",""); gsub("\t"," "); gsub(/\.\.\. at/, "  ")} {split($0,ac,", >")} ($0 ~ "%") {printf "%s\n", $0}' | \
awk 'BEGIN {FS=">"} {gsub(/[0-9]* *[0-9]*aa, >/," "); gsub(/[0-9]*.[0-9]*%/,""); gsub("  "," ")} {print}' | \
sort -u > $imageDir/$group/selected_$group.txt

#CP015397.1  NZ_CP006924.1... *
#CP017282.1... * CP017288.1  
#NC_003079.1... * NC_019102.1  
#NC_003114.1  NC_008444.1  NC_019987.1... *
#NC_003457.1  NC_011214.1... *
#NC_004989.1... * U59131.1  
#NC_005248.1  NC_022376.1... *
#NC_005324.1... * NC_011795.1  NC_008597.1  NC_019060.1  NC_011602.1  NC_017658.1  NZ_CP006639.1  
#NC_005327.1... * CP019444.1  
#NC_005970.1... * NC_016903.1  NC_017661.1  NC_017636.1  NC_017662.1

cat $imageDir/$group/selected_$group.txt | awk '{match($0,/^.*... \*/,a);split(a[0],name," ");gsub(/\.\.\.( | \*)/," "); gsub(/\.\.\./,"",name[length(name)-1])} {print name[length(name)-1], $0}' | \
awk '{a[$1]=a[$1]" "$0}END{for (i in a) print i a[i]}' > $imageDir/$group/representative_ac_appended_$group.txt


#NC_018652.1 NC_018652.1  NC_019059.1  NC_018652.1 
#NZ_CP010149.1 NZ_CP010149.1  NZ_CP010149.1  NZ_CP006635.1  
#NZ_CP013221.1 NZ_CP013221.1  NZ_CP015996.1  NZ_CP013221.1  NC_021813.2  
#NC_014231.1 NC_014231.1  NC_019087.1  NC_014231.1 
#NC_019058.1 NC_019058.1  NC_019058.1  NC_019053.1 

awk '{ while(++i<=NF) printf (!a[$i]++) ? $i FS : ""; i=split("",a); print "" }' $imageDir/$group/representative_ac_appended_$group.txt > $imageDir/$group/representative_nodupes_$group.txt

#NC_018652.1 NC_019059.1 
#NZ_CP010149.1 NZ_CP006635.1 
#NZ_CP013221.1 NZ_CP015996.1 NC_021813.2 
#NC_014231.1 NC_019087.1 
#NC_019058.1 NC_019053.1 
#NC_003486.1 NC_015515.1 
#NC_019040.1 NZ_CP006639.1 
#NZ_CP010199.1 NZ_CP010224.1 
#NC_006815.1 NZ_CP007488.1 
#NC_019043.1 NC_017718.1 NC_017675.1

#This command remove dupes but doesn't maintain the $1 intact
#awk '{delete a; for (i=1; i<=NF; i++) a[$i]++; for (i in a) printf i" "; print "" }' $imageDir/$group/representative_ac_appended_$group.txt > $imageDir/$group/representative_nodupes_$group.txt

awk 'NR==FNR{representante[$1]=$1;similares[$1]=substr($0,index($0,$2),length($0));next}{if ($1 == representante[$1]) {print $0, "\t",similares[$1];} else {print $0}}' \
$imageDir/$group/representative_nodupes_$group.txt $imageDir/$group/FINAL_TMP_$group.tsv > $imageDir/$group/FINAL_REPORT_$group.tsv

#AC_Number		Length	Species					Description						N		MPVECO111_S45	MPVECO7660_S55	MPVECO7676_S56	MPVECO7840_S47	MPVECO7896_S53	MPVECO8037_S51	MPVECO8183_S52	MPVECO8514_S57	MPVECO8516_S58	MPVECO8517_S44	MPVECO8557_S82	MPVK8386_S46	
#CP017282.1 	 58666 	 Klebsiella variicola 	 strain GJ1 plasmid pKPGJ-1c 	3	0	93.0999	0	94.1397	94.0681	0	0	0	0	0	0	0	 	 CP017288.1 
#CP017853.1 	 58666 	 Klebsiella variicola 	 strain GJ2 plasmid pKPGJ-2d 	3	0	93.0931	0	94.1414	94.0681	0	0	0	0	0	0	0	
#CP019214.1 	 289255	 Escherichia coli 	 strain WCHEC1613 plasmid pMCR_WCHEC1613 	2	0	0	0	0	0	0	0	81.7438	0	82.7758	0	0	
#CP019215.1 	 84697 	 Escherichia coli 	 strain WCHEC1613 plasmid pI1_WCHEC1613 	2	0	0	0	0	0	92.0812	0	0	0	85.7882	0	0	
#CP019443.1 	 233802  Salmonella enterica 	 subsp. enterica serovar Typhimurium strain 81741 plasmid unnamed1 	4	0	85.1528	87.4377	0	0	0	0	81.8004	0	81.7243	0	0	
#NC_002525.1 	 75582 	 Escherichia coli 	 K-12 plasmid R721 	1	0	0	0	95.0941	0	0	0	0	0	0	0	0	
#NC_004989.1 	 2120 	 Serratia marcescens 	 plasmid R478 	2	0	99.9528	0	0	0	0	0	0	0	99.9528	0	0	 	 U59131.1 
#NC_005246.1 	 60145 	 Erwinia amylovora 	 LebB66 plasmid pEL60 	3	0	83.1823	0	84.1699	84.2946	0	0	0	0	0	0	0	
#NC_006671.1 	 101375  Escherichia coli 	 A2363 plasmid pAPEC-O2-R 	1	0	0	0	0	0	0	0	0	0	85.2587	0	0	




echo "All done with $group , Check FINAL_REPORT_$group.tsv"



<<C
#extract all fields but $1
#awk '{sim=substr($0,index($0,$2),length($0));print sim}'

#cat dolar1.txt | awk '{a[$1]=a[$1]" "$0}END{for (i in a) print i a[i]}' > appended.txt
#cat $(find ../../../MAPPING/PLASMIDS/NLAB -name "*.clstr") | \
#awk 'BEGIN{RS=">Cluster "} {gsub("\n",""); gsub("\t"," "); gsub(/\.\.\. (\*|at)/, " ")} {split($0,ac,", >")} ($0 ~ "%") {printf "%s\n", $0}' | \
#awk 'BEGIN {FS=">"} {gsub(/[0-9]* *[0-9]*aa, >/,""); gsub(/[0-9]*.[0-9]*%/,""); gsub("  "," ")} {print}' | sort -u


#awk 'BEGIN{RS=" "};{gsub("\n","")}{printf "%s\n", $0}' | awk 'dupes[$1]++'


#cat MPVECO7676_S56.ac.covered.fasta_80.clstr| awk 'BEGIN{RS=">Cluster "} {gsub("\n",""); gsub("\t"," "); gsub(/\.\.\. (\*|at)/, " ")} {split($0,ac,", >")} ($0 ~ "%") {printf "%s\n", $0}' > MPVECO7676_S56.ac.covered.fasta_80.clstr.line

#cat MPVECO7676_S56.ac.covered.fasta_80.clstr.line | awk 'BEGIN {FS=">"} {gsub(/[0-9]* *[0-9]*aa, >/,""); gsub(/[0-9]*.[0-9]*%/,"")} {print}'
C