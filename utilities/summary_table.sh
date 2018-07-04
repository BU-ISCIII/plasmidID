#!/bin/bash

#set -e

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.0 
#CREATED: 12 June 2018
#REVISION:
#DESCRIPTION:Script creates a summary table fir plasmiID mapping results
#
#
#================================================================
# END_OF_HEADER
#================================================================



usage() {
	cat << EOF

summary table is a script that  creates a summary table fir plasmiID mapping results

usage : $0 <pID_group>


Output directory is placed in group folder, named 00-summary

example: ./summary_table.sh ECO_NDM 
		 

EOF
}

if [ $# = 0 ] ; then
 usage >&2
 exit 1
fi

group=$1
plasmidIdDir="$(pwd)"
mappedDir=$plasmidIdDir/$group/$sample/mapping
summaryDir=$plasmidIdDir/$group/00_summary

mkdir -p $summaryDir


if [ -f $plasmidIdDir/$group/00_summary/sampleid_found_$group.txt ];then \
	
	rm $plasmidIdDir/$group/00_summary/sampleid_found_$group.txt
fi

if [ -f $plasmidIdDir/$group/00_summary/sampleid_not_found_$group.txt ]; then \

	rm $plasmidIdDir/$group/00_summary/sampleid_not_found_$group.txt
	
fi


echo "checking samples within the group"

for sample_in_group in $(ls $plasmidIdDir/$group | grep -v 00_summary); do \

	final_mapping_file=$(find $plasmidIdDir/$group/$sample_in_group/mapping -type f -name "*_adapted_filtered_*" | awk '/_term.fasta_..$/' | sort | awk 'NR==1' | wc -l)

	if [ $final_mapping_file -gt 0 ]; then \

		echo $sample_in_group >> $plasmidIdDir/$group/00_summary/sampleid_found_$group.txt
	
	else 
	
		echo $sample_in_group >> $plasmidIdDir/$group/00_summary/sampleid_not_found_$group.txt
	fi

done

if [ -f $plasmidIdDir/$group/00_summary/sampleid_not_found_$group.txt ]; then \

	echo "WARNING"
	printf "%s," $(cat $plasmidIdDir/$group/00_summary/sampleid_not_found_$group.txt) | sed 's/,$//g'
	echo -e "\nNot found, using the rest for summary table"
	
fi



echo "obtaining list and description of all plasmids matching requisites in $group"

#check full description of each plasmid matching the requeriments
#Obtain full description of plasmids present in all samples in one file

awk '{split(FILENAME,sample,"/")} />/ {print sample[length(sample)-2],"\t",$0}' $plasmidIdDir/$group/*/mapping/*coverage_adapted_filtered_*_term.fasta_[0-9][0-9] > $summaryDir/sample_description_list_$group.txt

#A14_03   >NZ_CP007733.1 Klebsiella pneumoniae subsp. pneumoniae KPNIH27 plasmid pKPN-068, complete sequence
#A14_03   >NC_009649.1 Klebsiella pneumoniae subsp. pneumoniae MGH 78578 plasmid pKPN3, complete sequence
#A14_03   >NZ_CP011990.1 Klebsiella pneumoniae UHKPC33 plasmid pUHKPC33-162.533kb, complete sequence
#C01_03   >NZ_CP018443.1 Klebsiella pneumoniae strain Kp_Goe_822917 plasmid pKp_Goe_917-2, complete sequence
#C01_03   >NZ_CP018694.1 Klebsiella pneumoniae strain Kp_Goe_821588 plasmid pKp_Goe_588-2, complete sequence 


#Obtain a list of all common plasmids present in all samples within the group
#create an array with no overlap
#print array

echo "obtaining list of common plamsids in all samples $group"

awk '{gsub(">","")} {uniq[$2]++} END {for (i in uniq) print i}' $summaryDir/sample_description_list_$group.txt | sort > $summaryDir/sample_ac_common_list_$group.txt

#NZ_CP006801.1                                                                                                                                                        
#NC_019089.1                                                                                                                                                          
#NC_023027.1                                                                                                                                                          
#NC_019251.1                                                                                                                                                          
#NZ_CP009857.1                                                                                                                                                        

echo "obtaining percentages of common plasmids in all samples $group"

for i in $(find $plasmidIdDir/$group/ -type f -name "*_percentage"); do \
	awk 'BEGIN{OFS="\t"} split(FILENAME,sample,"/") {print sample[length(sample)-2],$1,$2}' $i ;
done > $summaryDir/sample_ac_percentage_$group.txt


#MPVECO8516_S58  NC_017630.1     84.0323                                                                                                                                                               
#MPVECO8516_S58  NC_011749.1     87.3656                                                                                                                                                               
#MPVECO8516_S58  NC_021363.1     91.4132                                                                                                                                                               
#MPVECO8516_S58  NC_021364.1     96.1917                                                                                                                                                               
#MPVECO8516_S58  NC_016840.1     92.2687                                                                                                                                                               

echo "obtaining length of common plamsids in all samples $group"


awk 'BEGIN {FS=="| "} /^>/ {if (seqlen) print seqlen;printf "%s\t", $1; seqlen=0; next} {seqlen+=length($0)} END {print seqlen}' $plasmidIdDir/$group/*/mapping/*coverage_adapted_filtered_*_term.fasta_[0-9][0-9] \
|sort -u | sed 's/^>//g' > $summaryDir/sample_ac_common_length_$group.txt

#NC_016138.1	123322
#NC_019202.1	103532
#NC_020280.1	5022
#NC_022242.1	28466
#NC_022739.1	321653

awk '{split($2,ac_description,">");split(ac_description[2],ac_only," ")} {print ac_only[1],"\t", $0}' $summaryDir/sample_description_list_$group.txt | sort > $summaryDir/sample_ac_common_join_$group.txt

#NC_016138.1 	 PAE1294 	 >NC_016138.1 Pseudomonas aeruginosa plasmid pUM505, complete sequence
#NC_016138.1 	 PAE1305 	 >NC_016138.1 Pseudomonas aeruginosa plasmid pUM505, complete sequence
#NC_016138.1 	 PAE1323 	 >NC_016138.1 Pseudomonas aeruginosa plasmid pUM505, complete sequence
#NC_016138.1 	 PAE1349 	 >NC_016138.1 Pseudomonas aeruginosa plasmid pUM505, complete sequence
#NC_016138.1 	 PAE1359 	 >NC_016138.1 Pseudomonas aeruginosa plasmid pUM505, complete sequence

echo "adding length info to summary file"

join $summaryDir/sample_ac_common_length_$group.txt $summaryDir/sample_ac_common_join_$group.txt > $summaryDir/ac_length_sample_description_$group.txt


#NC_005246.1 60145 A14_03 >NC_005246.1 Erwinia amylovora LebB66 plasmid pEL60, complete sequence
#NC_005246.1 60145 G00_07 >NC_005246.1 Erwinia amylovora LebB66 plasmid pEL60, complete sequence
#NC_005246.1 60145 K3087 >NC_005246.1 Erwinia amylovora LebB66 plasmid pEL60, complete sequence
#NC_005246.1 60145 K3167 >NC_005246.1 Erwinia amylovora LebB66 plasmid pEL60, complete sequence
#NC_009649.1 175879 A14_03 >NC_009649.1 Klebsiella pneumoniae subsp. pneumoniae MGH 78578 plasmid pKPN3, complete sequence

#awk '{print substr($0, index($0,$3))}'

#Adapt to sequence description that doesn't have all fields required (19 june 2018)

length_description=$(awk 'END{print NF}' $summaryDir/ac_length_sample_description_$group.txt)
echo "preparing table with plasmid description $group"
if [ $length_description -ge 7 ]; then

	awk '{gsub(/complete sequence/,"") gsub(",","")} {print $1, "\t", $2, "\t", $3, "\t", $5, $6, "\t",substr($0, index($0,$7))}' \
	$summaryDir/ac_length_sample_description_$group.txt | sort > $summaryDir/table_ac_length_sample_description_$group.tsv

	awk '{gsub(/complete sequence/,"") gsub(",","")} {print $1, "\t", $2, "\t", $5, $6, "\t", substr($0, index($0,$7))}' \
	$summaryDir/ac_length_sample_description_$group.txt | sort -u | cut -d "," -f 1 > $summaryDir/table_ac_length_NO_sample_description_$group.tsv
else
	awk '{gsub(/complete sequence/,"") gsub(",","")} {print $1, "\t", $2, "\t", $3, "\t", $5, $6, "\t",substr($0, index($0,$NF))}' \
	$summaryDir/ac_length_sample_description_$group.txt | sort > $summaryDir/table_ac_length_sample_description_$group.tsv

	awk '{gsub(/complete sequence/,"") gsub(",","")} {print $1, "\t", $2, "\t", $5, $6, "\t", substr($0, index($0,$NF))}' \
	$summaryDir/ac_length_sample_description_$group.txt | sort -u | cut -d "," -f 1 > $summaryDir/table_ac_length_NO_sample_description_$group.tsv
fi

#NZ_CP011990.1    162533  3950_S3         Klebsiella pneumoniae   UHKPC33 plasmid pUHKPC33-162.533kb, NZ_CP011990.1 
#NZ_CP011990.1    162533  K5496_S4        Klebsiella pneumoniae   UHKPC33 plasmid pUHKPC33-162.533kb, NZ_CP011990.1 
#NZ_CP011990.1    162533  K5761_S8        Klebsiella pneumoniae   UHKPC33 plasmid pUHKPC33-162.533kb, NZ_CP011990.1 
#NZ_CP012803.1    92739   2013-16_S18     Escherichia coli        O157:H7 strain WS4202 NZ_CP012803.1 plasmid pO157-WS4202
#NZ_CP014006.1    46161   K5761_S8        Klebsiella pneumoniae   subsp. pneumoniae strain NZ_CP014006.1 NUHL24835 plasmid unnamed2


echo "obtaining table with presence/absence $group"

echo "creating a file with all combinations"

for sampleid in $(cat $plasmidIdDir/$group/00_summary/sampleid_found_$group.txt);do \
  for ac in $(cat $summaryDir/sample_ac_common_list_$group.txt); do \
    echo -e $sampleid "\t" $ac
  done
done > $summaryDir/sample_ac_combination_list_$group.txt

#2013-16_S18    NC_023027.1
#2013-16_S18    NC_023314.1
#2013-16_S18    NC_024992.1
#2013-16_S18    NC_025007.1
#2013-16_S18    NZ_CP006924.1

echo "checking presence/absence"

awk 'NR==FNR{comb[$1" "$2];next}{table[$3" "$1]}END{for (i in comb) {for (j in table) {if (i ~ j) {print i, 1;} else {print i, 0}}}}' \
$summaryDir/sample_ac_combination_list_$group.txt $summaryDir/table_ac_length_sample_description_$group.tsv | sort -u \
> $summaryDir/table_ac_sample_presence_$group.tsv

#2013-16_S18 NZ_LN824135.1 0
#2013-16_S18 NZ_LN868944.1 0
#2013-16_S18 NZ_LN868944.1 1
#2719_S10 NC_004989.1 0
#2719_S10 NC_006855.1 0

echo "changing presence to percentage"

awk 'NR==FNR{comb[$1" "$2];next}{table[$1" "$2]}{percentage[$1" "$2]=$3}END{for (i in comb) {for (j in table) {if (i ~ j) {print i, percentage[i];} else {print i, 0}}}}' \
$summaryDir/sample_ac_combination_list_$group.txt $summaryDir/sample_ac_percentage_$group.txt | sort -u \
> $summaryDir/table_ac_sample_percentage_$group.tsv

#MPVECO111_S45 NC_025179.1 0
#MPVECO111_S45 NC_025179.1 80.9484
#MPVECO111_S45 NC_025180.1 0
#MPVECO111_S45 NZ_AFVX01000096.1 82.4395
#MPVECO111_S45 NZ_AGTD01000002.1 0

#Remove false negatives by merging them with true positives and format output as table


echo "sorting presence/absence table $group"


columnNumber=$(cat $plasmidIdDir/$group/00_summary/sampleid_found_$group.txt | wc -l)
rowNumber=$(cat $summaryDir/sample_ac_common_list_$group.txt | wc -l )

awk '{a[$1" "$2]++};END{for (i in a) print i, a[i]}' $summaryDir/table_ac_sample_presence_$group.tsv | sort > $summaryDir/table_ac_sample_presence_fixed_$group.tsv


awk -v nr=$rowNumber '{gsub("1","0") gsub("2","1")} {a[NR]=$3} END{for (i=1;i<=nr;i++) {for (j=i;j<=NR;j+=nr) printf "%s\t",a[j]; print""}} END{gsub("1","0"); gsub("2","1")}' \
$summaryDir/table_ac_sample_presence_fixed_$group.tsv > $summaryDir/table_presence.plasmids_$group.tsv

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

sort -k 3 -nr $summaryDir/table_ac_sample_percentage_$group.tsv | awk '(!first[$1$2]++)' | sort > $summaryDir/table_ac_sample_percentage_fixed_$group.tsv

awk -v nr=$rowNumber '{a[NR]=$3} END{for (i=1;i<=nr;i++) {for (j=i;j<=NR;j+=nr) printf "%s\t",a[j]; print""}}' \
$summaryDir/table_ac_sample_percentage_fixed_$group.tsv > $summaryDir/table_percentage.plasmids_$group.tsv


echo "merging all data in a final report table $group"

paste <(awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' $summaryDir/table_presence.plasmids_$group.tsv) <(cat $summaryDir/table_percentage.plasmids_$group.tsv)  > $summaryDir/table_report.plasmids_$group.tsv

printf "%s\t" AC_Number Length Species Description N $(cat $plasmidIdDir/$group/00_summary/sampleid_found_$group.txt | sort) > $summaryDir/header_column_$group.tsv

printf "\n" >> $summaryDir/header_column_$group.tsv

paste $summaryDir/table_ac_length_NO_sample_description_$group.tsv $summaryDir/table_report.plasmids_$group.tsv > $summaryDir/FINAL_REPORT_NO_samples_$group.tsv

cat $summaryDir/header_column_$group.tsv $summaryDir/FINAL_REPORT_NO_samples_$group.tsv > $summaryDir/FINAL_REPORT_$group.tsv

echo "All done with $group , Check FINAL_REPORT_$group.tsv"

#Clean all tmp files

for i in $(ls $summaryDir | grep -v FINAL_REPORT_$group.tsv)
do
	rm $summaryDir/$i
done