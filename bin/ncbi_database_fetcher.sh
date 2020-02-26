#!/bin/bash

# Exit immediately if a pipeline, which may consist of a single simple command, a list, 
#or a compound command returns a non-zero status: If errors are not handled by user
set -e
#set -x
#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.0 
#CREATED: 12 June 2018
#REVISION:
# 22 June 2018: include quite mode that avoid watching the progress
#		
#
#DESCRIPTION:Script that extract a database from ncbi database using terms
#AKNOWLEDGE: 
#		-Multiple arguments in one flag: https://stackoverflow.com/questions/7529856/retrieving-multiple-arguments-for-a-single-option-using-getopts-in-bash
#
#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

ncbi_database_fetcher is a script that extract sequences from NCBI by term

usage : $0 <(-y term1 -y term2 | -y "term1 term2")> [(-n term1 -n term2 | -n "term1 term2")] [-O <organism>][-d (nucleotide|protein)] [-f <filename>] [-o <directory>] [-q] [-v] [-h]

	-y list of key terms separated by space to be INCLUDED in sequences title
	-n list of key terms separated by space to be EXCLUDED in sequences title
	-O organism to filter
	-d database type, default nucleotide
	-o output directory (optional). By default the file is placed in cwd
	-f file name (optional). By default is the first term used as query
	-q quiet
	-v version
	-h display usage message

example: ./ncbi_database_fetcher.sh -y plasmid -n unnamed -n partial -O Archaea

EOF
}

#================================================================
# OPTION_PROCESSING
#================================================================
#Make sure the script is executed with arguments
if [ $# = 0 ] ; then
 usage >&2
 exit 1
fi

#DECLARE FLAGS AND VARIABLES
cwd="$(pwd)"
use_term_and=false
use_term_not=false
use_term_org=false
quiet=false
database_type=nucleotide
#PARSE VARIABLE ARGUMENTS WITH getops

options=":y:n:o:f:d:O:qvh"
while getopts $options opt; do
	case $opt in
		o )
			output_dir=$OPTARG
			;;
		O)
			terms_organism+=($OPTARG)
			use_term_org=true
			;;
		f )
			file_name=$OPTARG
			;;
		d )
			database_type=$OPTARG
			;;
		y )
			terms_and+=($OPTARG)
			use_term_and=true
			;;
		n )
			terms_not+=($OPTARG)
			use_term_not=true
			;;
		q )
			quiet=true
			;;
		h )
		  	usage
		  	exit 1
		  	;;
		v )
		  	echo $VERSION
		  	exit 1
		  	;;
		\?)  
			echo "Invalid Option: -$OPTARG" 1>&2
			usage
			exit 1
			;;
		: )
      		echo "Option -$OPTARG requires an argument." >&2
      		exit 1
      		;;
      	* ) 
			echo "Unimplemented option: -$OPTARG" >&2;
			exit 1
			;;

	esac
done
shift $((OPTIND-1))

#================================================================
# MAIN_BODY
#================================================================
##CHECK DEPENDENCIES, MANDATORY FIELDS, FOLDERS AND ARGUMENTS

echo -e "\n#Executing" $0 "\n"

if [ $use_term_and = false ]; then
	echo "Please, introduce at least one term to include search"
	usage
	exit 1
fi

#MANAGE OUTPUT DIRECTORY
if [ ! $output_dir ]; then
	output_dir=$cwd
	echo "Default output_dir is" $output_dir
	mkdir -p $output_dir
else
	echo "Output directory is" $output_dir
	mkdir -p $output_dir
fi

#MANAGE FILE NAME

if [ ! $file_name ]; then

	if [ "${#terms_and[@]}" -gt 1 ]; then
		file_name_value_one=$(echo ${terms_and[0]})
		file_name_value_two=$(echo ${terms_and[1]})

		file_name=$file_name_value_one"_"$file_name_value_two
		echo "Default file name is" $file_name
	else
		file_name=$terms_and".database"
		echo "Default file name is" $file_name
	fi
else
	echo "File name is" $file_name
fi


##PROCESS REGULAR EXPRESSION TERMS

list_terms_and=$(for term in "${terms_and[@]}"; do echo "$term"; done)
list_terms_org=$(for organism in "${terms_organism[@]}"; do echo "$organism"; done)

#echo "${#terms_and[@]}" "NUMBER OF TERMS"

list_terms_regexp_and=$(printf "%s[Title] AND " $list_terms_and | sed 's/ AND $//g')
list_terms_regexp_organism=$(printf "AND %s[organism] " $list_terms_org | sed 's/ $//g')

if [ $use_term_not = true ]; then

	list_terms_not=$(for term in "${terms_not[@]}"; do echo "$term"; done)
	list_terms_regexp_not=$(printf "NOT %s[Title] " $list_terms_not | sed 's/ $//g')
	final_list_terms_regexp=$(echo $list_terms_regexp_and" "$list_terms_regexp_not" "$list_terms_regexp_organism) #concat all regexp into one

else
	final_list_terms_regexp=$(echo $list_terms_regexp_and " "$list_terms_regexp_organism)
fi

echo $final_list_terms_regexp

########EUTILS COMMAND############
##################################

echo "$(date)"
echo "Obtaining seuences with terms:" $list_terms_and
echo "But not those terms:" $list_terms_not
if [ $use_term_org = true ]; then
	echo "Filtering by organisms:" $list_terms_org
fi
echo ""

base="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

##DETERMINE RETMAX
wget -q -O $output_dir/$file_name".count" $base"esearch.fcgi?db="$database_type"&term=""$final_list_terms_regexp"

counter=$(cat $output_dir/$file_name".count" | awk '/<Count>/' | head -n 1 | awk '/<Count>/ {split($0,counter_prev,"</Count>");split(counter_prev[1],counter,"<Count>")}END{print counter[length(counter)]}')
echo -e "FOUND" $counter "RECORDS\n"

if [ $counter -eq 0 ]; then
	echo "Try different terms"
	echo "EXIT"
	exit 1
fi

echo "Retrieving Id"

##OBTAIN TOTAL LIST OF ID
wget -q -O $output_dir/$file_name".id" $base"esearch.fcgi?db="$database_type"&term=""$final_list_terms_regexp""&RetMax="$counter

list_of_id=$(cat $output_dir/$file_name".id"| awk '{split($0,id_prev,"</Id>");split(id_prev[1],id,"<Id>")}/<Id>/{print id[length(id)]}')
array_of_id=($list_of_id)

echo "And sequences"
counter=1


##Checking previous DDBB
if [ -s $output_dir/$file_name".fasta" ]; then
	echo -e "\nFound a ddbb with the same name, Removing it\n"
	rm $output_dir/$file_name".fasta"
fi


##RETRIEVING FASTA SEQUENCE 

for i in $list_of_id
do 
	if [ $quiet = false ]; then

		echo $counter"/""${#array_of_id[@]}"
	fi

	((counter++))
	
	curl -s $base"efetch.fcgi?db="$database_type"&id="$i"&retmode=text&rettype=fasta" >> $output_dir/$file_name".fasta"
done


echo "$(date)"
echo "DONE obtaining seuences with terms supplied"

seq_number_post=$(cat $output_dir/$file_name".fasta" | grep ">" | wc -l)
echo "File with filtered sequences can be found in" $output_dir/$file_name".fasta"
echo "with" $seq_number_post "sequences"

rm $output_dir/$file_name".count"
rm $output_dir/$file_name".id"
