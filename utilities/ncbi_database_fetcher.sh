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
#		
#
#DESCRIPTION:Script that extract a database from ncbi database using terms
#AKNOWLEDGE: 
#		-Multiple arguments in one flag: https://stackoverflow.com/questions/7529856/retrieving-multiple-arguments-for-a-single-option-using-getopts-in-bash
#TODO:
#		-Add and remove sequences in the same execution
#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

ncbi_database_fetcher is a script that extract sequences by term, either by key or file with a list

usage : $0 <(-y term1 -y term2 | -y "term1 term2")> [(-n term1 -n term2 | -n "term1 term2")] [-d (nucleotide|protein)] [-f <filename>] [-o <directory>]  [-v] [-h]

	-y list of key terms separated by space to be INCLUDED in sequences title
	-n list of key terms separated by space to be EXCLUDED in sequences title
	-d database type, default nucleotide
	-o output directory (optional). By default the file is placed in cwd
	-f file name (optional). By default is the first term used as query
	-v version
	-h display usage message

example: ./ncbi_database_fetcher.sh -y plasmid -y Escherichia -n unnamed -n partial -f 

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
database_type=nucleotide
#PARSE VARIABLE ARGUMENTS WITH getops

options=":y:n:o:f:d:vh"
while getopts $options opt; do
	case $opt in
		o )
			output_dir=$OPTARG
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

	file_name=$terms_and".database"
	echo "Default file name is" $file_name
else
	echo "File name is" $file_name
fi


#PROCESS REGULAR EXPRESSION TERMS

list_terms_and=$(for term in "${terms_and[@]}"; do echo "$term"; done)

#echo "${#terms_and[@]}" "NUMBER OF TERMS"

list_terms_regexp_and=$(printf "%s[Title] AND " $list_terms_and | sed 's/ AND $//g')

if [ $use_term_not = true ]; then

	list_terms_not=$(for term in "${terms_not[@]}"; do echo "$term"; done)
	list_terms_regexp_not=$(printf "NOT %s[Title] " $list_terms_not | sed 's/ $//g')
	final_list_terms_regexp=$(echo $list_terms_regexp_and" "$list_terms_regexp_not) #concat all regexp into one

else
	final_list_terms_regexp=$(echo $list_terms_regexp_and)
fi

echo $final_list_terms_regexp

########EUTILS COMMAND############
##################################

echo "$(date)"
echo "Obtaining seuences with terms:" $list_terms_and
echo "But not those terms:" $list_terms_not
echo ""

base="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


wget -q -O $output_dir/$file_name".count" $base"esearch.fcgi?db="$database_type"&term=""$final_list_terms_regexp"

counter=$(cat $output_dir/$file_name".count" | awk '/<Count>/' | head -n 1 | awk '/<Count>/ {split($0,counter_prev,"</Count>");split(counter_prev[1],counter,"<Count>")}END{print counter[length(counter)]}')
echo -e "FOUND" $counter "RECORDS\n"
echo "Retrieving Id"

wget -q -O $output_dir/$file_name".id" $base"esearch.fcgi?db="$database_type"&term=""$final_list_terms_regexp""&RetMax="$counter

list_of_id=$(cat $output_dir/$file_name".id"| awk '{split($0,id_prev,"</Id>");split(id_prev[1],id,"<Id>")}/<Id>/{print id[length(id)]}')

echo "And sequences"

for i in $list_of_id
do 
	curl -s $base"efetch.fcgi?db="$database_type"&id="$i"&retmode=text&rettype=fasta"
done > $output_dir/$file_name".fasta"


echo "$(date)"
echo "DONE obtaining seuences with terms supplied"

seq_number_post=$(cat $output_dir/$file_name".fasta" | grep ">" | wc -l)
echo "File with filtered sequences can be found in" $output_dir/$file_name".fasta"
echo "with" $seq_number_post "sequences"

rm $output_dir/$file_name".count"
rm $output_dir/$file_name".id"