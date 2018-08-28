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
#CREATED: 06 June 2018
#REVISION:
#DESCRIPTION:rename_from_fasta script rename any field in a file by either providing two fasta files or a dictionary file

#================================================================
# END_OF_HEADER
#================================================================


usage() {
	cat << EOF

rename_from_fasta script rename any field in a file by either providing two fasta files or a dictionary file

usage : $0 <-i file_to_rename> [-1 <inputfile1(.fasta)>] [-2 <inputfile2(.fasta)>] [-d <dictionary>] [-o <directory>] [-f <file_name>] [-v] [-h]

	-i input file to rename
	-1 original fata file whose names will be finally printed
	-2 new fata file whose names will be replaced
	-o output directory (optional). By default the file is replaced in the same location
	-f output file name (".rename" will be added at the end)
	-d dictionary file to be used if fasta files are not supplied
	-v version
	-h display usage message

example: process_cluster_output.sh -i ecoli_clustered.fasta_70 -b ecoli.coverage

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

# Error handling
error(){
  local parent_lineno="$1"
  local script="$2"
  local message="$3"
  local code="${4:-1}"

	RED='\033[0;31m'
	NC='\033[0m'

  if [[ -n "$message" ]] ; then
    echo -e "\n---------------------------------------\n"
    echo -e "${RED}ERROR${NC} in Script $script on or near line ${parent_lineno}; exiting with status ${code}"
    echo -e "MESSAGE:\n"
    echo -e "$message"
    echo -e "\n---------------------------------------\n"
  else
    echo -e "\n---------------------------------------\n"
    echo -e "${RED}ERROR${NC} in Script $script on or near line ${parent_lineno}; exiting with status ${code}"
    echo -e "\n---------------------------------------\n"
  fi

  exit "${code}"
}

#DECLARE FLAGS AND VARIABLES
cwd="$(pwd)"
input_file="Input_file"

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:1:2:f:o:d:vh"
while getopts $options opt; do
	case $opt in
		i )
			input_file=$OPTARG
			;;
		1 )
			fasta_file_old=$OPTARG
			;;
		2 )
			fasta_file_new=$OPTARG
			;;
		d )
			dictionary_file_new=$OPTARG
			;;
		o )
			output_dir=$OPTARG
			;;
		f )
			file_name=$OPTARG
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

check_mandatory_files.sh $input_file

if [ ! $output_dir ]; then
	output_dir=$(dirname $input_file)
	echo "Default output directory is" $output_dir
	mkdir -p $output_dir
else
	echo "Output directory is" $output_dir
	mkdir -p $output_dir
fi


if [ ! $file_name ]; then
	file_name=$(basename $input_file | cut -d "." -f1,2)
fi

fasta_file_old_name=$(basename $fasta_file_old)
fasta_file_new_name=$(basename $fasta_file_new)

echo "$(date)"
echo "Renaming" $file_name

cat $fasta_file_old | awk '/>/ {print $1}'| sed 's/>//g' | sed 's/|/-/g' > $output_dir/$fasta_file_old_name".ac"
cat $fasta_file_new | awk '/>/ {print $1}'| sed 's/>//g' | sed 's/|/-/g' > $output_dir/$fasta_file_new_name".ac"
cat $input_file | sed 's/|/-/g' > $output_dir/$file_name".nopipe.tmp"


#Paste colums to relate names in a dictionary
awk 'NR==FNR{ac[NR]=$0;next}{print ac[FNR], "\t", $0 }' $output_dir/$fasta_file_old_name".ac" $output_dir/$fasta_file_new_name".ac" > $output_dir/dictionary.txt || error ${LINENO} $(basename $0) "AWK command failed in dictionary.txt creation. See $output_dir/logs for more information."

#Rename fields

#cat $output_dir/dictionary.txt | while read -r line; do word1=$(cut -f1); word2=$(cut -f2); echo "##########word 1="$word1;echo "###########word 2="$word2; sed 's/$word2/$word1/g' $input_file; done > $output_dir/$file_name".renamed"


awk 'FNR==NR {dict[$2]=$1; next} {for (i in dict) gsub(i, dict[i])}1' $output_dir/dictionary.txt $output_dir/$file_name".nopipe.tmp" > $output_dir/$file_name".renamed" || error ${LINENO} $(basename $0) "AWK command failed in $file_name\".renamed\" creation. See $output_dir/logs for more information."

#awk 'NR==FNR{dict[$2]=$1;next}{$1=dict[$1]}1' $output_dir/dictionary.txt $input_file #> $output_dir/$file_name".renamed"


rm $output_dir/$fasta_file_old_name".ac"
rm $output_dir/$fasta_file_new_name".ac"
rm $output_dir/$file_name".nopipe.tmp"
rm $output_dir/dictionary.txt

echo "$(date)"
echo "DONE renaming" $file_name
echo -e "Renamed file can be found at" $output_dir/$file_name".renamed"
