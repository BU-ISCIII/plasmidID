#!/bin/bash

# Exit immediately if a pipeline, which may consist of a single simple command, a list, 
#or a compound command returns a non-zero status: If errors are not handled by user
set -e

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.0 
#CREATED: 20 March 2018
#REVISION:
#DESCRIPTION:Script that convert a supplied SAM file into compressed binary indexed BAM
#AKNOWLEDGE: 
#		-Adapted from klashxx: https://stackoverflow.com/questions/23992646/sequence-length-of-fasta-file/23992773
#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

Calculate_sequlen script calculates a supplied FASTA length

usage : $0 <-i inputfile(.fasta)> [-o <directory>] [-n <string>] [-r] [-v] [-h]

	-i input file 
	-o output directory (optional). By default the file is replaced in the same location
	-n file name (optional). By default is the same name with .length extension
	-r remove ">" (greater-than) symbol from fasta header
	-v version
	-h display usage message

example: calculate_sequlen.sh -i ecoli.fasta 

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
remove_head=remove_head_false
cwd="$(pwd)"
file_name="file_name"
input_file="Input_file"

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:o:n:rvh"
while getopts $options opt; do
	case $opt in
		i )
			input_file=$OPTARG
			;;
		o )
			output_dir=$OPTARG
			;;
		n )
			file_name=$OPTARG
			;;
		r )
			remove_head="^>"
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

bash lib/check_mandatory_files.sh $input_file

if [ ! $output_dir ]; then
	output_dir=$(dirname $input_file)
	#echo "Default output directory is" $output_dir
	mkdir -p $output_dir
else
	#echo "Output directory is" $output_dir
	mkdir -p $output_dir
fi

if [ ! $filename ]; then
	filename=$(basename $input_file | cut -d. -f1)
fi

awk '
BEGIN {FS=="| "}
/^>/ {if (seqlen)
	print seqlen;printf "%s\t", $1; seqlen=0; next
	} 
{seqlen+=length($0)}
END {print seqlen}' $input_file | sed 's/'$remove_head'//g' \
>$output_dir/$filename".length"

