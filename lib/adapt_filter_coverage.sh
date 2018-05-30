#!/bin/bash

# Exit immediately if a pipeline, which may consist of a single simple command, a list, 
#or a compound command returns a non-zero status: If errors are not handled by user
set -e
#~set -x

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.0 
#CREATED: 21 March 2018
#REVISION:
#DESCRIPTION:adapt_filter_coverage script that adapt percentages and filter coverage info from bedtools genomecov output

#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

adapt_filter_coverage script that adapt percentages and filter coverage info from bedtools genomecov output

usage : $0 <-i inputfile(.fasta)> [-o <directory>] [-c <int(0-100)>] [-s <suffix>] [-v] [-h]

	-i input file 
	-o output directory (optional). By default the file is replaced in the same location
	-c percentage value to filter >= values. If not supplied, all records will be outputted
	-s string to ad at the end of the outputted file (list of accession numbers)
	-v version
	-h display usage message

example: adapt_filter_coverage.sh -i ecoli.coverage -c 70

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
input_file="Input_file"
coverage_cutoff_input=100

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:o:c:s:vh"
while getopts $options opt; do
	case $opt in
		i )
			input_file=$OPTARG
			;;
		o )
			output_dir=$OPTARG
			;;
		c )
			if [ $OPTARG -lt 0 ] || [ $OPTARG -gt 100 ]; then
				echo "please, provide a percentage between 0 and 100"
				usage
				exit 1
			else
				coverage_cutoff_input=$OPTARG
			fi
			;;
		s )
			suffix=$OPTARG
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


suffix="_adapted_filtered_"$coverage_cutoff_input
coverage_cutoff=$(echo "(1 - ($coverage_cutoff_input/100))" | bc -l)

#echo $coverage_cutoff

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


if [ -f $input_file"_adapted" ]; then
	echo "Found previous" $(basename $input_file"_adapted")", removing it"
	rm $input_file"_adapted"
fi

awk '
BEGIN{OFS="\t"}
(!x[$1]++) {if ($1 != "genome") 
				{if ($2 == 0) 
					{print $0} 
				else 
					{print $1, 0, $4, $4, 0.0000000001}
				}
			}
	' $input_file > $input_file"_adapted"

awk '
{if ($2 == 0 && $5 < '"${coverage_cutoff}"')
	 {print $1}
}
	' $input_file"_adapted" > $input_file$suffix


echo "$(date)"
echo "Done filtering sequences with" $coverage_cutoff_input"% and greater coverage"
echo "Those sequences can be found at" $input_file$suffix
echo -e $(cat $input_file$suffix | wc -l) mapped equals or more than $coverage_cutoff_input "\n"