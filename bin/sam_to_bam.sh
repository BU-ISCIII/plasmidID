#!/bin/bash

# Exit immediately if a pipeline, which may consist of a single simple command, a list,
#or a compound command returns a non-zero status: If errors are not handled by user
#set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion.
# An error message will be written to the standard error, and a non-interactive shell will exit
#set -u
#Print everything as if it were executed, after substitution and expansion is applied: Debug|log option
#set -x

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.0
#CREATED: 19 March 2018
#REVISION:
#DESCRIPTION:Script that convert a supplied SAM file into compressed binary indexed BAM

#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

Sam_to_bam script converts a supplied SAM file into compressed binary indexed BAM

usage : $0 <-i inputfile(.sam)> [-o <directory>] [-s sample_name] [-g group_name] [-T <int>] [-v] [-h]

	-i input file
	-o output directory (optional). By default the BAM file will replace SAM in the same location
	-s sample name
	-g group name (optional). If unset, samples will be gathered in NO_GROUP group
	-T number of threads
	-v version
	-h display usage message

example: sam_to_bam.sh -i ecoli.sam

EOF
}

#================================================================
# OPTION_PROCESSING
#================================================================
#Make sure the script is executed with arguments
if [ $? != 0 ] ; then
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
threads=1
cwd="$(pwd)"
group="NO_GROUP"
input_file="Input_file"

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:o:s:g:vh"
while getopts $options opt; do
	case $opt in
		i )
			input_file=$OPTARG
			;;
		o )
			output_dir=$OPTARG
			;;
		s )
			sample=$OPTARG
			;;
		g)
			group=$OPTARG
			;;

        T )
			threads=$OPTARG
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

check_dependencies.sh samtools


if [ ! $output_dir ]; then
	output_dir=$(dirname $input_file)
	echo "Default output directory is" $output_dir
	mkdir -p $output_dir
else
	echo "Output directory is" $output_dir
	mkdir -p $output_dir
fi

if [ ! $sample ]; then
	sample=$(basename $input_file | cut -d. -f1)
fi

########SAM_TO_BAM##########
############################


if [ -f $output_dir/$sample.sorted.bam -a -f $output_dir/$sample.sorted.bam.bai  ];then \
	echo "Found a sorted .BAM file for sample" $sample;
	echo "Omitting BAM to SAM convertion"
else
	echo "$(date)"
	echo "Converting SAM to sorted indexed BAM in $sample"

	samtools view \
	-Sb $input_file \
	-o $output_dir/$sample.bam || error ${LINENO} $(basename $0) "Samtools view command failed. See $output_dir/logs for more information."


	echo "$(date)"
	echo "Sorting BAM file in $sample"

	samtools sort \
	-T $output_dir/$sample".sorted.bam" \
	-o $output_dir/$sample".sorted.bam" \
	$output_dir/$sample.bam || error ${LINENO} $(basename $0) "Samtools sort command failed. See $output_dir/logs for more information."

	echo "$(date)"
	echo "Indexing BAM file in $sample"

	samtools index \
	$output_dir/$sample".sorted.bam" || error ${LINENO} $(basename $0) "Samtools index command failed. See $output_dir/logs for more information."


	echo "$(date)"
	echo "DONE Converting SAM to sorted indexed BAM in $sample"
fi

if [ -f $output_dir/$sample.sam ];then \

	echo $sample.sam "removed"
	rm $output_dir/$sample.sam

fi

if [ -f $output_dir/$sample.bam ];then \

	echo $sample.bam "removed"
	rm $output_dir/$sample.bam

fi

echo -e "\n"
