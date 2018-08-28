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
#CREATED: 12 April 2018
#REVISION:
#DESCRIPTION:process_cluster_output script obtain a list of ac from fasta, and estract their coverage value from a coverage file

#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

process_cluster_output script obtain a list of ac from fasta, and estract their coverage value from a coverage file

usage : $0 <-i inputfile(.fasta)> <-b coverage_file> [-o <directory>] [-c <int(0-100)>] [-s <suffix>] [-v] [-h]

	-i input file
	-b file with coverage info
	-o output directory (optional). By default the file is replaced in the same location
	-c percentage value to filter >= values. If not supplied, all records will be outputted
	-s string to ad at the end of the outputted file (list of accession numbers)
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
coverage_cutoff_input=100

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:b:o:c:s:vh"
while getopts $options opt; do
	case $opt in
		i )
			input_file=$OPTARG
			;;
		b )
			coverage_file=$OPTARG
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

check_mandatory_files.sh $input_file

suffix="_clustered"
coverage_cutoff=$(echo "(1 - ($coverage_cutoff_input/100))" | bc -l)

if [ ! $output_dir ]; then
	output_dir=$(dirname $input_file)
	#echo "Default output directory is" $output_dir
	mkdir -p $output_dir
else
	#echo "Output directory is" $output_dir
	mkdir -p $output_dir
fi


if [ ! $file_name ]; then
	file_name=$(basename $input_file)
	coverage_name=$(basename $coverage_file)
fi

echo "$(date)"
echo "extracting coverage info from clustered sequences in" $file_name

ac_input_file=$(cat $input_file | grep ">" | awk '{gsub(">","");print $1}')

for i in $ac_input_file ;do
	awk '
		/^'"$i"'/
		' $coverage_file
done > $output_dir/$coverage_name$suffix || error ${LINENO} $(basename $0) "Awk command error in $coverage_name$suffix creation. See $output_dir/logs for more information."


awk '
	{if ($2 == 0 && $5 <= '"${coverage_cutoff}"')
		{print $1}}
	' $output_dir/$coverage_name$suffix > $output_dir/$coverage_name$suffix"_ac" || error ${LINENO} $(basename $0) "Awk command error in $coverage_name$suffix\"_ac\" creation. See $output_dir/logs for more information."


awk '
	{if ($2 == 0 && $5 <= '"${coverage_cutoff}"')
	 	{print $1, ((1 - $5)*100)}
	}
	' $output_dir/$coverage_name$suffix > $output_dir/$coverage_name$suffix"_percentage" || error ${LINENO} $(basename $0) "Awk command error in $coverage_name$suffix\"_percentage\" creation. See $output_dir/logs for more information."

echo "$(date)"
echo "DONE extracting coverage info from clustered sequences in" $file_name
echo -e "Info can be found at" $coverage_name$suffix"_ac and" "\n" $coverage_name$suffix"_percentage" "\n"
