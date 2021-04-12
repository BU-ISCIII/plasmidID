#!/bin/bash

#set -e

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.0
#CREATED: 30 April 2018
#REVISION:
#12 June 2018: Handled cleaning process without hard coded paths
#
#DESCRIPTION:Script that uses prokka to annotate a FASTA file
#
#DOCUMENTATION
#
#Prokka outputs the fasta headers as:
# gnl|center|locustag_01
# gnl|center|locustag_02
#
#TO DO:
#Handle cleaning [v]
#
#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

Prokka_annotation is a script that uses prokka to annotate a FASTA file

usage : $0 <-i inputfile(FASTA)> <-p prefix> [-o <directory>] [-k <kingdom>]
		[-T <threads>] [-g group_name][-G genus] [-S species] [-c] [-v] [-h]

	-i input file in FASTA format
	-o output directory
	-p prefix for sample identification (mandatory) and output file name
	-k kingdom (Bacteria by default)
	-g group name (optional). If unset, samples will be gathered in NO_GROUP group
	-G sample genus in case is known by user
	-S sample species in case is known by user
	-c clean:remove files other than gff and renamed fasta
	-T number of threads
	-v version
	-h display usage message


Output directory is the same as input directory by default

example: prokka_annotation -i ecoli.fasta -p ECO -T 5


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
group="NO_GROUP"
input_file="Input_file"
kingdom="Bacteria"
clean=false
genus=""
species=""
threads=1

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:o:p:k:g:G:S:T:cvh"
while getopts $options opt; do
	case $opt in
		i )
			input_file=$OPTARG
			;;

		o )
			output_dir=$OPTARG
			;;
		p)
			prefix=$OPTARG
			file_name=$OPTARG
			;;
		k )
          	kingdom=$OPTARG
          	;;
        g )
          	group=$OPTARG
          	;;
        S )
          	species=$OPTARG
          	;;
        G)
          	genus=$OPTARG
          	;;
		c )
          	clean=true
          	;;
        T)
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

check_dependencies.sh prokka
echo "PREFIX" $prefix

if [ ! $prefix ]; then
	echo "please provide a prefix"
	exit 1
fi

if [ ! $output_dir ]; then
	output_dir=$(dirname $input_file)
	echo "Default output directory is" $output_dir
	mkdir -p $output_dir
else
	echo "Output directory is" $output_dir
	mkdir -p $output_dir
fi

if [ ! $file_name ]; then
	file_name=$(basename $input_file)
	echo "filename is" $file_name
fi


##PROKKA EXECUTION

echo "$(date)"
echo "Annotating $input_file with prokka"

prokka --force --outdir $output_dir \
--prefix $prefix \
--addgenes \
--kingdom $kingdom \
--genus $genus \
--species $species \
--usegenus \
--centre BU-ISCIII \
--locustag $prefix \
--addgenes \
--cpus $threads \
$input_file #|| error ${LINENO} $(basename $0) "Prokka command failed. See $output_dir/logs for more information."

echo "$(date)"
echo "done annotating $input_file with prokka"

##CLEAN FILES THAT WILL NOT BE USED IN PLASMIDID

if [ $clean = true ]; then

	echo "Removing unwanted files"
	for i in $(ls $output_dir/$prefix.??? | awk '!/fna|gff|log|err|gb/')
	do
		rm $i
	done
fi

echo -e "\n"
