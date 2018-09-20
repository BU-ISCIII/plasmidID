#!/bin/bash

set -e
#set -x

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.0
#CREATED: 21 May 2018
#REVISION:
#DESCRIPTION:Script that assemble illumina sequences using SPAdes
#
#
#================================================================
# END_OF_HEADER
#================================================================


usage() {
	cat << EOF

spades_assembly script that assemble illumina sequences using SPAdes

usage : $0 <-p R1_paired file> <-P R2_paired file> [-u <R1_unpaired>] [-U <R2_unpaired>] [-o <directory>]
		 [-k <int>][-s sample_name] [-g group_name] [-f <file_name>] [-T <int>] [q] [-c] [-v] [-h]

	-p R1_paired file (mandatory)
	-P R2_paired file (mandatory)
	-u R1_unpaired file
	-U R2_unpaired file
	-k kmers, supplied as numbers sepparated by number or one flag per number, default: 21,33,55,77,99,127
	-o output directory (optional)
	-f file name
	-s sample name (mandatory)
	-g group name (optional). If unset, samples will be gathered in NO_GROUP group
	-q quick_mode: look for files in a folder SUPPLIED with "paired" and "unpaired" term
	-c clean mode: remove unnecesary temporary folders
	-T threads, default 1
	-v version
	-h display usage message

example: ./spades_assembly.sh -p ecoli_R1_paired.fastq.gz -P ecoli_R2_paired.fastq.gz -c

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
r1_paired_file="R1_paired_file"
r2_paired_file="R2_paired_file"
r1_unpaired_command=""
r2_unpaired_command=""
threads=1
kmer_values_command="21,33,55,77,99,127"
kmer_option=false
quick_mode=false
clean_mode=false

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":p:P:u:U:o:f:d:a:s:g:k:T:q:cvh"
while getopts $options opt; do
	case $opt in
		p )
			r1_paired_file=$OPTARG
			;;
		P )
			r2_paired_file=$OPTARG
			;;
		u )
			r1_unpaired_file=$OPTARG
			r1_unpaired_command=$(echo "--s1" $r1_unpaired_file)
			;;
		U )
			r2_unpaired_file=$OPTARG
			r2_unpaired_command=$(echo "--s1" $r2_unpaired_file)
			;;
		o )
			output_dir=$OPTARG
			;;
		f )
			file_name=$OPTARG
			;;
		s )
			sample=$OPTARG
			;;
		k)
			kmer_value+=($OPTARG)
			kmer_option=true
			;;
		q)
			directory_reads=$OPTARG
			quick_mode=true
			;;
		l)
			minimus_length=$OPTARG
			;;
		g)
			group=$OPTARG
			;;
		c)
			clean_mode=true
			;;
		M )
			max_mem=$OPTARG
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

check_dependencies.sh spades.py


if [ ! $directory_reads ]; then
	directory_reads=$(dirname $r1_paired_file)
	echo "Reads directory is" $directory_reads
else
	echo "Reads directory for quick mode is" $directory_reads
	sample_dir=$(dirname $directory_reads)
	output_dir=$sample_dir"/assembly"
	mkdir -p $output_dir
fi


if [ ! $output_dir ]; then
	sample_dir=$(dirname $directory_reads)
	output_dir=$sample_dir"/assembly"
	echo "Default output directory is" $output_dir
	mkdir -p $output_dir
else
	echo "Output directory is" $output_dir
	mkdir -p $output_dir
fi


if [ $quick_mode = true ]; then
	echo "Entering QUICK MODE"
	r1_paired_file=$(find $directory_reads -name  "*1_paired.fastq.gz" -type f)
	r2_paired_file=$(find $directory_reads -name  "*2_paired.fastq.gz" -type f)
	r1_unpaired_file=$(find $directory_reads -name  "*1_unpaired.fastq.gz" -type f)
	r1_unpaired_command=$(echo "--s1" $r1_unpaired_file)
	r2_unpaired_file=$(find $directory_reads -name  "*2_unpaired.fastq.gz" -type f)
	r2_unpaired_command=$(echo "--s2" $r2_unpaired_file)
fi


check_mandatory_files.sh $r1_paired_file $r2_paired_file

if [ $kmer_option = true ]; then
	list_kmer_values=$(for value in "${kmer_value[@]}"; do echo "$value"; done)
	kmer_values_command=$(printf "%s," $list_kmer_values | sed 's/,$//g')
fi


echo "$(date)"
echo "Assembly:"
echo "R1 paired file = " $r1_paired_file
echo "R2 paired file = " $r2_paired_file
echo "R1 unpaired file = " $r1_unpaired_file
echo "R2 unpaired file = " $r2_unpaired_file


spades.py \
--careful \
-t $threads \
-k $kmer_values_command \
--pe1-1 $r1_paired_file \
--pe1-2 $r2_paired_file \
$r1_unpaired_command \
$r2_unpaired_command \
-o $output_dir || error ${LINENO} $(basename $0) "Spades command failed. See $output_dir/logs for more information."



echo "$(date)"
echo "DONE. Assembled contigs can be found at $output_dir/contigs.fasta:"
echo "DONE. Assembled scaffolds can be found at $output_dir/scaffolds.fasta:"

if [ $clean_mode = true ]; then
	echo "Removing unnecesary folders"
	rm -rf $(find $output_dir -maxdepth 1 -mindepth 1 -type d)
	echo "DONE removing unwanted folders"
fi

echo -e "\n"

























