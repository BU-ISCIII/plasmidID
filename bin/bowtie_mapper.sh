#!/bin/bash

set -e

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.0 
#CREATED: 15 March 2018
#REVISION:
#		19 March 2018: Complete usage info
#		19 March 2018: Check mandatory files. folders and variables
#DESCRIPTION:Script that index a database and map a supplied pair-end sequences
#TODO
#	-Handle files extensions for bowtie, now is fastq by default
#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

Bowtie_mapper script index a database and map a supplied pair-end sequences

usage : $0 [-i <inputfile>] [-o <directory>] <-d database(fasta)> <-s sample_name> <-1 R1> <-2 R2> 
		[-g group_name] [-f <int>] [-T <int>] [-a] [-v] [-h]

	-i input directory (optional)
	-o output directory (optional)
	-d database to map (.fasta)
	-s sample name
	-g group name (optional). If unset, samples will be gathered in NO_GROUP group
	-1 reads corresponding to paired-end R1
	-2 reads corresponding to paired-end R2
	-f offrate index for bowtie_build (optional). Default value 1. for quicker indexing use higher number
	-a use -a mapping (off by default)
	-T number of threads
	-v version
	-h display usage message

example: bowtie_mapper.sh -d database.fasta -s COLI -1 ecoli_1.fastq -2 ecoli_2.fastq -a

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
threads=1
offrate=1
cwd="$(pwd)"
a_mapping=""
group="NO_GROUP"
database="Database"
R1="R1"
R2="R2"

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:o:s:g:d:1:2:f:avh"
while getopts $options opt; do
	case $opt in
		i )
			input_dir=$OPTARG
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
		d )
			database=$OPTARG
			;;
		1 )
			R1=$OPTARG
			;;
		2 )
			R2=$OPTARG
			;;
		f )			
          	offrate=$OPTARG
      		;;
        T ) 
			threads=$OPTARG
            ;;
        a)
			a_mapping="-a"
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

check_dependencies.sh bowtie2-build bowtie2

check_mandatory_files.sh $database $R1 $R2

if [ ! $sample ]; then
	echo "ERROR: please, provide a sample name"
	usage
	exit 1
fi

if [ ! $output_dir ]; then
	output_dir=$cwd"/$group/$sample/mapping/"
	echo "Default output directory is" $output_dir
	mkdir -p $output_dir
else
	echo "Output directory is" $output_dir
	mkdir -p $output_dir
fi


########INDEXING############
############################

files_bt2=$(ls $database*bt2 2> /dev/null | wc -l)


if [ "$files_bt2" = "6" ];then \
	echo "Found an indexed ddbb for" $(basename $database);
	echo "Omitting indexing"
else
	echo "Building index of " $(basename $database);
	bowtie2-build \
	--offrate $offrate \
	$database $database
fi

########MAPPING#############
############################

if [ -f $mappedDir/$sample.sorted.bam -a -f $mappedDir/$sample.sorted.bam.bai -o -f $mappedDir/$sample.sam ];then \
	echo "Found a mapping file for sample" $sample;
	echo "Omitting mapping"
else
	echo "$(date)"
	echo mapping $R1
	echo mapping $R2

	bowtie2 \
	-1 $R1 \
	-2 $R2 \
	-S $output_dir/$sample.sam \
	-q \
	--local \
	$a_mapping \
	-p $threads \
	-x $database

	echo "$(date)"
	echo -e "DONE Mapping $sample of $group Group" "\n"
fi

