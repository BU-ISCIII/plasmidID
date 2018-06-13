#!/bin/bash

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
#DESCRIPTION:Script that uses bedtool to obtain coverage data from a BAMm file
#The default output format is as follows:
#
#chromosome (or entire genome)
#depth of coverage from features in input file
#number of bases on chromosome (or genome) with depth equal to column 2.
#size of chromosome (or entire genome) in base pairs
#fraction of bases on chromosome (or entire genome) with depth equal to column 2.
#
#chr1   0  980  1000  0.98
#chr1   1  20   1000  0.02
#chr2   1  500  500   1
#genome 0  980  1500  0.653333
#genome 1  520  1500  0.346667
#
#-p option is equivalent to -bga BEDGRAPH output
#
#chr1  0       554304  0
#chr1  554304  554309  5
#chr1  554309  554313  6
#chr1  554313  554314  1
#chr1  554314  554315  0
#chr1  554315  554316  6
#chr1  554316  554317  5
#chr1  554317  554318  1
#chr1  554318  554319  2
#chr1  554319  554321  6
#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

Get_coverage script uses bedtool to obtain coverage data from a BAMm file

usage : $0 <-i inputfile(sorted.bam)> [-o <directory>] [-d <database(fasta)>] [-s sample_name]
		 [-g group_name] [-m <int>] [p] [-v] [-h]

	-i input file in sorted BAM format
	-o output directory (optional)
	-d database to extract length. Fasta file used to map against
	-m max depth reported (default 500)
	-p reports genome coverage for all positions in BEDGRAPH format includig 0 positions.
		Default option is bedtools genomecov that needs the reference genome
	-s sample name
	-g group name (optional). If unset, samples will be gathered in NO_GROUP group
	-v version
	-h display usage message

example: get_coverage.sh -i ecoli.bam -d database.fasta
		 get_coverage.sh -i ecoli.bam -p -m 100

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
group="NO_GROUP"
input_file="Input_file"
database="Database"
positional=false
max_coverage=500

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:o:d:s:g:m:n:pvh"
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
		d )
			database=$OPTARG
			;;
		m )
			max_coverage=$OPTARG
			;;
		p )			
          	positional=true
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

bash lib/check_dependencies.sh bedtools

if [ ! $output_dir ]; then
	output_dir=$(dirname $input_file)
	echo "Default output directory is" $output_dir
	mkdir -p $output_dir
else
	echo "Output directory is" $output_dir
	mkdir -p $output_dir
fi

if [ ! $filename ]; then
	filename=$(basename $input_file | cut -d. -f1)
fi



if [ $positional = true ]; then 
	if [ -f $imageDir/$sample".plasmid.bedgraph" ];then \
		echo "Found a bedgraph file for sample" $sample;
		echo "Omitting bedgraph step"
	else
		echo "$(date)"
		echo "Obtaining coverage coordinates from sequences"

		bedtools genomecov -ibam $input_file -bga -max $max_coverage > $output_dir/$filename".bedgraph"

		echo "$(date)"
		echo "DONE obtaining coverage coordinates from sequences"
	fi
else


	bash lib/check_mandatory_files.sh $database

	if [ -f $database".length" ]; then
		echo "Found length file for" $(basename $database)
		echo "Omitting length calculation"
	else
		echo "$(date)"
		echo "Creating a length file for" $(basename $database)
		bash lib/calculate_seqlen.sh -r -i $database > $database".length"
	fi

	if [ -f $output_dir/$filename".coverage" ];then \
		echo "Found a coverage file for sample" $sample;
		echo "Omitting coverage calculation"
	else
		echo "$(date)"
		echo "Calculating coverage for every position that mapped $filename"

		bedtools genomecov -ibam $input_file -g $database".length" > $output_dir/$filename".coverage"

		echo "$(date)"
		echo "DONE Calculating coverage for every plamid that mapped $sample"
	fi
fi

echo -e "\n"




























