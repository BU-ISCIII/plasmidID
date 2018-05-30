#!/bin/bash

set -e

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.0 
#CREATED: 21 May 2018
#REVISION:
#DESCRIPTION:Script that execute trimmomatic to filter by quality
#
#
#================================================================
# END_OF_HEADER
#================================================================


usage() {
	cat << EOF

quality_trim script execute trimmomatic to filter by quality

usage : $0 <-1 R1 file> <-2 R2 file> [-o <directory>] [-d <trimmomatic_directory>] <-s sample_name>
		[-a adapter_file] [-g group_name] [-f <file_name>] [-l <int>] [-M <int>] [-T <int>][-v] [-h]

	-1 R1 file (mandatory)
	-2 R2 file (mandatory)
	-d directory where trimmomatic is installed, default: /opt/Trimmomatic/
	-a adapters to remove, default: TruSeq3-PE.fa
	-o output directory (optional)
	-f file name
	-l minimus length of trimmed reads (default 40)
	-s sample name (mandatory)
	-g group name (optional). If unset, samples will be gathered in NO_GROUP group
	-M RAM memmory (Gb), default 8
	-T threads, default 1
	-v version
	-h display usage message

example: ./quality_trim.sh -1 ecoli_R1.fastq.gz -2 ecoli_R2.fastq.gz -s ECO232 -g ENTERO -T 8

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
r1_file="R1_file"
r2_file="R2_file"
trimmomatic_directory=/opt/Trimmomatic/
adapter_file="TruSeq3-PE.fa"
minimus_length=40
max_mem=8
threads=1

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":1:2:o:f:d:a:s:g:l:n:M:T:vh"
while getopts $options opt; do
	case $opt in
		1 )
			r1_file=$OPTARG
			;;
		2 )
			r2_file=$OPTARG
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
		d)
			trimmomatic_directory=$OPTARG
			;;
		a)
			adapter_file=$OPTARG
			;;
		l)
			minimus_length=$OPTARG
			;;
		g)
			group=$OPTARG
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

bash lib/check_mandatory_files.sh $r1_file $r2_file

bash lib/check_dependencies.sh java

if [ ! $sample ]; then
	echo "Please include a sample name"
	exit 1
fi


if [ ! $output_dir ]; then
	output_dir="$group/$sample/trimmed"
	echo "Default output directory is" $output_dir
	mkdir -p $output_dir
else
	echo "Output directory is" $output_dir
	mkdir -p $output_dir
fi

if [ ! $filename ]; then
	filename=$sample
fi


trimmomatic_executable=$(find $trimmomatic_directory -type f -name "trimmomatic*.jar" | awk 'NR==1')
trimmomatic_adapter=$(find $trimmomatic_directory -type f -name $adapter_file | awk 'NR==1')


echo "$(date)"
echo "Quality trimming:"
echo "R1 = " $r1_file
echo "R2 = " $r2_file

java -jar "-Xmx"$max_mem"G" $trimmomatic_executable PE -threads $threads \
$r1_file \
$r2_file \
$output_dir/$sample"_1_paired.fastq.gz" \
$output_dir/$sample"_1_unpaired.fastq.gz" \
$output_dir/$sample"_2_paired.fastq.gz" \
$output_dir/$sample"_2_unpaired.fastq.gz" \
ILLUMINACLIP:$trimmomatic_adapter:2:30:10 SLIDINGWINDOW:4:20 MINLEN:$minimus_length

echo "$(date)"
echo "DONE quality trimming, file can be fount at:"
echo $output_dir/$sample"_1_paired.fastq.gz"
echo $output_dir/$sample"_1_unpaired.fastq.gz"
echo $output_dir/$sample"_2_paired.fastq.gz"
echo $output_dir/$sample"_2_unpaired.fastq.gz"
echo -e "\n"





























