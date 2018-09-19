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
#CREATED: 01 May 2018
#REVISION:
#		11 July 2018: Apply good practices bash
#						Include independent files
#						Include several databases
#		13 July 2018: Include log file
#						manage directories
#DESCRIPTION:Script that creates and execute a cicos config file for plasmidID
#
#
#
#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

draw_circos_image script that creates and execute a cicos config file for plasmidID

usage : $0 <-i input_directory> <-d config_files_directory> <-C config_file> <-s sample> <-g <group> <-o <output_directory> [-l <log_file>] [-V] [-c] [-v] [-h]

	-i input directory containing files to represent
	-d directory containing config files
	-C config file selected to draw
	-s sample
	-g group
	-l log file
	-o output directory to create config and pictures
	-c clean: remove config files
	-v version
	-V verbose
	-h display usage message

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
clean=false
verbose=false


#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:m:o:g:l:s:d:C:cVvh"
while getopts $options opt; do
	case $opt in
		i )
			input_dir=$OPTARG
			;;
		o )
			output_dir=$OPTARG
			;;
		d )
			config_dir=$OPTARG
			;;
		C )
			config_file_individual=$OPTARG
			;;
		l )
			log_file=$OPTARG
			;;
		g )
			group=$OPTARG
			;;
		s )
			sample=$OPTARG
			;;
		c )
			clean=true
			;;
        h )
		  	usage
		  	exit 1
		  	;;
		V )
		  	verbose=true
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

imageDir=$input_dir"/data"

if [ -f $log_file ]; then
	rm $log_file
fi

echo -e "\n#Executing" $0 "\n" &>> $log_file

cdsDdbb_file=$input_dir/database/$sample".gff.bed"
cdsDdbb_file_forward=$input_dir/database/$sample".gff.forward.bed"
cdsDdbb_file_reverse=$input_dir/database/$sample".gff.reverse.bed"


circos_conf_summary="$config_dir/circos_summary_1_3_3.conf"
circos_conf_individual="$config_dir/$config_file_individual"
circosDir="$output_dir"


plasmidMapped=$imageDir/$sample".coverage_adapted_clustered_ac"

karyotype_file_individual=$imageDir/$sample".karyotype_individual.txt"
karyotype_file_summary=$imageDir/$sample".karyotype_summary.txt"
annotation_text_file=$imageDir/pID_text_annotation.coordinates
annotation_highlights_file=$imageDir/pID_highlights.conf

coverage_file=$imageDir/$sample".bedgraph_term"
cds_contig_file=$imageDir/$sample".gff.coordinates"
cds_contig_file_forward=$imageDir/$sample".gff.forward.coordinates"
cds_contig_file_reverse=$imageDir/$sample".gff.reverse.coordinates"


contig_file=$imageDir/$sample".plasmids.bed"
contig_file_complete=$imageDir/$sample".plasmids.complete"
links_file=$imageDir/$sample".plasmids.links"

imageName=$sample"_summary.png"

mkdir -p $circosDir


echo "Creating individual config file for SAMPLE $sample using FILE $circos_conf_individual" &>> $log_file

awk '{gsub("PLASMID_KARYOTYPE","'$karyotype_file_individual'"); \
gsub("PLASMID_SPECIFIC_TEXT","'$annotation_text_file'"); \
gsub("PID_ALL_HIGHLIGHTS","'$annotation_highlights_file'"); \
gsub("PLASMID_COVERAGE_GRAPH","'$coverage_file'"); \
gsub("PLASMID_CDS_CONTIG","'$cds_contig_file'"); \
gsub("PLASMID_CDS_FORWARD","'$cds_contig_file_forward'"); \
gsub("PLASMID_CDS_REVERSE","'$cds_contig_file_reverse'"); \
gsub("PLASMID_CDS_DDBB","'$cdsDdbb_file'"); \
gsub("CDS_DDBB_FORWARD","'$cdsDdbb_file_forward'"); \
gsub("CDS_DDBB_REVERSE","'$cdsDdbb_file_reverse'"); \
gsub("PLASMID_CONTIGS_COMPLETE","'$contig_file_complete'"); \
gsub("PLASMID_CONTIGS","'$contig_file'"); \
gsub("OUTPUTDIR","'$circosDir'"); \
print $0}' $circos_conf_individual > $circosDir/$sample"_individual.circos.conf"

echo "DONE Creating config file for circos in SAMPLE $sample" &>> $log_file

echo "Executing circos in SAMPLE $sample" &>> $log_file



for i in $(cat $plasmidMapped)
do
	echo "Creating image for plasmid $i in sample $sample" &>> $log_file
	awk '{gsub("SAMPLE_SHOWN","'$i'"); \
	gsub("IMAGENAME_SAMPLE_PLASMID","'$sample'_'$i'.png"); \
	print $0}' $circosDir/$sample"_individual.circos.conf" > $circosDir/$sample"_"$i"_individual.circos.conf"
	if [ $verbose = true ];then
		$(circos -conf $circosDir/$sample"_"$i"_individual.circos.conf" |& tee -a $log_file) || error ${LINENO} $(basename $0) "Circos command for individual image has failed. See $output_dir/logs for more information"
	else
		$(circos -conf $circosDir/$sample"_"$i"_individual.circos.conf" &>> $log_file) || error ${LINENO} $(basename $0) "Circos command for individual image has failed. See $output_dir/logs for more information"
	fi
done


if [ -s $karyotype_file_summary ]; then

	echo "Creating summary image for in sample" $sample "from FILE" $circos_conf_summary &>> $log_file

	awk '{gsub("PLASMID_KARYOTYPE","'$karyotype_file_summary'"); \
	gsub("PLASMID_SPECIFIC_TEXT","'$annotation_text_file'"); \
	gsub("PID_ALL_HIGHLIGHTS","'$annotation_highlights_file'"); \
	gsub("PLASMID_COVERAGE_GRAPH","'$coverage_file'"); \
	gsub("PLASMID_CDS_CONTIG","'$cds_contig_file'"); \
	gsub("PLASMID_CDS_FORWARD","'$cds_contig_file_forward'"); \
	gsub("PLASMID_CDS_REVERSE","'$cds_contig_file_reverse'"); \
	gsub("PLASMID_CDS_DDBB","'$cdsDdbb_file'"); \
	gsub("PLASMID_CONTIGS","'$contig_file'"); \
	gsub("PLASMID_LINKS","'$links_file'"); \
	gsub("OUTPUTDIR","'$circosDir'"); \
	gsub("IMAGENAME","'$imageName'"); \
	print $0}' $circos_conf_summary > $circosDir/$sample"_summary.circos.conf"

	if [ $verbose = true ]; then
		circos -conf $circosDir/$sample"_summary.circos.conf" |& tee -a $log_file || exit 1
	else
		circos -conf $circosDir/$sample"_summary.circos.conf" &>> $log_file || exit 1
	fi

else

	echo "No plasmid matched requirements to draw the summary image"

fi


#Remove config files
if [ clean = true ];then
	for i in $(cat $plasmidMapped)
	do
		if [ -f $circosDir/$sample"_"$i"_individual.circos.conf" ]; then
			rm $circosDir/$sample"_"$i"_individual.circos.conf"
		fi
	done

		rm $circosDir/$sample"_summary.circos.conf"
		rm $circosDir/$sample"_individual.circos.conf"
fi

echo "DONE, files can be found at $circosDir"
