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
#
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

usage : $0 <-i input_directory> <-d config_files_directory> <-s sample> <-g <group> <-o <output_directory> [-l <log_file>] [-V] [-c] [-v] [-h]

	-i input directory containing files to represent
	-d directory containing config files
	-s sample
	-g group
	-l log file
	-o output directory to create config and pictures
	-c clean: remove config files
	-v version
	-V vervose
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


#DECLARE FLAGS AND VARIABLES

cwd="$(pwd)"
clean=false
vervose=false

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:m:o:g:l:s:d:cVvh"
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
		  	vervose=true
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
mappedDir=$input_dir"/mapping"

echo -e "\n#Executing" $0 "\n"

cdsDdbb_file=$input_dir/database/$sample".gff.bed"

circos_conf_summary="$config_dir/circos_summary.conf"
circos_conf_individual="$config_dir/circos_individual.conf"
circosDir="$output_dir"

plasmidMapped=$mappedDir/$sample".coverage_adapted_clustered_ac"

karyotype_file_individual=$imageDir/$sample".karyotype_individual.txt"
karyotype_file_summary=$imageDir/$sample".karyotype_summary.txt"
annotation_text_file=$imageDir/pID_text_annotation.coordinates
annotation_highlights_file=$imageDir/pID_highlights.conf

coverage_file=$imageDir/$sample".bedgraph_term"
cds_contig_file=$imageDir/$sample".gff.coordinates"

contig_file=$imageDir/$sample".plasmids.bed"
contig_file_complete=$imageDir/$sample".plasmids.complete"
links_file=$imageDir/$sample".plasmids.links"

imageName=$sample"_summary.png"

mkdir -p $circosDir


echo "Creating config file for circos in SAMPLE $sample FILE $circosDir/$sample.circos.conf"

awk '{gsub("PLASMID_KARYOTYPE","'$karyotype_file_individual'"); \
gsub("PLASMID_SPECIFIC_PLASMIDFINDER","'$annotation_text_file'"); \
gsub("PID_ALL_HIGHLIGHTS","'$annotation_highlights_file'"); \
gsub("PLASMID_COVERAGE_GRAPH","'$coverage_file'"); \
gsub("PLASMID_CDS_CONTIG","'$cds_contig_file'"); \
gsub("PLASMID_CDS_DDBB","'$cdsDdbb_file'"); \
gsub("PLASMID_CONTIGS_COMPLETE","'$contig_file_complete'"); \
gsub("PLASMID_CONTIGS","'$contig_file'"); \
gsub("OUTPUTDIR","'$circosDir'"); \
print $0}' $circos_conf_individual > $circosDir/$sample"_individual.circos.conf"

echo "DONE Creating config file for circos in SAMPLE $sample"

echo "Executing circos in SAMPLE $sample"

if [ -f $log_file ]; then
	rm $log_file
fi


for i in $(cat $plasmidMapped)
do
	echo "Creating image for plasmid" $i "in sample" $sample
	awk '{gsub("SAMPLE_SHOWN","'$i'"); \
	gsub("IMAGENAME_SAMPLE_PLASMID","'$sample'_'$i'.png"); \
	print $0}' $circosDir/$sample"_individual.circos.conf" > $circosDir/$sample"_"$i"_individual.circos.conf"
	if [ $vervose = true ]; then
		circos -conf $circosDir/$sample"_"$i"_individual.circos.conf" |& tee -a $log_file
	else
		circos -conf $circosDir/$sample"_"$i"_individual.circos.conf" &>> $log_file
	fi
done 


if [ -s $karyotype_file_summary ]; then

	echo "Creating summary image for in sample" $sample

	awk '{gsub("PLASMID_KARYOTYPE","'$karyotype_file_summary'"); \
	gsub("PLASMID_SPECIFIC_PLASMIDFINDER","'$annotation_text_file'"); \
	gsub("PID_ALL_HIGHLIGHTS","'$annotation_highlights_file'"); \
	gsub("PLASMID_COVERAGE_GRAPH","'$coverage_file'"); \
	gsub("PLASMID_CDS_CONTIG","'$cds_contig_file'"); \
	gsub("PLASMID_CDS_DDBB","'$cdsDdbb_file'"); \
	gsub("PLASMID_CONTIGS","'$contig_file'"); \
	gsub("PLASMID_LINKS","'$links_file'"); \
	gsub("OUTPUTDIR","'$circosDir'"); \
	gsub("IMAGENAME","'$imageName'"); \
	print $0}' $circos_conf_summary > $circosDir/$sample"_summary.circos.conf"

	if [ $vervose = true ]; then
		circos -conf $circosDir/$sample"_summary.circos.conf" |& tee -a $log_file
	else
		circos -conf $circosDir/$sample"_summary.circos.conf" &>> $log_file
	fi
else

	echo "No plasmid mathed requirements to draw the summary image"

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

echo "DONE, files can be found at" $circosDir