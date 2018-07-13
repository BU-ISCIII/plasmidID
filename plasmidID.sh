#!/bin/bash

# Exit immediately if a pipeline, which may consist of a single simple command, a list, 
#or a compound command returns a non-zero status: If errors are not handled by user
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion.

#Print everything as if it were executed, after substitution and expansion is applied: Debug|log option
#set -x

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.3.0
#CREATED: 15 March 2018
#
#ACKNOLEDGE: longops2getops.sh: https://gist.github.com/adamhotep/895cebf290e95e613c006afbffef09d7
#
#DESCRIPTION: plasmidID is a computational pipeline tha reconstruct and annotate the most likely plasmids present in one sample		
#
#
#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

plasmidID is a computational pipeline tha reconstruct and annotate the most likely plasmids present in one sample

usage : $0 <-1 R1> <-2 R2> <-d database(fasta)> <-s sample_name> [-g group_name] [options]

	Mandatory input data:
	-1 | --R1	<filename>	reads corresponding to paired-end R1 (mandatory)
	-2 | --R2	<filename>	reads corresponding to paired-end R2 (mandatory)
	-d | --database	<filename>	database to map and reconstruct (mandatory)
	-s | --sample	<string>	sample name (mandatory)

	Optional input data:
	-g | --group	<string>	group name (optional). If unset, samples will be gathered in NO_GROUP group
	-c | --contig	<filename>	file with contigs. If supplied, plasmidID will not assembly reads
	-a | --annotate <filename>	file with configuration file for specific annotation
	-o 		<output_dir>	output directory, by default is the current directory

	Pipeline options:
	-C | --coverage-cutoff	<int>	minimun coverage percentage to select a plasmid as scafold (0-100), default 80
	-S | --coverage-summary	<int>	minimun coverage percentage to include plasmids in summary image (0-100), default 90
	-f | --cluster		<int>	identity percentage to cluster plasmids into the same representative sequence (0-100), default 80
	-i | --alignment-identity <int>	minimun identity percentage aligned for a contig to annotate, default 90
	-l | --alignment-percentage <int>	minimun length percentage aligned for a contig to annotate, default 30
	-L | --length-total	<int>	minimun alignment length to filter blast analysis

	--explore	Relaxes default parameters to find less reliable relationships within data supplied and database
	--no-trim	Reads supplied will not be quality trimmed
	--only-reconstruct	Database supplied will not be filtered and all sequences will be used as scaffold [NOT TESTED]
	
	Additional options:

	-M | --memory	<int>	max memory allowed to use
	-T | --threads	<int>	number of threads
	-v | --version		version
	-h | --help		display usage message

example: ./plasmidID.sh -1 ecoli_R1.fastq.gz -2 ecoli_R2.fastq.gz -d database.fasta -s ECO_553 -G ENTERO
		./plasmidID.sh -1 ecoli_R1.fastq.gz -2 ecoli_R2.fastq.gz -d PacBio_sample.fasta -c scaffolds.fasta -C 60 -s ECO_60 -G ENTERO --no-trim

EOF
}

#================================================================
# OPTION_PROCESSING
#================================================================
#Make sure the script is executed with arguments
if [ $# = 0 ]; then
	echo "NO ARGUMENTS SUPPLIED"
	usage >&2
	exit 1
fi


# translate long options to short
reset=true
for arg in "$@"
do
    if [ -n "$reset" ]; then
      unset reset
      set --      # this resets the "$@" array so we can rebuild it
    fi
    case "$arg" in
##MANDATORY MINIMAL OPTIONS
    	--R1)		set -- "$@"	-1 ;;
		--R2)		set -- "$@"	-2 ;;
		--database)	set -- "$@"	-d ;;
		--sample)	set -- "$@"	-s ;;
##OPTIONAL
		--group)	set -- "$@"	-g ;;
		--contig)	set -- "$@"	-c ;;
		--annotate)	set -- "$@"	-a ;;
##PIPELINE
		--only-reconstruct)	set -- "$@" -R ;;
		--no-trim)	set -- "$@" -t ;;
		--trimmomatic_directory) set -- "$@" -X ;;
		--explore)	set -- "$@" -E ;;
		--alignment-identity)	set -- "$@" -i ;;

		--coverage-cutoff)	set -- "$@"	-C ;;
		--coverage-summary)	set -- "$@"	-S ;;
		--cluster)	set -- "$@"	-f ;;
		--alignment-percentage)	set -- "$@"	-l ;;
		--length-total)	set -- "$@"	-L ;;
##ADDITIONAL-
		--memory)		set -- "$@"	-M ;;
		--threads)		set -- "$@"	-T ;;
       	--help)    	set -- "$@" -h ;;
		--vervose) 	set -- "$@" -V ;;
       	--version) 	set -- "$@" -v ;;
       # pass through anything else
       *)         set -- "$@" "$arg" ;;
    esac
done

#DECLARE FLAGS AND VARIABLES
threads=1
max_memory=4000
cwd="$(pwd)"
group="NO_GROUP"
database="Database"
coverage_cutoff=80
coverage_summary=80
cluster_cutoff=80
alignment_identity=90
alignment_percentage=30
R1="R1_file"
R2="R2_file"
no_trim=false
only_reconstruct=false
explore=false
include_assembly=true
annotation=false
vervose_option_circos=""
is_vervose=false

#SET COLORS

YELLOW='\033[0;33m'
WHITE='\033[0;37m'
CYAN='\033[0;36m'
BLUE='\033[0;34m'
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":1:2:d:s:g:c:a:i:o:C:S:f:l:L:T:M:X:RVtvh"
while getopts $options opt; do
	case $opt in
		1 )
			r1_file=$OPTARG
			;;
		2 )
			r2_file=$OPTARG
			;;
		d )
			database=$OPTARG
			;;
		s )
			sample=$OPTARG
			;;
		g)
			group=$OPTARG
			;;
		c)
			contigs=$OPTARG
			include_assembly=false
			;;
		g)
			group=$OPTARG
			;;
		a )
			annotation_config_file=$OPTARG
			annotation=true
			;;
		t)
			no_trim=true
			;;
		X)
			trimmomatic_directory=$OPTARG
			;;
		E)
			explore=true
			;;
		R )
			only_reconstruct=true
			reconstruct_fasta=$database
			;;
		C )
			coverage_cutoff=$OPTARG
			;;
		S )
			coverage_summary=$OPTARG
			;;
		f )			
          	cluster_cutoff=$OPTARG
      		;;
      	i)
			alignment_identity=$OPTARG
			;;
      	l )			
          	alignment_percentage=$OPTARG
      		;;
      	L )			
          	alignment_total=$OPTARG
      		;;
        T ) 
			threads=$OPTARG
            ;;
        M)
			max_memory=$OPTARG
			;;
        o)
			output_dir=$OPTARG
			;;
		V)
			vervose_option_circos="-V"
			is_vervose=true
			log_file="/dev/stdout"
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
echo -e "CHECHKING DEPENDENCIES AND MANDATORY FILES"
check_dependencies.sh blastn bowtie2-build bowtie2 cd-hit-est bedtools prokka samtools circos

check_mandatory_files.sh $r1_file $r2_file $database 



if [ ! $sample ]; then
	echo "ERROR: please, provide a sample name"
	usage
	exit 1
fi

if [ ! $output_dir ]; then
	output_dir=$cwd
	echo "Default output directory is" $output_dir
	mkdir -p $output_dir
else
	echo "Output directory is" $output_dir
	mkdir -p $output_dir
fi

mkdir -p $output_dir/$group/$sample/logs

if [ $is_vervose = false ]; then
	log_file=$output_dir/$group/$sample/logs/plasmidID.log
	echo $log_file
	echo "LOG FILE PLASMIDID" > $log_file
	echo $(date) >> $log_file
fi


reconstruct_fasta=$output_dir/$group/$sample/mapping/$sample".coverage_adapted_filtered_"$coverage_cutoff"_term.fasta"_$cluster_cutoff

if [ $explore = true ]; then
	coverage_cutoff=60
	coverage_summary=70
	cluster_cutoff=70
	alignment_percentage=10
fi

echo ""

if [ $annotation = true ]; then
	number_of_annotation=$(cat $annotation_config_file | awk '!/^#/ &&  NF > 0' | wc -l)
	lines_to_annotate=$(cat $annotation_config_file | awk '!/^#/ &&  NF > 0 {print NR}')
	echo "Annotation requested"
	echo "Number of databases to annotate:" $number_of_annotation
	database_number=0

	for i in $lines_to_annotate
	do

		let database_number=database_number+1
		#check if config file exist with highlights and text

		#Declare variables
		ddbb_file=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $1}')
		ddbb_name=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $2}')
		percent_identity=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $3}')
		percent_aligment=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $4}')
		query_divisor=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $5}')
		query_side=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $6}')
		is_unique=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $7}')
		double_unique=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $8}')
		database_type=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $9}')
		color_highlight=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $10}')


		#Check each variable and its position in order to rectify the correct one
		if [[ "$percent_identity" -lt 0  &&  "$percent_identity" -gt 100 ]]; then
			echo -e "${RED}ERROR${NC}:Invalid value for P_IDENTITY. Please check $annotation_config_file\n"
			exit 1
		fi

		if [[ $percent_aligment -lt 0  &&  $percent_aligment -gt 100 ]]; then
			echo -e "${RED}ERROR${NC}:Invalid value for P_ALIGNMENT. Please check $annotation_config_file\n"
			exit 1
		fi

		counter_variable=0
		for i in "$ddbb_file" "$ddbb_name" "$percent_identity" "$percent_aligment" "$query_divisor" "$query_side" "$is_unique" "$double_unique" "$color_highlight"
		do
			let counter_variable=counter_variable+1
			if [[ "$i" != "$query_divisor" ]];then
			#echo $i " " $counter_variable
				if [ -z "$i" -a "$counter_variable" -eq 9 ];then
					echo "${RED}ERROR${NC}:An input filed is missing in annotation" $database_number "please check annotation file with all fields separated by comma"
					exit 1
				elif [[ -z "$i" ]];then
					echo "${RED}ERROR${NC}:Input" $counter_variable "in annotation" $database_number "is missing, please check annotation file with all fields separated by comma"
					exit 1
				fi
			fi
		done
	done
fi


command_log=$output_dir/$group/$sample/logs/00_commands.log

if [ -f $command_log ];then
	rm $command_log
fi


####START PLASMIDID########################################################################################################################
###########################################################################################################################################

if [ $only_reconstruct = false ]; then

	reconstruct_fasta=$output_dir/$group/$sample/mapping/$sample".coverage_adapted_filtered_"$coverage_cutoff"_term.fasta"_$cluster_cutoff

####TRIMMING#################################################################
#############################################################################
	
	if [ $no_trim = false ]; then

		r1_file_mapping=$r1_file
		r2_file_mapping=$r2_file

		echo -e "\n${CYAN}TRIMMING READS{NC}\n"
		echo "R1:" $r1_file
		echo "R2:" $r2_file

		filestrimmed=$(ls $output_dir/$group/$sample/trimmed 2> /dev/null | grep "paired" | wc -l)

		if [ "$filestrimmed" = "4" ];then
			echo "Found trimmed files for this sample" $sample;
			echo "Omitting trimming"
		else
			echo "####TRIMMING################################################################"
			quality_trim.sh -1 $r1_file -2 $r2_file -s $sample -g $group -o $output_dir/$group/$sample/trimmed -d $trimmomatic_directory -T $threads
		fi
	else
		echo "No trim selected, skipping trimming step"
		r1_file_mapping=$r1_file
		r2_file_mapping=$r2_file
	fi

	#group/sample/trimmed/sample_1_paired.fastq.gz
	#group/sample/trimmed/sample_1_unpaired.fastq.gz
	#group/sample/trimmed/sample_2_paired.fastq.gz
	#group/sample/trimmed/sample_2_unpaired.fastq.gz

####ASSEMLY##################################################################
#############################################################################

	if [ $include_assembly = true ]; then
		
		if [ -f $output_dir/$group/$sample/assembly/scaffolds.fasta ]; then
			echo "Found an assembled file for sample" $sample
			echo "Omitting assembly"
			contigs=$output_dir/$group/$sample/assembly/scaffolds.fasta
		else
			if [ $no_trim = true ]; then
				check_dependencies.sh spades.py
				echo "####ASSEMBLY################################################################"
				spades_assembly.sh -p $r1_file_mapping -P $r2_file_mapping -c -T $threads -o $output_dir/$group/$sample/assembly &>> $log_file
			else
				check_dependencies.sh spades.py
				echo "####ASSEMBLY################################################################"
				spades_assembly.sh -q $output_dir/$group/$sample/trimmed/ -c -T $threads &>> $log_file
			fi
			contigs=$output_dir/$group/$sample/assembly/scaffolds.fasta
		fi
	else
		echo "Contigs supplied, ommiting assembly step"
	fi

	#group/sample/assembly/scaffolds.fasta

####MAPPING##################################################################
#############################################################################
	if [ -f  $output_dir/$group/$sample/mapping/$sample.sorted.bam -a -f  $output_dir/$group/$sample/mapping/$sample.sorted.bam.bai -o -f  $output_dir/$group/$sample/mapping/$sample.sam ];then
		echo "Found a mapping file for sample" $sample;
		echo "Omitting mapping"
	else
		echo "####MAPPING#################################################################"
		bowtie_mapper.sh -a -d $database \
		-T $threads \
	 	-g $group \
	 	-s $sample \
	 	-1 $r1_file_mapping \
	 	-2 $r2_file_mapping \
		-o $output_dir/$group/$sample/mapping &>> $log_file

 #group/sample/mapping/sample.sam

	 	sam_to_bam.sh -i  $output_dir/$group/$sample/mapping/$sample.sam

 #group/sample/mapping/sample.bam
 #group/sample/mapping/sample.bam.bai
 	fi

####COVERAGE FILTERING#######################################################
#############################################################################
	if [ -f  $output_dir/$group/$sample/mapping/$sample".coverage" ];then \
		echo "Found a coverage file for sample" $sample;
		echo "Omitting coverage calculation"
	else
		echo "####COVERAGE FILTERING########################################################"
 		get_coverage.sh -i  $output_dir/$group/$sample/mapping/$sample".sorted.bam" -d $database
 	fi
 #group/sample/mapping/sample.coverage


	temporary_coverage_files=$(ls  $output_dir/$group/$sample/mapping/$sample".coverage_adapted_filtered_??" 2> /dev/null | wc -l)

 	if [ $temporary_coverage_files -gt 1 ]; then

		echo "Removing previous coverage filtered files"
		for i in $(ls  $output_dir/$group/$sample/mapping/$sample".coverage_adapted_filtered"*)
		do
			rm $i
		done
	fi

 	adapt_filter_coverage.sh -i  $output_dir/$group/$sample/mapping/$sample".coverage" -c $coverage_cutoff

 #sample.coverage_adapted
 #sample.coverage_adapted_filtered_80

	filter_fasta.sh -i $database -f  $output_dir/$group/$sample/mapping/$sample".coverage_adapted_filtered_"$coverage_cutoff

#sample.coverage_adapted_filtered_50_term.fasta

	
	if [ ! -s  $output_dir/$group/$sample/mapping/$sample".coverage_adapted_filtered_"$coverage_cutoff ]; then \
		echo -e "${RED}ERROR${NC}"
		echo "NO PLASMIDS MATCHED MAPPING REQUERIMENTS, PLEASE, TRY WITH LOWER COVERAGE CUTOFF"
		echo "################################################################################"
		exit 1
	fi



	if [ -f  $output_dir/$group/$sample/mapping/$sample".coverage_adapted_filtered_"$coverage_cutoff"_term.fasta_"$cluster_cutoff ];then \
		echo "Found a clustered file for sample" $sample;
		echo "Omitting clustering"
	else

		cdhit_cluster.sh -i  $output_dir/$group/$sample/mapping/$sample".coverage_adapted_filtered_"$coverage_cutoff"_term.fasta" -c $cluster_cutoff -M $max_memory -T $threads
	fi



	process_cluster_output.sh -i $reconstruct_fasta -b  $output_dir/$group/$sample/mapping/$sample".coverage_adapted" -c $coverage_cutoff
#$group/$sample/mapping/$sample".coverage_adapted_filtered_"$coverage_cutoff"_term.fasta"_$cluster_cutoff >> FINAL CLUSTERED AND FILTERED FASTA FILE TO USE AS SCAFFOLD
########################################################################################################

#sample.coverage_adapted_filtered_50_term.fasta_80
#sample.coverage_adapted_filtered_50_term.fasta_80.clstr

fi


##sample.coverage_adapted_clustered
##sample.coverage_adapted_clustered_percentage
##sample.coverage_adapted_clustered_ac
############################################################

###############################################################################################################################################
######################## DONE WITH MAPPING AND CLUSTERING #####################################################################################
###############################################################################################################################################

echo -e "\n${YELLOW} STARTING CONTIG ALIGNMENT and ANNOTATON${NC}\n"

echo -e "\n${CYAN}OBTAINING KARYOTYPE TRACKS${NC}\n \
A file with the informatin of putative plasmid and its length will be generated.\n"

build_karyotype.sh -i  $output_dir/$group/$sample/mapping/$sample".coverage_adapted_clustered" -K $coverage_summary -k $coverage_cutoff -o  $output_dir/$group/$sample/data &>> $log_file

#group/sample/data
#sample.karyotype_individual.txt
#sample.karyotype_individual.txt

echo -e "\n${CYAN}OBTAINING COVERAGE TRACK${NC}\n \
A bedgraph file containing mapping information for filtered plasmids will be generated.\n"

if [ -f  $output_dir/$group/$sample/data/$sample".bedgraph" ];then
	echo "Found a coverage file for sample" $sample;
	echo "Omitting coverage calculation"
else
	get_coverage.sh -i  $output_dir/$group/$sample/mapping/$sample".sorted.bam" -p -o  $output_dir/$group/$sample/data &>> $log_file
fi
#sample.bedgraph

filter_fasta.sh -i  $output_dir/$group/$sample/data/$sample".bedgraph" -f  $output_dir/$group/$sample/mapping/$sample".coverage_adapted_clustered_ac" -G &>> $log_file

#sample.bedgraph_term

echo -e "\n${CYAN}ANNOTATING CONTIGS${NC}\n \
A file including all automatic annotations on contigs will be generated.\n"

if [ -f  $output_dir/$group/$sample/data/$sample".fna" -a -f  $output_dir/$group/$sample/data/$sample".gff" ];then
	echo "Found annotations files for sample" $sample;
	echo "Omitting automatic annotation"
else
	echo "prokka_annotation.sh -i $contigs -p $sample -T $threads -o  $output_dir/$group/$sample/data -c &>> $log_fil" >> $command_log
	prokka_annotation.sh -i $contigs -p $sample -T $threads -o  $output_dir/$group/$sample/data -c &>> $log_file
fi
#sample.fna
#sample.gff


echo -e "\n${CYAN}ALIGNING CONTIGS TO FILTERED PLASMIDS${NC}\n \
Contigs are aligned to filtered plasmids and those are selected by alignment identity and alignment percentage \
in order to create links, full length and annotation tracks\n"

if [ -f $output_dir/$group/$sample/logs/contig_alignment.log ]; then
	rm $output_dir/$group/$sample/logs/contig_alignment.log
fi

blast_align.sh -i  $output_dir/$group/$sample/data/$sample".fna" -d $reconstruct_fasta -o  $output_dir/$group/$sample/data -p plasmids &>> $log_file

#sample.plasmids.blast


blast_to_bed.sh -i  $output_dir/$group/$sample/data/$sample".plasmids.blast" -b $alignment_identity -l 0 -L 500 -d - -q _ -Q r -I &>> $log_file

#sample.plasmids.bed

blast_to_complete.sh -i  $output_dir/$group/$sample/data/$sample".plasmids.blast" -l $alignment_percentage &>> $log_file

#sample.plasmids.complete

blast_to_link.sh -i  $output_dir/$group/$sample/data/$sample".plasmids.blast" -I -l $alignment_percentage &>> $log_file

#sample.plasmids.links
#sample.plasmids.blast.links

gff_to_bed.sh -i  $output_dir/$group/$sample/data/$sample".gff" -u -L &>> $log_file

#sample.gff.bed

coordinate_adapter.sh -i  $output_dir/$group/$sample/data/$sample".gff.bed" -l  $output_dir/$group/$sample/data/$sample".plasmids.blast.links" -p -n 1500 &>> $log_file

#sample.gff.coordinates


if [ $annotation = true ]; then

	echo -e "\n${CYAN}ANNOTATING SPECIFIC DATABASES SUPPLIED${NC}\n \
	Each database supplied will be locally aligned against contigs and the coordinates will be adapted for image representation\n"


	if [ -f $output_dir/$group/$sample/logs/speciffic_annotation.log ]; then
		rm $output_dir/$group/$sample/logs/speciffic_annotation.log
	fi

	number_of_annotation=$(cat $annotation_config_file | awk '!/^#/ &&  NF > 0' | wc -l)
	lines_to_annotate=$(cat $annotation_config_file | awk '!/^#/ &&  NF > 0 {print NR}')
	echo "Number of databases to annotate:" $number_of_annotation
	database_number=0

	if [ -f $output_dir/$group/$sample/data/pID_highlights.conf ]; then
		rm $output_dir/$group/$sample/data/pID_highlights.conf
	fi

	if [ -f $output_dir/$group/$sample/data/pID_text_annotation.coordinates ]; then
		rm $output_dir/$group/$sample/data/pID_text_annotation.coordinates
	fi

	for i in $lines_to_annotate
	do

		let database_number=database_number+1
		
		#check if config file exist with highlights and text

		#Declare variables
		ddbb_file=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $1}')
		ddbb_name=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $2}')
		percent_identity=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $3}')
		percent_aligment=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $4}')
		query_divisor=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $5}')
		query_side=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $6}')
		is_unique=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $7}')
		double_unique=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $8}')
		database_type=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $9}')
		color_highlight=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $10}')

		echo "Annotating database" $database_number":" $ddbb_name

		if [ "$is_unique" == "y" ];then
			is_unique_command="-u"
		elif [ "$is_unique" == "n" ];then
			is_unique_command=""
		fi

		if [ "$double_unique" == "n" ];then
			double_unique_command=""
		else
			double_unique_command="-U ${double_unique}"
		fi

		#echo "blast_align.sh -i $ddbb_file -d $output_dir/$group/$sample/data/$sample".fna" -o $output_dir/$group/$sample/data -p $ddbb_name -f $sample -t $database_type &>> $log_file >> $command_log
		blast_align.sh -i $ddbb_file -d $output_dir/$group/$sample/data/$sample".fna" -o $output_dir/$group/$sample/data -p $ddbb_name -f $sample -t $database_type &>> $log_file
		#sample.annotation.blast

		#echo "		blast_to_bed.sh -i $output_dir/$group/$sample/data/$sample"."$ddbb_name".blast" -b $percent_identity -l $percent_aligment -d _ -D r -q "$query_divisor" -Q $query_side $double_unique_command &>> $log_file"
		blast_to_bed.sh -i $output_dir/$group/$sample/data/$sample"."$ddbb_name".blast" -b $percent_identity -l $percent_aligment -d _ -D r -q "$query_divisor" -Q $query_side $double_unique_command &>> $log_file
		#sample.annotation.bed

		#echo "coordinate_adapter.sh -i $output_dir/$group/$sample/data/$sample"."$ddbb_name".bed" -l $output_dir/$group/$sample/data/$sample".plasmids.blast.links" $is_unique_command &>> $log_file"
		coordinate_adapter.sh -i $output_dir/$group/$sample/data/$sample"."$ddbb_name".bed" -l $output_dir/$group/$sample/data/$sample".plasmids.blast.links" $is_unique_command &>> $log_file
		#sample.annotation.coordinates


		coordinates_file=$(echo  $output_dir/$group/$sample/data/$sample"."$ddbb_name".coordinates")
		z_value="10${database_number}"


		printf '%s\n' "<highlight>" "file = ${coordinates_file}" "z= ${z_value}" "r1 = 0.90r" "r0 = 0.69r" "fill_color = ${color_highlight}" "</highlight>" >>  $output_dir/$group/$sample/data/pID_highlights.conf
		
		cat $output_dir/$group/$sample/data/$sample"."$ddbb_name".coordinates" >> $output_dir/$group/$sample/data/pID_text_annotation.coordinates

	done

else
	touch $output_dir/$group/$sample/data/pID_highlights.conf
	touch $output_dir/$group/$sample/data/pID_text_annotation.coordinates
	
fi




echo -e "\n${CYAN}ANNOTATING PLASMIDS USED AS SCAFFOLD${NC}\n \
A file including all automatic annotations on plasmids that matched requeriments will be generated.\n"

#Annotate plasmids selected as database
if [ -f $output_dir/$group/$sample/database/$sample".fna" -a -f $output_dir/$group/$sample/database/$sample".gff" ];then
	echo "Found annotations files for reconstruct plasmids i sample" $sample;
	echo "Omitting automatic annotation"
else
	echo "Executing prokka for plasmids from database"

	echo "prokka_annotation.sh -i $reconstruct_fasta -p $sample -o $output_dir/$group/$sample/database -c" >> $command_log
	prokka_annotation.sh -i $reconstruct_fasta -p $sample -T $threads -o $output_dir/$group/$sample/database -c &>> $log_file

	#database/sample.fna
	#database/sample.gff
	echo "rename_from_fasta.sh -i $output_dir/$group/$sample/database/$sample.gff -1 $reconstruct_fasta -2 $output_dir/$group/$sample/database/$sample.fna" >> $command_log

	rename_from_fasta.sh -i $output_dir/$group/$sample/database/$sample".gff" -1 $reconstruct_fasta -2 $output_dir/$group/$sample/database/$sample".fna" &>> $log_file

	#sample.gff.renamed

	echo "gff_to_bed.sh -i $output_dir/$group/$sample/database/$sample.gff.renamed -q " " -u -L" >> $command_log

	gff_to_bed.sh -i $output_dir/$group/$sample/database/$sample".gff.renamed" -q " " -u -L &>> $log_file

	#database/sample.gff.bed
fi




echo -e "\n${CYAN}DRAWING CIRCOS IMAGES${NC}\n \
An image per putative plasmid will be drawn having into account all data supplied.\n \
Additionally a summary image will be created to determine redundancy within remaining plasmids\n"

echo "draw_circos_images.sh -i $output_dir/$group/$sample \
-d config_files \
-o $output_dir/$group/$sample/images \
-g $group -s $sample -l $output_dir/$group/$sample/logs/draw_circos_images.log -c $vervose_option_circos" >> $command_log 

draw_circos_images.sh -i $output_dir/$group/$sample \
-d config_files \
-o $output_dir/$group/$sample/images \
-g $group -s $sample -l $output_dir/$group/$sample/logs/draw_circos_images.log -c $vervose_option_circos


echo -e "\n${YELLOW}ALL DONE WITH plasmidID${NC}\n \
Thank you for using plasmidID\n"
