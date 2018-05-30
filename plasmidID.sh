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
VERSION=Beta 
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

usage : $0 [-1 <R1>] [-2 <R2>] <-d database(fasta)> <-s sample_name> [-g group_name] [options]

	Mandatory input data:
	-1 | --R1	<filename>	reads corresponding to paired-end R1 (mandatory)
	-2 | --R2	<filename>	reads corresponding to paired-end R2 (mandatory)
	-d | --database	<filename>	database to map and reconstruct (mandatory)
	-s | --sample	<string>	sample name (mandatory)

	Optional input data:
	-g | --group	<string>	group name (optional). If unset, samples will be gathered in NO_GROUP group
	-c | --contig	<filename>	file with contigs. If supplied, plasmidID will not assembly reads
	-a | --annotate <filename>	file with sequences to draw in the final images
	-o 		<output_dir>	output directory, by default is the current directory

	Pipeline options:
	-C | --coverage-cutoff	<int>	minimun coverage percentage to select a plasmid as scafold (0-100), default 80
	-S | --coverage-summary	<int>	minimun coverage percentage to include plasmids in summary image (0-100), default 90
	-f | --cluster		<int>	identity percentage to cluster plasmids into the same representative sequence (0-100), default 80
	-i | --alignment-identity <int>	minimun identity percentage aligned for a contig to annotate
	-l | --alignment-percentage <int>	minimun length percentage aligned for a contig to annotate
	-L | --length-total	<int>	minimun alignment length to filter blast analysis

	--explore	Relaxes default parameters to find less reliable relationships within data supplied and database
	--no-trim	Reads supplied will not be quality trimmed
	--only-reconstruct	Database supplied will not be filtered and all sequences will be used as scaffold 
	
	Additional options:

	-M | --memory	<int>	max memory allowed to use
	-T | --threads	<int>	number of threads
	-v | --version		version
	-h | --help		display usage message

example: ./plasmidID.sh -1 ecoli_R1.fastq.gz -2 ecoli_R2.fastq.gz -d database.fasta -s ECO_553 -G ENTERO
		./plasmidID.sh -1 ecoli_R1.fastq.gz -2 ecoli_R2.fastq.gz -d PacBio_sample.fasta -c scaffolds.fasta -C 60 -s ECO_60 -G ENTERO --only_reconstruct

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
coverage_summary=90
cluster_cutoff=80
alignment_identity=90
alignment_percentage=30
R1="R1_file"
R2="R2_file"
no_trim=false
only_reconstruct=false
explore=false
include_assembly=true

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":1:2:d:s:g:c:a:i:o:C:S:f:l:L:T:M:Rtvh"
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
			annotation+=($OPTARG)
			;;
		t)
			no_trim=true
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

reconstruct_fasta=$group/$sample/mapping/$sample".coverage_adapted_filtered_"$coverage_cutoff"_term.fasta"_$cluster_cutoff
#contigs=$group/$sample/assembly/scaffolds.fasta



lib/check_mandatory_files.sh $r1_file $r2_file $database 

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


if [ $explore = true ]; then
	coverage_cutoff=60
	coverage_summary=70
	cluster_cutoff=70
	alignment_percentage=20
fi

if [ $only_reconstruct = false ]; then

	reconstruct_fasta=$group/$sample/mapping/$sample".coverage_adapted_filtered_"$coverage_cutoff"_term.fasta"_$cluster_cutoff

####TRIMMING#################################################################
#############################################################################
	
	if [ $no_trim = false ]; then

		r1_file_mapping=$group/$sample/trimmed/$sample"_1_paired.fastq.gz"
		r2_file_mapping=$group/$sample/trimmed/$sample"_2_paired.fastq.gz"



		filestrimmed=$(ls $group/$sample/trimmed 2> /dev/null | grep "paired" | wc -l)

		if [ "$filestrimmed" = "4" ];then
			echo "Found trimmed files for this sample" $sample;
			echo "Omitting trimming"
		else
			echo "####TRIMMING################################################################"
			lib/quality_trim.sh -1 $r1_file -2 $r2_file -s $sample -g $group -T $threads
		fi
	else

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
		echo "####ASSEMBLY################################################################"
		lib/spades_assembly.sh -q $group/$sample/trimmed/ -c -T $threads
		contigs=$group/$sample/assembly/scaffolds.fasta
	else
		echo "Contigs supplied, ommiting assembly step"
	fi

	#group/sample/assembly/scaffolds.fasta

####MAPPING##################################################################
#############################################################################
	if [ -f $group/$sample/mapping/$sample.sorted.bam -a -f $group/$sample/mapping/$sample.sorted.bam.bai -o -f $group/$sample/mapping/$sample.sam ];then
		echo "Found a mapping file for sample" $sample;
		echo "Omitting mapping"
	else
		echo "####MAPPING#################################################################"
		lib/bowtie_mapper.sh -d $database \
	 	-g $group \
	 	-s $sample \
	 	-1 $r1_file_mapping \
	 	-2 $r2_file_mapping

 #group/sample/mapping/sample.sam

	 	lib/sam_to_bam.sh -i $group/$sample/mapping/$sample.sam

 #group/sample/mapping/sample.bam
 #group/sample/mapping/sample.bam.bai
 	fi

####COVERAGE FILTERING#######################################################
#############################################################################
	if [ -f $group/$sample/mapping/$sample".coverage" ];then \
		echo "Found a coverage file for sample" $sample;
		echo "Omitting coverage calculation"
	else
		echo "####COVERAGE FILTERING########################################################"
 		lib/get_coverage.sh -i $group/$sample/mapping/$sample".sorted.bam" -d $database
 	fi
 #group/sample/mapping/sample.coverage

 	lib/adapt_filter_coverage.sh -i $group/$sample/mapping/$sample".coverage" -c $coverage_cutoff

 #sample.coverage_adapted
 #sample.coverage_adapted_filtered_80

	lib/filter_fasta.sh -i $database -f $group/$sample/mapping/$sample".coverage_adapted_filtered_"$coverage_cutoff

#sample.coverage_adapted_filtered_50_term.fasta
	if [ -f $group/$sample/mapping/$sample".coverage_adapted_filtered_"$coverage_cutoff"_term.fasta_"$cluster_cutoff ];then \
		echo "Found a clustered file for sample" $sample;
		echo "Omitting clustering"
	else

		lib/cdhit_cluster.sh -i $group/$sample/mapping/$sample".coverage_adapted_filtered_"$coverage_cutoff"_term.fasta" -c $cluster_cutoff -M $max_memory -T $threads
	fi
#$group/$sample/mapping/$sample".coverage_adapted_filtered_"$coverage_cutoff"_term.fasta"_$cluster_cutoff >> FINAL CLUSTERED AND FILTERED FASTA FILE TO USE AS SCAFFOLD
########################################################################################################

#sample.coverage_adapted_filtered_50_term.fasta_80
#sample.coverage_adapted_filtered_50_term.fasta_80.clstr

	lib/process_cluster_output.sh -i $reconstruct_fasta -b $group/$sample/mapping/$sample".coverage_adapted" -c $coverage_cutoff

fi


##sample.coverage_adapted_clustered
##sample.coverage_adapted_clustered_percentage
##sample.coverage_adapted_clustered_ac
############################################################

###############################################################################################################################################
######################## DONE WITH MAPPING AND CLUSTERING #####################################################################################
###############################################################################################################################################

echo "####CONTIG and ANNOTATON########################################################"

lib/build_karyotype.sh -i $group/$sample/mapping/$sample".coverage_adapted_clustered" -K $coverage_summary -k $coverage_cutoff -o $group/$sample/data

#group/sample/data
#sample.karyotype_individual.txt
#sample.karyotype_individual.txt

if [ -f $group/$sample/data/$sample".bedgraph" ];then
	echo "Found a coverage file for sample" $sample;
	echo "Omitting coverage calculation"
else
	lib/get_coverage.sh -i $group/$sample/mapping/$sample".sorted.bam" -p -o $group/$sample/data
fi
#sample.bedgraph

lib/filter_fasta.sh -i $group/$sample/data/$sample".bedgraph" -f $group/$sample/mapping/$sample".coverage_adapted_clustered_ac" -G

#sample.bedgraph_term

if [ -f $group/$sample/data/$sample".fna" -a -f $group/$sample/data/$sample".gff" ];then
	echo "Found annotations files for sample" $sample;
	echo "Omitting automatic annotation"
else
	lib/prokka_annotation.sh -i $contigs -p $sample -o $group/$sample/data -c
fi
#sample.fna
#sample.gff


lib/blast_align.sh -i $group/$sample/data/$sample".fna" -d $reconstruct_fasta -o $group/$sample/data -p plasmids

#sample.plasmids.blast


lib/blast_to_bed.sh -i $group/$sample/data/$sample".plasmids.blast" -b $alignment_identity -l 0 -L 500 -d - -q _ -Q r -I

#sample.plasmids.bed

lib/blast_to_complete.sh -i $group/$sample/data/$sample".plasmids.blast"

#sample.plasmids.complete

lib/blast_to_link.sh -i $group/$sample/data/$sample".plasmids.blast" -I

#sample.plasmids.links
#sample.plasmids.blast.links

lib/gff_to_bed.sh -i $group/$sample/data/$sample".gff" -L

#sample.gff.bed

lib/coordinate_adapter.sh -i $group/$sample/data/$sample".gff.bed" -l $group/$sample/data/$sample".plasmids.blast.links" -p -n 5000

#sample.gff.coordinates


echo "####SPECIFIC ANNOTATON########################################################"


######################### ABR _ INCLUDE FILENAME

lib/blast_align.sh -i /processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PSMP_T/ANALYSIS/PLASMIDID/references/ARGannot.r1.pID.fasta -d $group/$sample/data/$sample".fna" -o $group/$sample/data -p abr -f $sample

#sample.abr.blast

lib/blast_to_bed.sh -i $group/$sample/data/$sample".abr.blast" -b 100 -l 90 -d _ -D r -q " " -Q r

#sample.abr.bed

lib/coordinate_adapter.sh -i $group/$sample/data/$sample".abr.bed" -l $group/$sample/data/$sample".plasmids.blast.links" -u

#sample.abr.coordinates

###################### INC

lib/blast_align.sh -i /processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PSMP_T/ANALYSIS/PLASMIDID/references/plasmidFinder_01_26_2018.fsa -d $group/$sample/data/$sample".fna" -o $group/$sample/data -p inc -f $sample

#sample.inc.blast

lib/blast_to_bed.sh -i $group/$sample/data/$sample".inc.blast" -b 95 -l 80 -d _ -D r -q _ -Q l
#sample.inc.bed

lib/coordinate_adapter.sh -i $group/$sample/data/$sample".inc.bed" -l $group/$sample/data/$sample".plasmids.blast.links" -u

#sample.inc.coordinates


lib/draw_circos_images.sh $group $sample