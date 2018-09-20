#!/bin/bash

set -e

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.0
#CREATED: 6 April 2018
#REVISION:
#DESCRIPTION:Script that uses cd-hit/psi-cd-hit to clusterize a FASTA file
#
#DOCUMENTATION
#
#
#Compare floats in BASH
#
#if [ $(echo "$cluster_cutoff > 0.7"|bc -l) -eq 1 ]; then
#	echo "YES"
#else
#	echo "NO"
#fi
#
#-d length of description in .clstr file, default 20. if set to 0,
#	it takes the fasta defline and stops at first space
#-s length difference cutoff, default 0.0
#	if set to 0.9, the shorter sequences need to be
#	at least 90% length of the representative of the cluster
#-B 1 or 0, default 0, by default, sequences are stored in RAM
#	if set to 1, sequence are stored on hard drive
#	it is recommended to use -B 1 for huge databases
#-g 1 or 0, default 0
#	By cd-hit’s default algorithm, a sequence is clustered to the first
#	cluster that meet the threshold (fast mode). If set to 1, the program
#	will cluster it into the most similar cluster that meet the threshold
#	(accurate but slow mode)
#
#	PSI-CD-HIT
#-G (1/0) use global identity? default 1, sequence identity
#	calculated as total identical residues of local alignments
#	length of shorter sequence
#
#-n 5 for thresholds 0.7 ~ 1.0
#-n 4 for thresholds 0.6 ~ 0.7
#-n 3 for thresholds 0.5 ~ 0.6
#-n 2 for thresholds 0.4 ~ 0.5

#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

Cdhit_cluster script uses cd-hit/psi-cd-hit to clusterize a FASTA file

usage : $0 <-i inputfile(FASTA)> [-o <directory>] [-n <filename>] [-c <percentage>]
		[-T <threads>] [-g group_name] [-s <int>] [-M <int>][-C <(0|1)>] [-G <(0|1)>] [-b <blast_prog>] [p] [-v] [-h]

	-i input file in FASTA format
	-c percentage threshold to cluster, default 80
	-M max available memory (Mbyte), default 400
	-n file name
	-s length difference cutoff, default 0.8
	-g group name (optional). If unset, samples will be gathered in NO_GROUP group
	-p runs psi-cd-hit instead of cd-hit
	-C psi-cd-hit only: circular sequences, default 1. If set to 0 sequence is assumed lineal
	-G psi-cd-hit only: gobal identity, -G 0 only takes the first local alignment for clustering
	-b psi-cd-hit only: choose blast program, default blastn
	-T number of threads
	-v version
	-h display usage message


Output directory is the same as input directory

example: cdhit_cluster -i ecoli.fasta -c 90 -M 50000 -T 0


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
cluster_cutoff=0.8
max_memory=400
length_cutoff=0.8
cd_hit_command=cd-hit-est
is_circle=1
global_psi_cd_hit=1
psi_cd_hit_program=blastn
word_size=0
threads=0

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:o:c:M:n:s:g:C:G:b:T:pvh"
while getopts $options opt; do
	case $opt in
		i )
			input_file=$OPTARG
			;;

		c )
			cluster_cutoff_input=$OPTARG
			;;
		g)
			group=$OPTARG
			;;
		M )
			max_memory=$OPTARG
			;;
		n )
			file_name=$OPTARG
			;;
		s )
          	length_cutoff=$OPTARG
          	;;
        p )
          	cd_hit_command=psi-cd-hit.pl
          	;;
        C )
          	is_circle=$OPTARG
          	;;
        G)
          	global_psi_cd_hit=$OPTARG
          	;;
        T)
          	threads=$OPTARG
          	;;
        b)
          	psi_cd_hit_program=$OPTARG
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

check_dependencies.sh cd-hit-est
 #psi-cd-hit.pl



# Set word size (parameter -n for cd-hit) as author recomends
#according to clustering percentage


cluster_cutoff=$(echo "$cluster_cutoff_input / 100" | bc -l | sed 's/0\{1,\}$//')
#cluster_cutoff=${cluster_cutoff%.*} #Remove float value


if [[ "$cluster_cutoff_input" -gt 70  &&  "$cluster_cutoff_input" -le 100 ]]; then
	word_size=5
elif [[ "$cluster_cutoff_input" -gt 60  &&  "$cluster_cutoff_input" -le 70 ]]; then
	word_size=4
elif [[ "$cluster_cutoff_input" -gt 50  &&  "$cluster_cutoff_input" -le 60 ]]; then
	word_size=3
elif [[ "$cluster_cutoff_input" -ge 40  &&  "$cluster_cutoff_input" -le 50 ]]; then
	word_size=2
else
	echo "please introduce a valid cluster percentage value between 0.4 and 1"
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

##CD-HIT EXECUTION

echo "$(date)"
echo "Clustering sequences with identity" $cluster_cutoff_input"% or higher"
echo "Using" $cd_hit_command "with file" $input_file
seq_number_prev_clstr=$(cat $input_file | grep ">" | wc -l)

cd $(dirname $input_file)

if [ -f $output_dir/$file_name""_""$cluster_cutoff_input ]; then \
	echo "Found a clustered file for sample" $file_name;
	echo "Omitting clustering process calculation"
	exit 1
else
	if [ $cd_hit_command  == "psi-cd-hit.pl" ]; then

		check_dependencies.sh psi-cd-hit.pl
		$cd_hit_command -i $(basename $input_file) -o $file_name""_""$cluster_cutoff_input -c $cluster_cutoff -G $global_psi_cd_hit -g 1 -prog $psi_cd_hit_program -circle $is_circle -core $threads || error ${LINENO} $(basename $0) "PSI-CD-HIT command failed. See $output_dir/logs for more information."

	else

		$cd_hit_command -i $(basename $input_file) -o $file_name""_""$cluster_cutoff_input -c $cluster_cutoff -n $word_size -d 0 -s $length_cutoff -B 1 -M $max_memory -T $threads|| error ${LINENO} $(basename $0) "CD-HIT command failed. See $output_dir/logs for more information"

	fi
fi

seq_number_post_clstr=$(cat $file_name""_""$cluster_cutoff_input | grep ">" | wc -l)

echo "$(date)"
echo "DONE Clustering sequences with identity" $cluster_cutoff_input"% or higher"
echo "fasta file can be found in" $output_dir/$file_name""_""$cluster_cutoff_input
echo "Previous number of sequences=" $seq_number_prev_clstr
echo -e "Number of sequences after clustering=" $seq_number_post_clstr "\n"
cd $cwd
