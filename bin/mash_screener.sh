#!/bin/bash

set -e

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.0
#CREATED: 27 November 2019
#REVISION:

#DESCRIPTION:Script that screen reads over a database using kmers and estract sequences ids with higher values
#TODO
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
	-d database to screen (.fasta)
	-s sample name
	-g group name (optional). If unset, samples will be gathered in NO_GROUP group
	-1 reads corresponding to paired-end R1
	-2 reads corresponding to paired-end R2
	-f threshold identity value to retieve sequence ids with at least this value (default 0.9)
	-w use winner takes it all
	-T number of threads
	-v version
	-h display usage message

example: mash_screener.sh -d database.fasta -s COLI -1 ecoli_1.fastq -2 ecoli_2.fastq

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
threads=1
offrate=1
filter_identity=0.9
cwd="$(pwd)"
w_winner=""
group="NO_GROUP"
database="Database"
R1="R1"
R2="R2"

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:o:s:g:d:1:2:f:T:avwh"
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
          	filter_identity=$OPTARG
      		;;
		w)
			w_winner="-w"
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

check_dependencies.sh bash mash

check_mandatory_files.sh $database $R1

if [ ! $sample ]; then
	echo "ERROR: please, provide a sample name"
	usage
	exit 1
fi

if [ ! $output_dir ]; then
	output_dir=$cwd"/$group/$sample/kmer/"
	echo "Default output directory is" $output_dir
	mkdir -p $output_dir
else
	echo "Output directory is" $output_dir
	mkdir -p $output_dir
fi


########SKETCH##############
############################

if [ -f $output_dir/database.msh ]; then \
	echo "Found a sketch ddbb for" $(basename $database);
	echo "Omitting sketching"
else
	echo "creating sketch of " $(basename $database);
	mash sketch -i -k 32 -s 1000 -p $threads -o $output_dir/database $database || error ${LINENO} $(basename $0) "mash screen command failed. See $output_dir/logs for more information"
fi

########SCREEN##############
############################

if [ -f $output_dir/database.screen.tab ];then \
	echo "Found a mash screen file for sample" $sample;
	echo "Omitting screening"
else
	echo "$(date)"
	echo screening $R1

	mash screen $w_winner -p $threads $output_dir/database.msh $R1 > $output_dir/database.screen.tab || error ${LINENO} $(basename $0) "Bowtie2 command failed. See $output_dir/logs for more information"


	echo "$(date)"
	echo -e "DONE Screening $sample of $group Group" "\n"
fi

######PARSE_RESULT##########
############################

output_mash_id=$output_dir/database.filtered_$filter_identity

echo "Retrieving sequences matching more than $filter_identity identity"

cat $output_dir/database.screen.tab | awk '($1 >= '"${filter_identity}"') {print $5}' > $output_mash_id


#####FILTER SEQUENCES#######
############################
if [ $(cat $output_mash_id | wc -l | cut -d " " -f 1) -gt 0 ]
then
	filter_fasta.sh -i $database -f $output_mash_id
else
	echo "No plasmids have passed the mash identity filter!! Exiting!!"
	exit 0
fi
