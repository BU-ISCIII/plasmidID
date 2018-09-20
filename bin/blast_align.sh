#!/bin/bash

set -e

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.0
#CREATED: 1 May 2018
#REVISION:
#DESCRIPTION:Script that blast a query against a database
#
#DOCUMENTATION
#
#Blast output6 with aditions:
#1	 	Query label.(qseqid)
#2	 	Target or subject(database sequence or cluster centroid) label. (sseqid)
#3	 	Percent identity. (pident)
#4	 	Alignment length. (length)
#5	 	Number of mismatches. (mismatch)
#6	 	Number of gap opens. (gapopen)
#7	 	Start position in query. Query coordinates start with 1 at the first base in the sequence as it appears in the input file. For translated searches (nucleotide queries, protein targets), query start<end for +ve frame and start>end for -ve frame. (qstart)
#8	 	End position in query. (qend)
#9	 	Start position in target. Target coordinates start with 1 at the first base in sequence as it appears in the database. For untranslated nucleotide searches, target start<end for plus strand, start>end for a reverse-complement alignment. (sstart)
#10	 	End position in target. (send)
#11	 	E-value calculated using Karlin-Altschul statistics. (evalue)
#12	 	Bit score calculated using Karlin-Altschul statistics. (bitscore)
#13		Lenght of query (qlen)
#14		Length of target (slen)
#
#
#TO DO:
#
#Handle all types of blast: blastn, blastp...
#
#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

blast_align is a script that blast a query against a database

usage : $0 <-i inputfile(query)> <-d inputfile(database)> [-p <prefix>] [-o <directory>] [-t <nucl|prot>]
		[-T <threads>] [-e <evalue>] [-v] [-h]

	-i query file in FASTA format
    -d database to blast against
	-o output directory, default same directory as query
	-p prefix for blast identification (mandatory) and output file name
	-t type of database, nucl by default
    -e evalue for blast analysis, default 0.0001
	-T number of threads
	-v version
	-h display usage message

Output directory is the same as input directory by default

example: blast_align -i ecoli.fasta -d plasmid_ddbb.fasta -p plasmid


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
database="Database"
database_type="nucl"
evalue=0.0001
threads=1

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:o:p:f:d:t:e:T:vh"
while getopts $options opt; do
	case $opt in
		i )
			input_file=$OPTARG
			;;
        d )
			database=$OPTARG
			;;
		o )
			output_dir=$OPTARG
			;;
		p)
			prefix=$OPTARG
			;;
		f)
			file_name=$OPTARG
			;;
		t )
          	database_type=$OPTARG
          	;;
        g )
          	group=$OPTARG
          	;;
		e )
          	evalue=$OPTARG
          	;;
        T)
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

check_mandatory_files.sh $input_file $database

#check_dependencies.sh blastn


if [ ! $prefix ]; then
	echo "please provide a prefix to identify this blast analysis"
	exit 1
fi

if [ $database_type == "prot" ] || [ $database_type == "nucl" ]; then
	echo "database type selected as" $database_type
else
	echo "please provide a proper database type"
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
	file_name=$(basename $input_file | cut -d. -f1)
	echo "filename is" $file_name
fi

database_name=$(basename $database)
database_dir=$(dirname $database)

##BLAST EXECUTION

echo "$(date)"
echo "Blasting" $file_name "agaist" $database_name

makeblastdb -in $database -out $database_dir/$database_name".blast.tmp" -dbtype $database_type || error ${LINENO} $(basename $0) "Makeblastdb command failed. See $output_dir/logs for more information."


blastn -query $input_file \
-db $database_dir/$database_name".blast.tmp" \
-out $output_dir/$file_name"."$prefix".blast" \
-evalue $evalue \
-num_threads $threads \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" || error ${LINENO} $(basename $0) "Blastn command failed. See $output_dir/logs for more information"


echo "$(date)"
echo "Done blasting" $file_name "agaist" $database_name
echo -e "blasted file can be found in" $output_dir/$file_name"."$prefix".blast" "\n"
