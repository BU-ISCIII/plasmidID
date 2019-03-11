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
#AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
VERSION=1.3.3
#CREATED: 15 March 2018
#
#ACKNOLEDGE: longops2getops.sh: https://gist.github.com/adamhotep/895cebf290e95e613c006afbffef09d7
#
#DESCRIPTION: test.sh uses test data for testing plasmidID installation.
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

usage : $0

	-v | --version		version
	-h | --help		display usage message

example: ./test.sh

EOF
}

#================================================================
# OPTION_PROCESSING
#================================================================
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

# translate long options to short
reset=true
for arg in "$@"
do
    if [ -n "$reset" ]; then
      unset reset
      set --      # this resets the "$@" array so we can rebuild it
    fi
    case "$arg" in
       	--help)    	set -- "$@" -h ;;
       	--version) 	set -- "$@" -v ;;
       # pass through anything else
       *)         set -- "$@" "$arg" ;;
    esac
done

#DECLARE FLAGS AND VARIABLES
script_dir=$(dirname $(readlink -f $0))
R1=KPN_TEST_R1.fastq.gz
R2=KPN_TEST_R2.fastq.gz
database=plasmids_TEST_database.fasta
contigs=contigs_KPN_TEST.fasta

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":1:2:d:s:g:c:a:i:o:C:S:f:l:L:T:M:X:y:Y:RVtvh"
while getopts $options opt; do
	case $opt in
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

## Execute plasmidID with test data.
echo "Executing:../plasmidID.sh -1 $R1 -2 $R2 -d $database -c $contigs -s KPN --no-trim"
echo "Forward reads: $R1"
echo "Reverse reads: $R2"
echo "PlasmidDatabase: $database"
echo "Contigs: $contigs"
echo "Options: --no-trim"

export PATH=$PATH:$script_dir/bin
$script_dir/../plasmidID.sh -1 $script_dir/$R1 -2 $script_dir/$R2 -d $script_dir/$database -c $script_dir/$contigs -s KPN --no-trim


echo "ALL DONE. TEST COMPLETED SUCCESSFULLY YOUR INSTALLATION SHOULD BE CORRECT."
