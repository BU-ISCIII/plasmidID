#!/bin/bash

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.0 
#CREATED: 19 March 2018
#REVISION:
#DESCRIPTION:Short function to evaluate if files exist

#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

Check_mandatory_files Short function to evaluate if files exist

usage : $0 <file1> [file2] ...

example: lib/check_mandatory_files.sh foo.txt bar.fasta

EOF
}

if [ $# = 0 ] ; then
 usage >&2
 exit
fi

#DECLARE FLAGS AND VARIABLES
missing_files=0

for file in "$@"; do
	if [ ! -f $file ]; then
		echo "$(basename $file)" "not supplied, please, introduce a valid file" >&2
		let missing_files++
	fi
done

if [ $missing_files -gt 0 ]; then 
	echo "ERROR: $missing_files missing files, aborting execution" >&2
	exit 1
fi