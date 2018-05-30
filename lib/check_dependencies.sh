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
#DESCRIPTION:Short function to evaluate if programs are on path

#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

Check_dependencies Short function to evaluate if files exist

usage : $0 <program_name1> [program_name2] ...

example: lib/check_dependencies.sh foo bar

EOF
}

if [ $# = 0 ] ; then
 usage >&2
 exit 1
fi

#DECLARE FLAGS AND VARIABLES
missing_dependencies=0

for command in "$@"; do
	if ! [ -x "$(which $command 2> /dev/null)" ]; then
		echo "Error: Please install $command or make sure it is in your path" >&2
		let missing_dependencies++
	else
		echo "$command installed"
	fi
done

if [ $missing_dependencies -gt 0 ]; then 
	echo "ERROR: $missing_dependencies missing dependencies, aborting execution" >&2
	exit 1
fi

