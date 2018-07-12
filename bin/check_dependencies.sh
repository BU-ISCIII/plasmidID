#!/bin/bash

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.0 
#CREATED: 19 March 2018
#REVISION: 12 July 2018: add formated output and colors
#AKNOWLEDGE: Colored text: https://stackoverflow.com/questions/5947742/how-to-change-the-output-color-of-echo-in-linux
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

#SET COLORS

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

printf '\n%s\t%30s\n' "DEPENDENCY" "STATUS"
printf '%s\t%30s\n'   "----------" "------"

for command in "$@"; do
	#dependency_version=$($command --version)
	length_command=$(echo $command | wc -m)
	distance_table=$((30 - $length_command))
	distance_expression=$(echo "%${distance_table}s")
	
	printf '%s' $command
	if ! [ -x "$(which $command 2> /dev/null)" ]; then
		
		
		printf $distance_expression
		printf "${RED}NOT INSTALLED${NC} \n"
		let missing_dependencies++
	else
		printf $distance_expression
		printf "${GREEN}INSTALLED${NC} \n"
	fi
done

if [ $missing_dependencies -gt 0 ]; then 
	printf "${RED}ERROR${NC}: $missing_dependencies missing dependencies, aborting execution\n" >&2
	exit 1
fi

