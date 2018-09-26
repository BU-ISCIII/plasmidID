#!/bin/bash

# Exit immediately if a pipeline, which may consist of a single simple command, a list,
#or a compound command returns a non-zero status: If errors are not handled by user
#set -e
#set -x

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.0
#CREATED: 13 May 2018
#
#DESCRIPTION:blast_to_complete script obtain full length of sequences from blast and adapt it to circos
#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

blast_to_complete is a script that obtain  full length of sequences from blast and adapt it to circos

usage : $0 <-i inputfile(.blast)> <-b id cutoff> [-o <directory>] [-b <int(0-100)>] [-l <int(0-100)>]
		[-p <prefix>] [-d <delimiter>] [-D (l|r)] [-q <delimiter>] [-Q (l|r)] [-I] [-u] [-v] [-h]

	-i input file
	-b blast identity cutoff (0 - 100), default 90
	-l blast length percentage cutoff (0 - 100), default 50, use 90 for genes
	-o output directory (optional). By default the file is replaced in the same location
	-q database chraracter delimiter, default "_"
	-Q query field to retrieve (l=left, r=right), default left
	-d database chraracter delimiter, default "_"
	-D database field to retrieve (l=left, r=right), default right
	-I contig mode
	-u unique. Outputs only one query entry per database entry
	-v version
	-h display usage message

example: blast_to_complete.sh -i ecoli_prefix.blast
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
input_file="Input_file"
blast_id_cutoff=90
blast_len_percentage=15
database_delimiter="-"
database_field=r
query_delimiter="_"
query_field=r
unique=false
suffix=""
id_circos=false
id_output=""

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:b:q:Q:d:D:o:l:Iuvh"
while getopts $options opt; do
	case $opt in
		i )
			input_file=$OPTARG
			;;
		b )
			if [ $OPTARG -lt 0 ] || [ $OPTARG -gt 100 ]; then
				echo "please, provide a percentage between 0 and 100"
				exit 1
			else
				blast_id_cutoff=$OPTARG
			fi
			;;
		o )
			output_dir=$OPTARG
			;;
		l )
			if [ $OPTARG -lt 0 ] || [ $OPTARG -gt 100 ]; then
				echo "please, provide a percentage between 0 and 100"
				exit 1
			else
				blast_len_percentage=$OPTARG
			fi
			;;
		d )
			database_delimiter=$OPTARG
			;;
		D )
			database_field=$OPTARG
			;;
        q )
			query_delimiter=$OPTARG
			;;
		Q )
			query_field=$OPTARG
			;;
        u )
			unique=true
            suffix=".unique.tmp"
			;;
        I)
			id_circos=true
            id_output=",\"id=\"database_name[length(database_name)]"
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


blast_len_percentage_value=$(echo "($blast_len_percentage/100)" | bc -l)
#blast_len_percentage_decimal=$(echo $blast_len_percentage_value | sed 's/0\{1,\}$//')


if [ ! $output_dir ]; then
	output_dir=$(dirname $input_file)
	#echo "Default output directory is" $output_dir
	mkdir -p $output_dir
else
	#echo "Output directory is" $output_dir
	mkdir -p $output_dir
fi


if [ ! $file_name ]; then
	file_name=$(basename $input_file | cut -d. -f1,2)
fi

##CHECK FIELDS TO RETRIEVE

if [ "$database_field" == "l" ] || [ "$database_field" == "r" ]; then

	if [ $database_field == l ]; then
		database_field="1"
	else
		database_field="length(database_name)"
	fi

else
	echo "Please introduce r or l for database"
	exit 1
fi

if [ $query_field == "l" ] || [ $query_field == "r" ]; then

	if [ $query_field == l ]; then
		query_field="1"
	else
		query_field="length(query_name)"
	fi

else

	echo "Please introduce 0 or 1 for query"
	exit 1
fi

echo "$(date)"
echo "Adapting blast to complete using" $(basename $input_file) "with:"
echo "Blast identity=" $blast_id_cutoff
echo "Min len percentage=" $blast_len_percentage


cat $input_file |\
awk '
	BEGIN{OFS="\t"}
	{split($1, query_name, "'"${query_delimiter}"'")
	split($2,database_name, "'"${database_delimiter}"'")}
	(($3 >= '"${blast_id_cutoff}"') && (($4/$13)>='"${blast_len_percentage_value}"') && (!x[$1$2]++)) \
	{{isInverted=($10-$9)
		ext2=($13-$8)}
		{if (isInverted < 0)
			{pos1 = $10
			pos2 = $9}
		else
			{pos1 =$9
			pos2 = $10}
		{if ((isInverted < 0) && (($14 - pos2) > $7))
			{coordChr2 = (pos2 + $7)}
		else if ((isInverted < 0) && (($14 - pos2) <= $7))
			{coordChr2=$14}
		{if ((isInverted < 0) && (ext2 <= pos1))
			{coordChr1= pos1 - ext2;}
		else if ((isInverted < 0) && (ext2 > pos1))
			{coordChr1= 1}
		{if ((isInverted > 0) && (pos1 > $7))
			{coordChr1=(pos1 - $7)}
		else if ((isInverted > 0) && (pos1 <= $7))
			{coordChr1=1}
		{if ((isInverted > 0) && (ext2 > ($14-pos2)))
			{coordChr2= $14;}
		else if ((isInverted > 0) && (ext2 <= ($14-pos2)))
			{coordChr2= (pos2 + ext2)}
	{print database_name['"$database_field"'], coordChr1, coordChr2, query_name['"$query_field"'], "id="$13} }}}}}}
	' \
	>$output_dir/$file_name".complete"|| error ${LINENO} $(basename $0) "Awk command parsing blast output for circos input in $file_name\".complete\" creation failed. See $output_dir/logs for more information"


cat $input_file |\
awk '
	BEGIN{OFS="\t"}
	{split($1, query_name, "'"${query_delimiter}"'")
	split($2,database_name, "'"${database_delimiter}"'")}
	(($3 >= '"${blast_id_cutoff}"') && (($4/$13)>='"${blast_len_percentage_value}"') && (!x[$2$1]++)) \
	{{isInverted=($10-$9)
	ext2=($13-$8)}
	{if (isInverted < 0)
		{pos1=$10
		pos2=$9}
	else
		{pos1 =$9
		pos2=$10}; \
	{if ((isInverted < 0) && (($14 - pos2) < $7))
		{coordChr1=1
		coordChr2=($7-($14-pos2))
		{print database_name['"$database_field"'], coordChr1, coordChr2, query_name['"$query_field"'], "id="$13}}
	{if ((isInverted < 0) && (ext2 > pos1))
		{coordChr1=($14-(ext2-pos1))
		coordChr2=$14
		{print database_name['"$database_field"'], coordChr1, coordChr2, query_name['"$query_field"'], "id="$13}}
	{if ((isInverted > 0) && (pos1 < $7))
		{coordChr1=($14-($7-pos1))
		coordChr2=$14
		{print database_name['"$database_field"'], coordChr1, coordChr2, query_name['"$query_field"'], "id="$13}}
	{if ((isInverted > 0) && (ext2 > ($14-pos2)))
		{coordChr1=1
		coordChr2=(ext2-($14-pos2))
		{print database_name['"$database_field"'], coordChr1, coordChr2, query_name['"$query_field"'], "id="$13}
		}
	}}}}}}
	' \
	>>$output_dir/$file_name".complete" || error ${LINENO} $(basename $0) "Awk command parsing blast output for circos input in $file_name\".complete\" second step creation failed. See $output_dir/logs for more information"



echo "$(date)"
echo "DONE adapting blast to complete"
echo -e "File can be found at" $output_dir/$file_name".complete" "/n"
