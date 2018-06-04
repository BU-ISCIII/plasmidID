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
#CREATED: 4 May 2018
#REVISION: 
#06 May 2018: add id optiopn in bed output
#04 June 2018: add an option for an aditional division mostly for ABR sort
#
#DESCRIPTION:blast_to_bed script obtain a BED file with coordinates of local blast alignments matching some given conditions
#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

blast_to_bed is a script than obtain a BED file with coordinates of local blast alignments matching some given conditions

usage : $0 <-i inputfile(.blast)> <-b id cutoff> [-o <directory>] [-b <int(0-100)>] [-l <int(0-100)>] [-L <int>]
		[-p <prefix>] [-d <delimiter>] [-D (l|r)] [-q <delimiter>] [-Q (l|r)] [-U <delimiter>] [-I] [-u] [-v] [-h]

	-i input file 
	-b blast identity cutoff (0 - 100), default 90
	-l blast length percentage cutoff (0 - 100), default 20, use 90 for genes
	-L blast length alignment cutoff, default 0, use 200 or 500 for contigs
	-o output directory (optional). By default the file is replaced in the same location
	-q database chraracter delimiter, default "_"
	-Q query field to retrieve (l=left, r=right), default left
	-d database chraracter delimiter, default "_"
	-D database field to retrieve (l=left, r=right), default right
	-I contig mode
	-u unique. Outputs only one query entry per database entry
	-U unique mode with delimiter. Outputs only one delimited query per database entry
	-v version
	-h display usage message

example: blast_to_bed.sh -i ecoli_prefix.blast -b 80 -l 50 -q - -Q r

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


#DECLARE FLAGS AND VARIABLES
cwd="$(pwd)"
input_file="Input_file"
blast_id_cutoff=90
blast_len_percentage=10
blast_len_alignment=0
database_delimiter="_"
database_field=r
query_delimiter="_"
query_field=l
unique=false
unique_divider=false
divider_delimiter="-"
suffix=""
id_circos=false
id_output=""

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:b:q:Q:d:D:o:l:L:U:Iuvh"
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
		L )
			blast_len_alignment=$OPTARG
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
		U )
			unique_divider=true
			suffix=".unique.divider.tmp"
			divider_delimiter=$OPTARG
			;;
        I)
			id_circos=true
            id_output=",\"id=\"query_name[length(query_name)]"
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

bash lib/check_mandatory_files.sh $input_file


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
echo "Adapting blast to bed using" $(basename $input_file) "with:"
echo "Blast identity=" $blast_id_cutoff
echo "Min length aligned=" $blast_len_alignment
echo "Min len percentage=" $blast_len_percentage


cat $input_file | sort -k3 -nr | \
awk '
	{OFS="\t"
	split($2, database_name, "'"${database_delimiter}"'")
	split($1, query_name, "'"${query_delimiter}"'")}
	(($3 >= '"${blast_id_cutoff}"')&&(($4/$13) >= '"${blast_len_percentage_value}"')&&($4 >= '"${blast_len_alignment}"')) \
	{print database_name['"$database_field"'], $9, $10, query_name['"$query_field"']'"$id_output"'}
	' \
> $output_dir/$file_name".bed"$suffix


if [ "$unique" == "true" ]; then
	echo "unique option enabled"
    awk '
        (!x[$1$4]++)
    	' $output_dir/$file_name".bed"$suffix \
> $output_dir/$file_name".bed"
fi


if [ "$unique_divider" == "true" ]; then
	echo "unique delimiter option enabled"
    awk '
    	{split($4,query,"'"${divider_delimiter}"'")}
        (!x[query[1]$1]++)
    	' $output_dir/$file_name".bed"$suffix \
	> $output_dir/$file_name".bed"
fi

rm $output_dir/$file_name".bed"$suffix

echo "$(date)"
echo "DONE adapting blast to bed"
echo -e "File can be found at" $output_dir/$file_name".bed" "\n"