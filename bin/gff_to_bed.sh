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
#CREATED: 16 May 2018
#REVISION:
#DESCRIPTION:gff_to_bed script obtain a list of genes with name from a gff. [Tested with prokka output]
#
#TO DO:
#Test in a regular ncbi gff
#
#
#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

gff_to_bed script obtain a list of genes with name from a GFF file. [Tested with prokka output]

usage : $0 <-i inputfile(.fasta)> [-o <directory>] [-C] [-L] [-q <character>] [-Q (l|r)] [-s <suffix>] [-u] [-v] [-h]

	-i input file 
	-o output directory (optional). By default the file is placed in the same location as input
	-C include a supplied word in cds
	-L include locus tag in cds
	-q database chraracter delimiter, default "_"
	-Q query field to retrieve (l=left, r=right), default right
	-u uniq mode. Remove duplicates
	-s string to ad at the end of the outputted file
	-v version
	-h display usage message

example: ./gff_to_bed.sh -i ecoli.gff -C CDS_ -L

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
query_delimiter="_"
query_field=r
unique=false
suffix=""
cds_word=""
cds_locus=""
#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:b:o:q:Q:C:Luvh"
while getopts $options opt; do
	case $opt in
		i )
			input_file=$OPTARG
			;;
		b )
			coverage_file=$OPTARG
			;;
		o )
			output_dir=$OPTARG
			;;
		C )
			cds_word=\"$OPTARG\"
			;;
		L )
			cds_locus="locusname[length(locusname)]"
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
echo "Getting bed file from GFF in" $file_name

#Filter Gff(3) file from prokka and create a bed file with coordinates with annotated genes (WITH NAMES)
awk '
	BEGIN{OFS="\t"}
	{split($1, query_name, "'"${query_delimiter}"'")
	split($9,description,"Name=")
	split(description[2],name,";")
	split($9,locustag,"locus_tag=")
	split(locustag[2],locus,";")
	split(locus[1],locusname,"_")
	split(name[1],nameLowBar,"_")}
	{if ($3 == "gene" && $1 != "#" && $9 ~ /Name/)
		{print query_name['"$query_field"'],$4,$5, name[1]}
	{if ($3 == "gene" && $1 != "#" && $9 !~ /Name/)
		{print query_name['"$query_field"'],$4,$5,'"${cds_word}"' "" '"${cds_locus}"'}
	}
	}
	' \
	$input_file \
	> $output_dir/$file_name".bed"$suffix


if [ "$unique" == "true" ]; then
    awk '
        (!x[$1$4]++)
    	' $output_dir/$file_name".bed"$suffix \
		> $output_dir/$file_name".bed"$suffix".name"

	awk '
		BEGIN{OFS="\t"}
		{split($4, namelowbar, "_")} 
		{$4=($4 !~ /CDS/) ? namelowbar[1] : $4}1
		' $output_dir/$file_name".bed"$suffix".name" \
		> $output_dir/$file_name".bed"

		rm $output_dir/$file_name".bed"$suffix
		rm $output_dir/$file_name".bed"$suffix".name"
fi

echo "$(date)"
echo "DONE getting bed file from GFF"
echo "File can be found at" $output_dir/$file_name".bed"
echo -e "\n"