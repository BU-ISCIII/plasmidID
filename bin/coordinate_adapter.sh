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
#CREATED: 17 May 2018
#REVISION:
#DESCRIPTION:coordinate_adapter script adapt coordinates obtained with a bed file to a reference sequences in a link file
#
#
#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

coordinate_adapter script adapt coordinates obtained with a bed file to a reference sequences in a link file

usage : $0 <-i inputfile(.bed)> <-l link_file> [-o <directory>] [-n <number>] [-f <file_name>] [-u] [-v] [-h]

	-i input file in bed format
	-l link file with coordinates relationship within bed file ddbb and link reference
	-o output directory (optional). By default the file is placed in the same location as input
	-n length to extend annotation, default 5000
	-f file name
	-u uniq mode. Remove duplicates
	-p prokka mode. Remove suffix of prokka 
	-v version
	-h display usage message

example: ./coordinate_adapter.sh -i genes.bed -l ecoli.links -n 10000

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
input_file="Bed_file"
link_file="Link_file"
number_extension=5000
unique=false
prokka_mode=false
suffix=""

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:l:n:f:puvh"
while getopts $options opt; do
	case $opt in
		i )
			input_file=$OPTARG
			;;
		l )
			link_file=$OPTARG
			;;
		o )
			output_dir=$OPTARG
			;;
		n )
			number_extension=$OPTARG
			;;
		f)
			file_name=$OPTARG
			;;
		u )
			unique=true
			suffix=".unique.tmp"
			;;
		p )
			prokka_mode=true
			suffix=".prokka.tmp"
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

check_mandatory_files.sh $input_file $link_file

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


echo "$(date)"
echo "adapting coordinates from" $input_file and $link_file
echo "file name is:" $file_name
#Create a dictionary file with all posibilities: Column 1 and 5 must have some common terms


awk 'NR==FNR{a[NR]=$1;b[NR]=$0;next}{for(i = 1; i <= NR; ++i){if (a[i] == $1) print b[i],"\t", $0}}' \
$input_file $link_file > $output_dir/$file_name".coordinates.tmp"


awk '(($2 >= $6 - '"${number_extension}"' && $2 <= $7) || ($3 >= $6 && $3 <= $7 + '"${number_extension}"')) {{isInverted=($10-$9); \
genelength=($3-$2)};{if (isInverted < 0) {coordChr1=(($7-$3)+$10);} else {coordChr1=(($2-$6)+$9)}}; \
coordChr2=(coordChr1+genelength); {print $8, coordChr1, coordChr2, $4}}' $output_dir/$file_name".coordinates.tmp" > $output_dir/$file_name".coordinates.negatives"

#resulting in a bed file with coordinated of plasmid bur refering to contig annotation:
#NZ_CP010574.1 34820 33528 arsB_1
#NZ_CP008930.1 90527 89235 arsB_1
#NZ_CP006927.1 44969 43677 arsB_1
#NZ_CP010574.1 81021 82508 ltrA_1
#NZ_CP008930.1 144220 145707 ltrA_1


#Remove duplicate of several matches 

awk '($2 > 0) && ($3 > 0)' $output_dir/$file_name".coordinates.negatives" \
> $output_dir/$file_name".coordinates"$suffix


if [ "$unique" == "true" ]; then
    awk '
        (!x[$1$4]++)
    	' $output_dir/$file_name".coordinates"$suffix \
		> $output_dir/$file_name".coordinates"

		rm $output_dir/$file_name".coordinates"$suffix
fi

if [ "$prokka_mode" == "true" ]; then

	awk '
		(!uniq[$1$4]++)
		' $output_dir/$file_name".coordinates"$suffix \
		> $output_dir/$file_name".coordinates.prokka.unique.tmp"

	awk '
		BEGIN{OFS="\t"}{split($4, namelowbar, "_")} {$4=($4 !~ /CDS/) ? namelowbar[1] : $4}1
		' $output_dir/$file_name".coordinates.prokka.unique.tmp" \
		> $output_dir/$file_name".coordinates"

	rm $output_dir/$file_name".coordinates.prokka.unique.tmp"
	rm $output_dir/$file_name".coordinates"$suffix
    
fi

rm $output_dir/$file_name".coordinates.tmp"
rm $output_dir/$file_name".coordinates.negatives"


echo "$(date)"
echo -e "Coordinates adapted to file" $output_dir/$file_name".coordinates" "\n"