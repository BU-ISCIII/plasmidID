annotation_config_file=$1

number_of_annotation=$(cat $annotation_config_file | awk '!/^#/ &&  NF > 0' | wc -l)
lines_to_annotate=$(cat $annotation_config_file | awk '!/^#/ &&  NF > 0 {print NR}')
echo "Number of databases to annotate:" $number_of_annotation
database_number=0

for i in $lines_to_annotate
do

	let database_number=database_number+1
	echo "##Annotating database" $database_number
	#check if config file exist with highlights and text

	#Declare variables
	ddbb_file=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $1}')
	ddbb_name=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $2}')
	percent_identity=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $3}')
	percent_aligment=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $4}')
	query_divisor=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $5}')
	query_side=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $6}')
	is_unique=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $7}')
	color_highlight=$(cat $annotation_config_file | awk 'BEGIN{FS=","} NR == "'$i'" {print $8}')


	#Check each variable and its position in order to rectify the correct one
	if [[ "$percent_identity" -lt 0  &&  "$percent_identity" -gt 100 ]]; then
		echo -e "Invalid value for P_IDENTITY. Please check $annotation_config_file\n"

		exit 1
	fi

	if [[ $percent_aligment -lt 0  &&  $percent_aligment -gt 100 ]]; then
		echo -e "Invalid value for P_ALIGNMENT. Please check $annotation_config_file\n"
		exit 1
	fi

	counter_variable=0
	for i in "$ddbb_file" "$ddbb_name" "$percent_identity" "$percent_aligment" "$query_divisor" "$query_side" "$is_unique" "$color_highlight"
	do
		let counter_variable=counter_variable+1
		if [[ "$i" != "$query_divisor" ]];then
			#echo $i " " $counter_variable
			if [ -z "$i" -a "$counter_variable" -eq 8 ];then
				echo "An input filed is missing in annotation" $database_number "please check annotation file with all fields separated by comma"
				exit 1
			elif [[ -z "$i" ]];then
				echo "EMPTY"
				echo "Input" $counter_variable "in annotation" $database_number "is missing, please check annotation file with all fields separated by comma"
				exit 1
			fi
		fi
	done

	#blast_align.sh -i $annotation_file -d $group/$sample/data/$sample".fna" -o $group/$sample/data -p annotation -f $sample
	#sample.annotation.blast

	#blast_to_bed.sh -i $group/$sample/data/$sample".annotation.blast" -b 95 -l 85 -d _ -D r -q _ -Q l
	#sample.annotation.bed

	#coordinate_adapter.sh -i $group/$sample/data/$sample".annotation.bed" -l $group/$sample/data/$sample".plasmids.blast.links"
	#sample.annotation.coordinates
done

#DDBBFILE	NANE	P_IDENTITY	P_ALIGNMENT	Q_DIVISOR	Q_SIDE_LR	IS_UNIQUE	COLOR