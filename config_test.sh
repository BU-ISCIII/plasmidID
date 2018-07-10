annotation_config_file=$1

number_of_annotation=$(cat $annotation_config_file | awk '!/^#/ &&  NF > 0' | wc -l)
lines_to_annotate=$(cat $annotation_config_file | awk '!/^#/ &&  NF > 0 {print NR}')
echo "Number of databases to annotate:" $number_of_annotation


for i in $lines_to_annotate
do
	#check if config file exist with highlights and text

	#Declare variables

	ddbb_file=$(cat $annotation_config_file | awk 'BEGIN{FS="\t"} NR == "'$i'" {print $1}')
	ddbb_name=$(cat $annotation_config_file | awk 'BEGIN{FS="\t"} NR == "'$i'" {print $2}')
	percent_identity=$(cat $annotation_config_file | awk 'BEGIN{FS="\t"} NR == "'$i'" {print $3}')
	percent_aligment=$(cat $annotation_config_file | awk 'BEGIN{FS="\t"} NR == "'$i'" {print $4}')
	


	ddbb_name=$(cat $annotation_config_file | awk 'BEGIN{FS="\t"} NR == "'$i'" {print $2}')
	ddbb_name=$(cat $annotation_config_file | awk 'BEGIN{FS="\t"} NR == "'$i'" {print $2}')
	ddbb_name=$(cat $annotation_config_file | awk 'BEGIN{FS="\t"} NR == "'$i'" {print $2}')
	ddbb_name=$(cat $annotation_config_file | awk 'BEGIN{FS="\t"} NR == "'$i'" {print $2}')




	echo $ddbb_file
	echo $ddbb_name
	echo $percent_identity
	echo $percent_aligment

	#blast_align.sh -i $annotation_file -d $group/$sample/data/$sample".fna" -o $group/$sample/data -p annotation -f $sample
	#sample.annotation.blast

	#blast_to_bed.sh -i $group/$sample/data/$sample".annotation.blast" -b 95 -l 85 -d _ -D r -q _ -Q l
	#sample.annotation.bed

	#coordinate_adapter.sh -i $group/$sample/data/$sample".annotation.bed" -l $group/$sample/data/$sample".plasmids.blast.links"
	#sample.annotation.coordinates
done

#DDBBFILE	NANE	P_IDENTITY	P_ALIGNMENT	Q_DIVISOR	Q_SIDE_LR	IS_UNIQUE	COLOR