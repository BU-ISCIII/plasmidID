#!/bin/bash


echo -e "\n#Executing" $0 "\n"

group=$1
sample=$2
cdsDdbb_file=$3



mappedDir=$group/$sample/mapping
imageDir=$group/$sample/data

circos_conf_summary="config_files/circos_summary.conf"
circos_conf_individual="config_files/circos_individual.conf"
circosDir=$group/$sample/images

plasmidMapped=$mappedDir/$sample".coverage_adapted_clustered_ac"

karyotype_file_individual=$imageDir/$sample".karyotype_individual.txt"
karyotype_file_summary=$imageDir/$sample".karyotype_summary.txt"
abr_file=$imageDir/$sample".abr.coordinates"
replisome_file=$imageDir/$sample".inc.coordinates"
additional_file=$imageDir/$sample".annotation.coordinates"
coverage_file=$imageDir/$sample".bedgraph_term"
cds_contig_file=$imageDir/$sample".gff.coordinates"

contig_file=$imageDir/$sample".plasmids.bed"
contig_file_complete=$imageDir/$sample".plasmids.complete"
links_file=$imageDir/$sample".plasmids.links"

imageName=$sample"_summary.png"

mkdir -p $circosDir


echo "Creating config file for circos in SAMPLE $sample FILE $circosDir/$sample.circos.conf"

awk '{gsub("PLASMID_KARYOTYPE","'$karyotype_file_individual'"); \
gsub("PLASMID_ANTIBIOTIC_RESISTANCE","'$abr_file'"); \
gsub("PLASMID_REPLISOME_PLASMIDFINDER","'$replisome_file'"); \
gsub("PLASMID_ANNOTATION_USER","'$additional_file'"); \
gsub("PLASMID_COVERAGE_GRAPH","'$coverage_file'"); \
gsub("PLASMID_CDS_CONTIG","'$cds_contig_file'"); \
gsub("PLASMID_CDS_DDBB","'$cdsDdbb_file'"); \
gsub("PLASMID_CONTIGS_COMPLETE","'$contig_file_complete'"); \
gsub("PLASMID_CONTIGS","'$contig_file'"); \
gsub("PLASMID_LINKS","'$links_file'"); \
gsub("OUTPUTDIR","'$circosDir'"); \
print $0}' $circos_conf_individual > $circosDir/$sample"_individual.circos.conf"

echo "DONE Creating config file for circos in SAMPLE $sample"

echo "Executing circos in SAMPLE $sample FILE $circosDir/$sample.circos.conf"




for i in $(cat $plasmidMapped)
	do
		awk '{gsub("SAMPLE_SHOWN","'$i'"); \
		gsub("IMAGENAME_SAMPLE_PLASMID","'$sample'_'$i'.png"); \
		print $0}' $circosDir/$sample"_individual.circos.conf" > $circosDir/$sample"_"$i"_individual.circos.conf"
		circos -conf $circosDir/$sample"_"$i"_individual.circos.conf"

	done 


if [ -s $karyotype_file_summary ]; then

	echo "Creating config file for circos in SAMPLE $sample FILE $circosDir/$sample.circos.conf"

	awk '{gsub("PLASMID_KARYOTYPE","'$karyotype_file_summary'"); \
	gsub("PLASMID_ANTIBIOTIC_RESISTANCE","'$abr_file'"); \
	gsub("PLASMID_REPLISOME_PLASMIDFINDER","'$replisome_file'"); \
	gsub("PLASMID_ANNOTATION_USER","'$additional_file'"); \
	gsub("PLASMID_COVERAGE_GRAPH","'$coverage_file'"); \
	gsub("PLASMID_CDS_CONTIG","'$cds_contig_file'"); \
	gsub("PLASMID_CDS_DDBB","'$cdsDdbb_file'"); \
	gsub("PLASMID_CONTIGS","'$contig_file'"); \
	gsub("PLASMID_LINKS","'$links_file'"); \
	gsub("OUTPUTDIR","'$circosDir'"); \
	gsub("IMAGENAME","'$imageName'"); \
	print $0}' $circos_conf_summary > $circosDir/$sample"_summary.circos.conf"

	echo "DONE Creating config file for circos in SAMPLE $sample"


	circos -conf $circosDir/$sample"_summary.circos.conf"

else

	echo "No plasmid mathed requirements to draw the summary image"

fi


#Remove config files
<<R 
for i in $(cat $plasmidMapped)
do
	if [ -f $circosDir/$sample"_"$i"_individual.circos.conf" ]; then
		rm $circosDir/$sample"_"$i"_individual.circos.conf"
	fi
	
done

	rm $circosDir/$sample"_summary.circos.conf"
	rm $circosDir/$sample"_individual.circos.conf"
R
echo "ALL DONE"
