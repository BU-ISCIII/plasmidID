# plasmidID <img align="left" src="https://github.com/BU-ISCIII/plasmidID/blob/master/img/plasmidID_logo.png" alt="Logo" width="100"> 

<br>
<br>

![#f03c15](https://placehold.it/15/f03c15/000000?text=+) **PlasmidID is already available to download** ![#f03c15](https://placehold.it/15/f03c15/000000?text=+)

PlasmidID is a mapping-based, assembly-assisted plasmid identification tool that analyzes and gives graphic solution for plasmid identification.

PlasmidID is a **computational pipeline** implemented in **BASH** that maps Illumina reads over plasmid database sequences. The most covered sequences are clustered by identity to avoid redundancy and the longest are used as scaffold for plasmid reconstruction. Reads are assembled and annotated by automatic and specific annotation. All information generated from mapping, assembly, annotation and local alignment analyses is gathered and accurately represented in a **circular image** which allow user to determine plasmidic composition in any bacterial sample.

This image sumarizes PlasmidID pipeline, including the most important steps.
For furder details, including database download, dependencies and execution, please visit: <br> [**PLASMIDID WIKI**](https://github.com/BU-ISCIII/plasmidID/wiki)

<img align="center" src="https://github.com/BU-ISCIII/plasmidID/blob/master/img/Short_pipeline.png"  width="300">


# Quick usage
Example:
```Bash
docker run buisciii/plasmidid plasmidID.sh \
     -1 TEST_DATA/KPN_TEST_R1.fastq.gz  \
     -2 TEST_DATA/KPN_TEST_R2.fastq.gz \
     -d TEST_DATA/plasmids_TEST_database.fasta \
     -c TEST_DATA/contigs_KPN_TEST.fasta \
     --no-trim \ 
     -s KPN
```
Or you can use Singularity instead:
```Bash
singularity exec buisciii/plasmidid plasmidID.sh \
     -1 TEST_DATA/KPN_TEST_R1.fastq.gz  \
     -2 TEST_DATA/KPN_TEST_R2.fastq.gz \
     -d TEST_DATA/plasmids_TEST_database.fasta \
     -c TEST_DATA/contigs_KPN_TEST.fasta \
     --no-trim \ 
     -s KPN
```
