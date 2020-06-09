
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/fusion-report/README.html)
[![CircleCI Build Status](https://circleci.com/gh/circleci/circleci-docs.svg?style=shield)](https://circleci.com/gh/BU-ISCIII/plasmidID) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Scif](https://img.shields.io/badge/Filesystem-Scientific-brightgreen.svg)](https://sci-f.github.io)

# plasmidID <img align="left" src="https://github.com/BU-ISCIII/plasmidID/blob/develop/img/plasmidID_logo.png" alt="Logo" width="100">

<br>
<br>

* [Introduction](#introduction)
* [Requirements](#requirements)
    * [Software](#software)
    * [Plasmid database](#plasmid-database)
* [Installation](#installation)
    * [Install from source](#install-from-source)
    * [Install using conda](#install-using-conda)
* [Quick usage](#quick-usage)
* [Usage](#usage)
* [Output](#output)
* [Annotation file](#annotation-file)
* [Illustrated pipeline](#illustrated-pipeline)
* [Docker](#docker)

## Introduction

PlasmidID is a mapping-based, assembly-assisted plasmid identification tool that analyzes and gives graphic solution for plasmid identification.

PlasmidID is a **computational pipeline** implemented in **BASH** that maps Illumina reads over plasmid database sequences. The k-mer filtered, most covered sequences are clustered by identity to avoid redundancy and the longest are used as scaffold for plasmid reconstruction. Reads are assembled and annotated by automatic and specific annotation. All information generated from mapping, assembly, annotation and local alignment analyses is gathered and accurately represented in a **circular image** which allow user to determine plasmidic composition in any bacterial sample.

## Requirements

#### Software

* [Python >=3.6](https://www.python.org/)
* [Trimmomatic v0.33](http://www.usadellab.org/cms/?page=trimmomatic)(Optional)
* [Spades v3.8](http://cab.spbu.ru/software/spades/) (Optional)
* [Perl v5.26.0](https://www.perl.org/get.html)
* [NCBI_blast + v2.2.3](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [Bedtools v2.25](http://bedtools.readthedocs.io/en/latest/)
* [Bowtie 2 v2.2.4](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [SAMtools v1.2](http://samtools.sourceforge.net/)
* [prokka v1.12](http://www.vicbioinformatics.com/software.prokka.shtml)
* [cd-hit v4.6.6](http://weizhongli-lab.org/cd-hit/) (no longer needed since v1.6)
* [circos v0.69.3](http://circos.ca/software/download/circos/)
* [mash v2.2](https://github.com/marbl/Mash)

#### Plasmid database

Since version v1.5.1 plasmid database can be downloaded with the following command:

```Bash
 download_plasmid_database.py -o FOLDER
```

## Installation

#### Install from source

Install all dependencies and add them to $PATH

git clone https://github.com/BU-ISCIII/plasmidID.git

Add plasmidID and ./bin to $PATH

#### Install using conda

This option is recomended.

Install [Anaconda3](https://www.anaconda.com/distribution/)

```
conda install -c conda-forge -c bioconda plasmidid
```
Wait for the environment to solve

Ignore warnings/errors

## Quick usage

Illumina paired-end
```
plasmidID \
-1 SAMPLE_R1.fastq.gz  \
-2 SAMPLE_R2.fastq.gz \
-d YYYY-MM-DD_plasmids.fasta \
-c SAMPLE_assembled_contigs.fasta \
--no-trim \
-s SAMPLE
```

SMRT sequencing (only contigs)
```
plasmidID \
-d YYYY-MM-DD_plasmids.fasta \
-c SAMPLE_contigs.fasta \
-s SAMPLE
```

Annotate any fasta you want
```
plasmidID \
-d YYYY-MM-DD_plasmids.fasta \
-c SAMPLE_assembled_contigs.fasta \
-a annotation_file \
-s SAMPLE
```
More info about [annotation file](#annotation-file)

## Usage

```
usage : plasmidID <-1 R1> <-2 R2> <-d database(fasta)> <-s sample_name> [-g group_name] [options]

	Mandatory input data:
	-1 | --R1	<filename>	reads corresponding to paired-end R1 (mandatory)
	-2 | --R2	<filename>	reads corresponding to paired-end R2 (mandatory)
	-d | --database	<filename>	database to map and reconstruct (mandatory)
	-s | --sample	<string>	sample name (mandatory), less than 37 characters

	Optional input data:
	-g | --group	<string>	group name (optional). If unset, samples will be gathered in NO_GROUP group
	-c | --contigs	<filename>	file with contigs. If supplied, plasmidID will not assembly reads
	-a | --annotate <filename>	file with configuration file for specific annotation
	-o 		<output_dir>	output directory, by default is the current directory

	Pipeline options:
	--explore	Relaxes default parameters to find less reliable relationships within data supplied and database
	--only-reconstruct	Database supplied will not be filtered and all sequences will be used as scaffold
						This option does not require R1 and R2, instead a contig file can be supplied
	-w 			Undo winner takes it all algorithm when clustering by kmer - QUICKER MODE
	Trimming:
	--trimmomatic-directory Indicate directory holding trimmomatic .jar executable
	--no-trim	Reads supplied will not be quality trimmed

	Coverage and Clustering:
	-C | --coverage-cutoff	<int>	minimun coverage percentage to select a plasmid as scafold (0-100), default 80
	-S | --coverage-summary	<int>	minimun coverage percentage to include plasmids in summary image (0-100), default 90
	-f | --cluster	<int>	kmer identity to cluster plasmids into the same representative sequence (0 means identical) 		(0-1), default 0.5
	-k | --kmer	<int>	identity to filter plasmids from the database with kmer approach (0-1), default 0.95

	Contig local alignment
	-i | --alignment-identity <int>	minimun identity percentage aligned for a contig to annotate, default 90
	-l | --alignment-percentage <int>	minimun length percentage aligned for a contig to annotate, default 20
	-L | --length-total	<int>	minimun alignment length to filter blast analysis
	--extend-annotation <int>	look for annotation over regions with no homology found (base pairs), default 500bp

	Draw images:
	--config-directory <dir>	directory holding config files, default config_files/
	--config-file-individual <file-name> file name of the individual file used to reconstruct
	Additional options:

	-M | --memory	<int>	max memory allowed to use
	-T | --threads	<int>	number of threads
	-v | --version		version
	-h | --help		display usage message

example: ./plasmidID.sh -1 ecoli_R1.fastq.gz -2 ecoli_R2.fastq.gz -d database.fasta -s ECO_553 -G ENTERO
	./plasmidID.sh -1 ecoli_R1.fastq.gz -2 ecoli_R2.fastq.gz -d PacBio_sample.fasta -c scaffolds.fasta -C 60 -s ECO_60 -G ENTERO --no-trim
```

## Examples

Under construction

## Output

Since v1.6, the more relevant output is located in GROUP/SAMPLE folder:

- **SAMPLE_final_results.html(.tab)**
	- id: Name of the accession number of reference
	- length: length of the reference sequence
	- species: species of the reference sequence
	- description: rest of reference fasta header
	- contig_name: number of the contigs that align the minimun required for complete contig track
	- SAMPLE:
		- Image of the reconstructed plasmid (click to open in new tab)
		- MAPPING % (percentage): percentage of reference covered with reads
			- X for contig mode (gray colour)
			- Orientative colouring (the closer to 100% the better)
		- ALIGN FR (fraction_covered): total length of contigs aligned (complete) / reference sequence length
			- Orientative colouring (the closer to 1 the better)
			

## Annotation file

Under construction

## Illustrated pipeline

This image sumarizes PlasmidID pipeline, including the most important steps.
For furder details, including:
- [Results interpretation](https://github.com/BU-ISCIII/plasmidID/wiki/Understanding-the-image:-track-by-track)
- and more, please visit: [**PLASMIDID WIKI**](https://github.com/BU-ISCIII/plasmidID/wiki)

<p align="center"><img src="https://github.com/BU-ISCIII/plasmidID/blob/master/img/Short_pipeline.png" alt="workflow_small"  width="500">







## Docker

Example:
Clone the repo:
```Bash
git clone git@github.com:BU-ISCIII/plasmidID.git
cd plasmidID
git submodule init
git submodule update
```
Run it with the test data using docker or singularity:

**Notice that the input files MUST be in your present working directory or in any folder inside it. For example, if I execute this command in `/home/smonzon`, my folder with the files would be in `/home/smonzon/TEST_DATA`.**

```Bash
docker run -v $PWD:$PWD -w $PWD buisciii/plasmidid plasmidID.sh \
     -1 TEST_DATA/KPN_TEST_R1.fastq.gz  \
     -2 TEST_DATA/KPN_TEST_R2.fastq.gz \
     -d TEST_DATA/plasmids_TEST_database.fasta \
     -c TEST_DATA/contigs_KPN_TEST.fasta \
     --no-trim \
     -s KPN
```

Or you can use Singularity instead:
```Bash
singularity exec docker://buisciii/plasmidid plasmidID.sh \
     -1 TEST_DATA/KPN_TEST_R1.fastq.gz  \
     -2 TEST_DATA/KPN_TEST_R2.fastq.gz \
     -d TEST_DATA/plasmids_TEST_database.fasta \
     -c TEST_DATA/contigs_KPN_TEST.fasta \
     --no-trim \
     -s KPN
```

You can also build the image in your local machine cloning the repository and building it using the Singularity recipe provided:

```Bash
git clone git@github.com:BU-ISCIII/plasmidID.git
cd plasmidID
git submodule init
git submodule update
singularity build plasmidid.simg Singularity
singularity exec plasmidid.simg plasmidID.sh \
     -1 TEST_DATA/KPN_TEST_R1.fastq.gz  \
     -2 TEST_DATA/KPN_TEST_R2.fastq.gz \
     -d TEST_DATA/plasmids_TEST_database.fasta \
     -c TEST_DATA/contigs_KPN_TEST.fasta \
     --no-trim \
     -s KPN
```

If you want to use trimming functionality you have to supply trimmomatic-directory parameter:
```
docker run -v $PWD:$PWD -w $PWD buisciii/plasmidid plasmidID.sh \
     -1 file_R1.fastq.gz  \
     -2 file_R2.fastq.gz \
     -d database.fasta \
     -c contigs.fasta \
     --trimmomatic-directory /scif/apps/trimmomatic/Trimmomatic-0.38 \
     -s KPN
```
> **Note:** You MUST skip trimming filtering in order to use TEST_DATA as it is already filtered!
