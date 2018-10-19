
[![CircleCI Build Status](https://circleci.com/gh/circleci/circleci-docs.svg?style=shield)](https://circleci.com/gh/BU-ISCIII/plasmidID) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![CircleCi Community](https://img.shields.io/badge/community-CircleCI%20Discuss-343434.svg)](https://discuss.circleci.com) [![Scif](https://img.shields.io/badge/Filesystem-Scientific-brightgreen.svg)](https://sci-f.github.io)

# plasmidID <img align="left" src="https://github.com/BU-ISCIII/plasmidID/blob/develop/img/plasmidID_logo.png" alt="Logo" width="100">

<br>
<br>

![#f03c15](https://placehold.it/15/f03c15/000000?text=+) **PlasmidID is already available to download** ![#f03c15](https://placehold.it/15/f03c15/000000?text=+)

PlasmidID is a mapping-based, assembly-assisted plasmid identification tool that analyzes and gives graphic solution for plasmid identification.

PlasmidID is a **computational pipeline** implemented in **BASH** that maps Illumina reads over plasmid database sequences. The most covered sequences are clustered by identity to avoid redundancy and the longest are used as scaffold for plasmid reconstruction. Reads are assembled and annotated by automatic and specific annotation. All information generated from mapping, assembly, annotation and local alignment analyses is gathered and accurately represented in a **circular image** which allow user to determine plasmidic composition in any bacterial sample.

This image sumarizes PlasmidID pipeline, including the most important steps.
For furder details, including:
- [Database download](https://github.com/BU-ISCIII/plasmidID/wiki/Plasmid-Database)
- [Dependencies](https://github.com/BU-ISCIII/plasmidID/wiki/Installation-and-Dependencies)
- [Execution](https://github.com/BU-ISCIII/plasmidID/wiki/Execution)
- [Results interpretation](https://github.com/BU-ISCIII/plasmidID/wiki/Understanding-the-image:-track-by-track)
- and more, please visit: [**PLASMIDID WIKI**](https://github.com/BU-ISCIII/plasmidID/wiki)

<p align="center"><img src="https://github.com/BU-ISCIII/plasmidID/blob/master/img/Short_pipeline.png" alt="workflow_small"  width="500">

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
```
