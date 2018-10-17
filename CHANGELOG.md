# Changelog
All notable changes to this project will be documented in this file.


## 1.4.2 - 2018-09-29
### Added
- Specific config file for only reconstruct parameter

###Fixed
- Protein databases can be properly used

## 1.4 - 2018-09-20
### Added
- Automatically annotated genes/cds are displayed differently depending on whether they are located in forward or reverse 
- Psi-cd-hit and blast now handle threads
- Improved error handling
- Doocker/Singularity compatibility
- One multifasta file per reference plasmid is generated with all the similar contigs from the sample 
- Quick staus of values applied to plasmid reconstruction

###Fixed
- Some plasmids from the database were not annotated
- Limit sample name to 37 characters, capped by prokka
- Bug in complete contig track generator that took the wrong value and couldn't draw sequences that matched the position 0 of plasmid



## 1.3.0 - 2018-07-11
### New
- Summary table can be generated with new utility
- Several databases can be now annotated filling annotation_config_file.txt
- --only-reconstruct is now implemented if user only needs to reconstruct and annotate contigs with small known databases
### Fixed
- circos dependency is now checked
- Output is now correctly redirected with -o
### Added
- trimmomatic directory containing .jar can no be especified with --trimmomatic-directory
- Vervose mode included. By default a log file will be created
- Friendly terminal output

## 1.2.2 - 2018-06-22
### Fixed
- ***IMPORTANT***: PlasmidID maps with -a mode NOW, as it should have allways been. A bug on mapping script is now solved
- Number of threads are now implemented on mapping
- Some cumulative clustering temporary files are now removed

## 1.2.1 - 2018-06-14
### Fixed
- All dependencies are now checked at the beggining
- Path to scripts are no longer hard coded paths
- Links should be now displayed on summary image

### Added
- Added first utility ***ncbi_database_fetcher.sh***, a script to download FASTA databases from terms
- Short scripts now moved to /bin has to be added to PATH


## 1.1.1 - 2018-06-11
### Fixed
- Additional database will not be required for circos executios, even though the file will be created
- Fixed an issue when no plasmid matches mapping requeriments
- Fixed an issue when circos will trow an error message when no plasmids met mapping requeriments


## 1.1.0 - 2018-06-06
### Added
- Database plasmids used as scaffold are annotated after filtering. User doesn't need to annotate the initial huge plasmid database.
- User can add ONE nucleotide FASTA file wi that will be specifically annotated on final plasmids with a light blue color

## Unreleased

- Create config files as required by user and include visual parameters
- Test and adapt the --only-reconstruct option