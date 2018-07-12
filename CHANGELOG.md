# Changelog
All notable changes to this project will be documented in this file.

## 1.3.0 - 2018-07-11
### New
- Summary table can be generated with new utility
- Several databases can be now annotated filling annotation_config_file.txt
### Fixed
- circos dependency is now checked

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