#!/usr/bin/env python

# Standard library imports
import os
import sys
import logging

# Third party imports
import argparse
import datetime
import pandas as pd
import Bio
from Bio import Entrez
from Bio import SeqIO

logger = logging.getLogger()

"""
=============================================================
HEADER
=============================================================
FUNCTION: Download up to date plasmid database from https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/plasmids.txt.
            Remove those sequences with terms not related to complete plasmid such: gene, protein, partial, putative or hypothetical

INSTITUTION:CNM-ISCIII
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
d^v^b
VERSION=0.1
CREATED: 26 February 2020
REVISION: 

TODO:
    add user defined terms
    filter by record size (len(record))
================================================================
END_OF_HEADER
================================================================
"""


def check_create_dir(path):
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)


def main():

    def get_arguments():

        parser = argparse.ArgumentParser(
            prog='download_plasmid_database.py', description='Download up to date plasmid database from ncbi ftp')

        parser.add_argument('-o', '--output', type=str, required=True,
                            help='REQUIRED. Output directory to extract plasmid database')

        arguments = parser.parse_args()

        return arguments

    args = get_arguments()

    output_dir = os.path.abspath(args.output)

    check_create_dir(output_dir)

    # LOGGING
    # Create log file with date and time
    today = str(datetime.date.today())
    right_now_full = "".join(today.split("-"))

    log_filename = 'plasmidID_database' + "_" + right_now_full + ".log"
    log_full_path = os.path.join(output_dir, log_filename)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s:%(message)s')

    file_handler = logging.FileHandler(log_full_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    # stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    #####################START PIPELINE################

    logger.debug(args)

    plasmid_text_file = today + "_plasmids.txt"
    plasmid_text_path = os.path.join(output_dir, plasmid_text_file)

    plasmid_fasta_file = today + "_plasmids.fasta"
    plasmid_fasta_path = os.path.join(output_dir, plasmid_fasta_file)

    plasmid_failed_file = today + "failed_plasmids.txt"
    plasmid_failed_path = os.path.join(output_dir, plasmid_failed_file)

    try:
        df = pd.read_csv(
            'https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/plasmids.txt', sep='\t')
    except:
        logger.info('there was a problem accessing the ftp')
        sys.exit(1)

    df.to_csv(plasmid_text_path, sep='\t', index=False)

    plasmid_reference = df['RefSeq'][df.RefSeq !=
                                     "-"].tolist() + df['INSDC'][df.RefSeq == "-"].tolist()

    # remove duplicates
    plasmid_reference = set(plasmid_reference)
    # Set terms to exclude
    terms_to_exclude = ['gene ', 'protein',
                        'partial', 'putative', 'hypothetical']
    # Dictionary with erroneous accession numbers to determine the reason
    erroneous = {}

    Entrez.email = "A.N.Other@example.com"

    total_sequences = len(plasmid_reference)
    current_record = 1
    logger.info("")
    logger.info("Starting plasmid database download script: " +
                str(total_sequences) + " will be downloaded")
    logger.info("This will take a while.\nCheck progress in " + log_full_path)

    with open(plasmid_fasta_path, 'w+') as output_handle:
        for plasmid_accnumber in plasmid_reference:
            try:
                handle = Entrez.efetch(
                    db="nucleotide", id=plasmid_accnumber, rettype="fasta", retmode="text")
                record = SeqIO.read(handle, "fasta")
                terms_present = [
                    x in record.description for x in terms_to_exclude]
                handle.close()
                if sum(terms_present) > 0:
                    terms_true = [terms_to_exclude[i]
                                  for i, x in enumerate(terms_present) if x == True]
                    erroneous[record.id] = "Include terms: " + \
                        ', '.join(terms_true) + " => " + record.description
                    logger.debug(" %s/%s Invalid terms in record %s" %
                                 (current_record, total_sequences, record.id))
                else:
                    logger.debug(" %s/%s Downloading record %s" %
                                 (current_record, total_sequences, record.id))
                    SeqIO.write(record, output_handle, "fasta")
            except:
                logger.debug(" %s/%s Failed to download %s" %
                             (current_record, total_sequences, record.id))
                erroneous[record.id] = "failed to download"
            current_record = current_record + 1

    if len(erroneous) > 0:
        with open(plasmid_failed_path, 'w+') as ferror:
            for acc, reason in erroneous.items():
                ferror.write(acc + ": " + reason + "\n")

    logger.info("ALL DONE\nFASTA file is available in: " + plasmid_fasta_path)


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logger.exception(e)
        raise
