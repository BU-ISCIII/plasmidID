#!/usr/bin/env python

# Standard library imports
import os
import sys
import re
import logging
import subprocess

# Third party imports
import argparse
import datetime
import pandas as pd
import numpy as np 
from Bio import Entrez
from Bio import SeqIO

logger = logging.getLogger()

"""
=============================================================
HEADER
=============================================================
FUNCTION: Reduces redundancy in multifasta files using kmer mash distance,
takes the longest sequence per cluster as representative

INSTITUTION:CNM-ISCIII
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
d^v^b
VERSION=0.1
CREATED: 27 May 2020
REVISION: 

TODO:

================================================================
END_OF_HEADER
================================================================
"""

END_FORMATTING = '\033[0m'
WHITE_BG = '\033[0;30;47m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
BLUE =  '\033[34m'
CYAN = '\033[36m'
YELLOW = '\033[93m'
DIM = '\033[2m'

def check_create_dir(path):
    #exists = os.path.isfile(path)
    #exists = os.path.isdir(path)
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)

def mash_dist(input_file, output_dir, threads=10):
    # mash dist -i database.filtered_0.9_term.fasta database.filtered_0.9_term.fasta > database.filtered_0.9_term.mash.distances.tab
    #'reference_ID', 'query_ID', 'distance', 'p_value', 'shared_hashes'

    input_file = os.path.abspath(input_file)
    file_prefix = input_file.split('/')[-1]
    prefix = ('.').join(file_prefix.split('.')[0:-1])

    output_filename = prefix + ".mash.distances.tab"
    
    output_file = os.path.join(output_dir, output_filename)

    cmd = ["mash", "dist", "-i", "-p", str(threads), input_file, input_file]

    prog = cmd[0]
    param = cmd[1:]

    try:
        with open(output_file, "w+") as outfile:
            #calculate mash distance and save it in output file
            command = subprocess.run(cmd,
            stdout=outfile, stderr=subprocess.PIPE, universal_newlines=True)
        if command.returncode == 0:
            logger.debug(GREEN + "Program %s successfully executed" % prog + END_FORMATTING)
        else:
            logger.error(RED + BOLD + "Command %s FAILED\n" % prog + END_FORMATTING
                + BOLD + "WITH PARAMETERS: " + END_FORMATTING + " ".join(param) + "\n"
                + BOLD + "EXIT-CODE: %d\n" % command.returncode +
                "ERROR:\n" + END_FORMATTING + command.stderr)
    except OSError as e:
        sys.exit(RED + BOLD + "failed to execute program '%s': %s" % (prog, str(e)) + END_FORMATTING)

    return output_file

def mash_dist_to_pairwise(distance_file, distance_type='hash_distance'):
    df = pd.read_csv(distance_file, sep='\t', names=['reference_ID', 'query_ID', 'distance', 'p_value', 'shared_hashes'])
    df[['hash_1', 'hash_2']] = df['shared_hashes'].str.split('/', expand=True)
    df.hash_1 = df.hash_1.astype(float)
    df.hash_2 = df.hash_2.astype(float)
    df['hash_distance'] = 1 - (df.hash_1 / df.hash_2)
    dfpair = df[['reference_ID', 'query_ID', distance_type]]
    
    return dfpair

def pairwise_to_cluster(pw,threshold = 0.5):
    groups = {}
    columns = pw.columns.tolist()
    sorted_df = pw[(pw[columns[0]] != pw[columns[1]]) & (pw[columns[2]] <= threshold)].sort_values(by=[columns[2]])

    print(pw.head())
    print(sorted_df.shape)
    
    def rename_dict_clusters(cluster_dict):
        reordered_dict = {}
        for i, k in enumerate(list(cluster_dict)):
            reordered_dict[i] = cluster_dict[k]
        return reordered_dict
    
    def regroup_clusters(list_keys, groups_dict, both_samples_list):
        #sum previous clusters
        list_keys.sort()
        new_cluster = sum([groups_dict[key] for key in list_keys], [])
        #add new cluster
        cluster_asign = list(set(new_cluster + both_samples_list))
        #Remove duped cluster
        first_cluster = list_keys[0]
        groups_dict[first_cluster] = cluster_asign
        rest_cluster = list_keys[1:]
        for key in rest_cluster:
            del groups_dict[key]
        groups_dict = rename_dict_clusters(groups_dict)
        return groups_dict
        
    for _, row in sorted_df.iterrows():
        group_number = len(groups)
        sample_1 = str(row[0])
        sample_2 = str(row[1])
        both_samples_list = row[0:2].tolist()
                
        if group_number == 0:
            groups[group_number] = both_samples_list
        
        all_samples_dict = sum(groups.values(), [])
                
        if sample_1 in all_samples_dict or sample_2 in all_samples_dict:
            #extract cluster which have the new samples
            key_with_sample = {key for (key,value) in groups.items() if (sample_1 in value or sample_2 in value)}
            
            cluster_with_sample = list(key_with_sample)
            cluster_with_sample_name = cluster_with_sample[0]
            number_of_shared_clusters = len(key_with_sample)
            if number_of_shared_clusters > 1:
                groups = regroup_clusters(cluster_with_sample, groups, both_samples_list)
            else:
                groups[cluster_with_sample_name] = list(set(groups[cluster_with_sample_name] + both_samples_list))
        else:
            groups[group_number] = both_samples_list
            
    for _, row in pw[(pw[pw.columns[0]] != pw[pw.columns[1]]) & (pw[pw.columns[2]] > threshold)].iterrows():
        sample_1 = str(row[0])
        sample_2 = str(row[1])
        all_samples_dict = sum(groups.values(), [])
        if sample_1 not in all_samples_dict:
            group_number = len(groups)
            groups[group_number] = [sample_1]
        
        if sample_2 not in all_samples_dict:
            group_number = len(groups)
            groups[group_number] = [sample_2]
            
    cluster_df = pd.DataFrame(groups.values(),index=list(groups))
    
    cluster_df_return = cluster_df.stack().droplevel(1).reset_index().rename(columns={'index': 'group', 0: 'id'})

    cluster_df_return.to_csv('/processing_Data/antibioticos/mperezv/ANALYSIS/Polimixinas_OTA/plasmidID/NO_GROUP/8c/mapping/test_cluster.tab', sep='\t')

    return cluster_df_return

def big_pairwise_to_cluster(pw,threshold = 0.5):
    
    def rename_dict_clusters(cluster_dict):
        reordered_dict = {}
        for i, k in enumerate(list(cluster_dict)):
            reordered_dict[i] = cluster_dict[k]
        return reordered_dict
    
    def regroup_clusters(list_keys, groups_dict, both_samples_list):
        #sum previous clusters
        list_keys.sort()
        new_cluster = sum([groups_dict[key] for key in list_keys], [])
        #add new cluster
        cluster_asign = list(set(new_cluster + both_samples_list))
        #Remove duped cluster
        first_cluster = list_keys[0]
        groups_dict[first_cluster] = cluster_asign
        rest_cluster = list_keys[1:]
        for key in rest_cluster:
            del groups_dict[key]
        groups_dict = rename_dict_clusters(groups_dict)
        return groups_dict
    
    groups = {}
    
    with open(pw, "r") as f:
        for line in f:
            line_split = line.split('\t')
            sample_1 = line_split[0]
            sample_2 = line_split[1]
            hash1, hash2 = line_split[4].split('/')
            hash_distance = 1 - (int(hash1) / int(hash2))

            if hash_distance <= threshold:
                group_number = len(groups)

                both_samples_list = [sample_1,sample_2]

                if group_number == 0:
                    groups[group_number] = both_samples_list

                all_samples_dict = sum(groups.values(), [])

                if sample_1 in all_samples_dict or sample_2 in all_samples_dict:
                    #extract cluster which have the new samples
                    key_with_sample = {key for (key,value) in groups.items() if (sample_1 in value or sample_2 in value)}

                    cluster_with_sample = list(key_with_sample)
                    cluster_with_sample_name = cluster_with_sample[0]
                    number_of_shared_clusters = len(key_with_sample)
                    if number_of_shared_clusters > 1:
                        groups = regroup_clusters(cluster_with_sample, groups, both_samples_list)
                    else:
                        groups[cluster_with_sample_name] = list(set(groups[cluster_with_sample_name] + both_samples_list))
                else:
                    groups[group_number] = both_samples_list
            else:
                if sample_1 not in all_samples_dict:
                    group_number = len(groups)
                    groups[group_number] = [sample_1]

                if sample_2 not in all_samples_dict:
                    group_number = len(groups)
                    groups[group_number] = [sample_2]
            
    cluster_df = pd.DataFrame(groups.values(),index=list(groups))
    
    cluster_df_return = cluster_df.stack().droplevel(1).reset_index().rename(columns={'index': 'group', 0: 'id'})

    cluster_df_return.to_csv('/processing_Data/bioinformatics/references/plasmidID/plasmid_ddbb/20200203/test_ddbb.csv')
            
    return cluster_df_return
            
def calculate_seqlen(fasta_file):
    len_df = pd.DataFrame(columns=['id','length'])
    index = 0
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        len_df.loc[index, 'id'] = seq_record.id
        len_df.loc[index, 'length'] = len(seq_record)
        index = index + 1
    len_df['length'] = len_df.length.astype('int')
    return len_df
    
def extract_representative(row):
    row.clustered.remove(row.id)
    return

def extract_length(row, final_cluster):
    lengths = [final_cluster['length'][final_cluster.id == idclust].tolist()[0] for idclust in row.clustered]
    return lengths

def extract_distance_legacy(reprensetative, list_clustered, pwdist):
    distances = [round(pwdist[pwdist.columns[2]][(pwdist[pwdist.columns[0]] == reprensetative) & (pwdist[pwdist.columns[1]] == idclust)].tolist()[0],2) for idclust in list_clustered]
    return distances

def extract_distance(reprensetative, list_clustered, mash_file):
    distances = []
    for idclust in list_clustered:
        with open(mash_file, "r") as f:
            for line in f:
                line_split = line.split('\t')
                sample_1 = line_split[0]
                sample_2 = line_split[1]
                if sample_1 == reprensetative and sample_2 == idclust:
                    hash1, hash2 = line_split[4].split('/')
                    hash_distance = 1 - (int(hash1) / int(hash2)) 
                    distances.append(hash_distance)
    return distances

def retrieve_fasta_cluster(fasta_file, final_cluster, output_dir, mash_file, kmerdist, quantity_id=1, save_clustered=False):
    input_file = os.path.abspath(fasta_file)
    file_prefix = input_file.split('/')[-1]
    prefix = ('.').join(file_prefix.split('.')[0:-1])
    
    output_representative = os.path.join(output_dir, prefix + '.' + str(kmerdist) + '.representative.fasta')
    output_clustered = os.path.join(output_dir, prefix +  '.' + str(kmerdist) + '.clustered.fasta')
    output_summary = os.path.join(output_dir, prefix + '.' + str(kmerdist) + '.clusters.tab')
    
    representative_id = final_cluster.sort_values(by=['group', 'length'], ascending=[True, False]).groupby('group').head(quantity_id)
    summary_id_grouped = final_cluster.groupby('group')['id'].apply(list).reset_index(name='clustered')
    representative_list = representative_id.id.tolist()
    representative_and_sumary = representative_id.merge(summary_id_grouped, on='group', how='left')
    #Use function extract_representative to remove the repr. from column
    representative_and_sumary.apply(extract_representative, axis=1)
    representative_and_sumary['lengths_clustered'] = representative_and_sumary.apply(lambda x: extract_length(x, final_cluster), axis=1)
    representative_and_sumary['distance_clustered'] = representative_and_sumary.apply(lambda x: extract_distance(x.id, x.clustered, mash_file), axis=1)
    #read the fasta and retrieve representative sequences
    representative_records = []
    clustered_records = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        if seq_record.id in representative_list:
            representative_records.append(seq_record)
        else:
            clustered_records.append(seq_record)
        
    SeqIO.write(representative_records, output_representative, "fasta")
    
    if not save_clustered == False:
        SeqIO.write(clustered_records, output_clustered, "fasta")
        
    representative_and_sumary.to_csv(output_summary, sep='\t', index=False)

    previous_seq = final_cluster.shape[0]
    post_seq = representative_and_sumary.shape[0]

    logger.info(MAGENTA + "%s sequences clustered into %s" % (previous_seq, post_seq) + END_FORMATTING)



    

def main():

    def get_arguments():

        parser = argparse.ArgumentParser(prog = 'common_mash_reference.py', description= 'Search for all mash files and find the representative reference')
        
        parser.add_argument('-i', '--input', dest="input_file", metavar="input_directory", type=str, required=True, help='REQUIRED.Input FASTA file')
        parser.add_argument('-o', '--output', type=str, required=False, default=False, help='Output directory to extract clusteres FASTA')
        parser.add_argument('-d', '--distance', type=float, required=False, default=0.5, help='Threshold distance to cluster sequences[0-1] 0(identical) 1(unrelated) (default 0.5)')
        parser.add_argument('-g', '--output_grouped', required=False,  action='store_true', help='Output clustered (non representative sequences')

        arguments = parser.parse_args()

        return arguments

    args = get_arguments()

    input_file = os.path.abspath(args.input_file)

    if args.output == False:
        output_dir = ('/').join(input_file.split('/')[0:-1])
    else:
        output_dir = os.path.abspath(args.output)
    
    check_create_dir(output_dir)


    #LOGGING
    #Create log file with date and time
    right_now = str(datetime.date.today())
    right_now_full = "_".join(right_now.split(" "))

    log_filename = 'mashclust' + "_" + right_now_full + ".log"
    log_full_path = os.path.join(output_dir, log_filename)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s:%(message)s')

    file_handler = logging.FileHandler(log_full_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    #stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    #####################START PIPELINE################

    logger.info(args)
    #CALCULATE MASH DISTANCE
    logger.info('Obtaining mash distance')
    mash_file = mash_dist(input_file, output_dir, threads=10)
    logger.info('Obtaining cluster from distance')
    #pairwise_distance = mash_dist_to_pairwise(mash_file)
    #cluster_df = pairwise_to_cluster(pairwise_distance, threshold=args.distance)
    cluster_df = big_pairwise_to_cluster(mash_file, threshold=args.distance)
    logger.info('Calculating length')
    len_df = calculate_seqlen(input_file)
    final_cluster = cluster_df.merge(len_df, on='id', how='left')
    logger.info('Filtering representative fasta')
    retrieve_fasta_cluster(input_file, final_cluster, output_dir, mash_file, args.distance, quantity_id=1, save_clustered=args.output_grouped)
    logger.info('DONE')

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logger.exception(e)
        raise