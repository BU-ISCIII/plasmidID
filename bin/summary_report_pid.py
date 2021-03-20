#!/usr/bin/env python

# Standard library imports
import os
import sys
import re
import logging
import subprocess
import html


# Third party imports
import argparse
import datetime
import pandas as pd
import numpy as np
from Bio import Entrez
from Bio import SeqIO
from tabulate import tabulate

logger = logging.getLogger()

"""
=============================================================
HEADER
=============================================================
FUNCTION: Creates a summary report in tsv and hml from plasmidID execution

INSTITUTION:CNM-ISCIII
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
d^v^b
VERSION=0.1
CREATED: 04 June 2020
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

def check_file_exists(file_name):
    """
        Check file exist and is not 0 Kb, if not program exit.
    """
    #file_info = os.stat(file_name) #Retrieve the file info to check if has size > 0
    #or file_info.st_size == 0:
    if not os.path.isfile(file_name):
        logger.info(RED + BOLD + "File: %s not found or empty\n" % file_name + END_FORMATTING)
        sys.exit(1)
    return os.path.isfile(file_name)

def extract_files(folder):
    percentage_file = ""
    complete_file = ""
    representative_file = ""
    for root, _, files in os.walk(folder):
        for name in files:
            if name.endswith(".coverage_adapted_clustered_percentage"):
                percentage_file = os.path.join(root, name)
            elif name.endswith(".plasmids.complete"):
                complete_file = os.path.join(root, name)
            elif "representative" in name and name.endswith(".fasta"):
                representative_file = os.path.join(root, name)

    return percentage_file, complete_file, representative_file

def percentage_to_df(percentage_file):
    if not percentage_file == "":
        df = pd.read_csv(percentage_file, sep=" ", names=['id', 'percentage'])
        df['percentage'] = df['percentage'].round(2)
        return df
    else:
        return pd.DataFrame(columns=['id','percentage'])

def len_description_to_df(representative_file):
    df = pd.DataFrame(columns=['id','length','species', 'description'])
    index = 0
    for seq_record in SeqIO.parse(representative_file, "fasta"):
        df.loc[index, 'id'] = seq_record.id
        df.loc[index, 'length'] = len(seq_record)
        description_split = seq_record.description.split(' ')
        df.loc[index, 'species'] = (' ').join(description_split[1:3])
        df.loc[index, 'description'] = (' ').join(description_split[2:])
        index = index + 1
    df['length'] = df['length'].astype(int)

    return df

def complete_report_df(complete_file, len_description_df, percentage_df):

    def set_to_list(row):
        listed_set = list(row.contig_name)
        listed_set.sort()
        return listed_set

    #CP029217.1	176762	288994	9	id=170244
    dfc = pd.read_csv(complete_file, sep="\t", names=['id', 'start', 'end', 'contig_name', 'contig_id'])
    dfc['len_covered'] = dfc.end - dfc.start
    covered_df = dfc.groupby('id')['len_covered'].sum().reset_index()
    contigs_df = dfc.groupby('id')['contig_name'].apply(set).reset_index()#Merge all dataframes
    #Merge all dataframes
    df = len_description_df.merge(covered_df, on='id', how='left')
    df['fraction_covered'] = round(df.len_covered / df.length, 2)
    del df['len_covered']
    df = df.merge(contigs_df, on='id', how='left')
    df = df.dropna()
    df['contig_name'] = df.apply(lambda x: set_to_list(x), axis=1)
    df = df.merge(percentage_df, on='id', how='left')
    df = df.sort_values(by=['length'], ascending=False).reset_index(drop=True)
    df = df.fillna('X')
    return df

def include_images(sample_folder, summary_df):
    sample = sample_folder.split("/")[-1]
    def image_finder(row, sample_folder):
        for root, _, files in os.walk(sample_folder):
            for name in files:
                if 'images' in root and row.id in name and name.endswith('.png'):
                    return os.path.relpath(os.path.join(root, name), sample_folder)

    summary_df['images'] = summary_df.apply(lambda x: image_finder(x, sample_folder), axis=1)
    summary_df.to_csv(sample_folder + '/' + sample +  '_final_results.tab', sep='\t', index=False)
    return summary_df

html_template = """
    <!doctype html>

    <html lang="en">

    <style type="text/css">

    body {
        font: normal 20px "Courier New", Arial, sans-serif;
        padding: auto;
        margin: auto;
    }
    img {
        display: block;
        max-width: 350px;
        height: auto;
    }

    .summary{
        display: flex;
        flex-direction: column-reverse;
    }

    .numeric-values {
        display: flex;
        flex-direction: row;
        justify-content: space-around;
    }

    .numeric-values div {
        flex-grow: 1;
    }

    .neutral {
        background-color: lightgray;
    }

    .likely {
        background-color: limegreen;
    }

    .unlikely {
        background-color: sandybrown;
    }

    .unprobable {
        background-color: brown;
    }

    header {
        display: flex;
        flex-direction: row;
        justify-content: flex-start;
        align-items: flex-end;
        font-size: 2.7em;
    }

    table {
        position: relative;
        border-collapse: collapse;
    }

    header img {
        width: 70px;
    }

    th {
        font-size: 1.3em;
        background-color: skyblue;
        position: sticky;
        top: 0;
    }

    tr td {
        font-size: 1.1em;
        text-align: center;
    }

    footer {
        height: 3em;
        padding: 1;
    }

    tr:nth-child(even) {background-color: snow;}
    tr:hover {background-color:azure;}

    </style>

    <head>
      <meta charset="utf-8">

      <title>PlasmidID Report</title>
      <meta name="description" content="https://github.com/BU-ISCIII/plasmidID">
      <meta name="author" content="pedroscampoy@gmail.com">

      <link rel="stylesheet" href="css/styles.css?v=1.0">
      <link rel="shortcut icon" type="image/png" href="https://raw.github.com/BU-ISCIII/plasmidID/master/img/plasmidID_logo.png">

    </head>

    <body>
      <header>
         <img src="https://raw.github.com/BU-ISCIII/plasmidID/master/img/plasmidID_logo.png"' alt="header-image">
         PlasmidID REPORT
      </header>
      <div>
      TABLESUMMARY
      </div>
      <footer>
          More information at <a href="https://github.com/BU-ISCIII/plasmidID" target="_blank"> https://github.com/BU-ISCIII/plasmidID </a>
          <br>
          Author: pedroscampoy@gmail.com (BU-ISCIII)
      </footer>
    </body>
    </html>


    \n"""

def summary_to_html(sample_folder, final_individual_dataframe, html_template):
    df = final_individual_dataframe.copy()
    sample = sample_folder.split("/")[-1]
    html_filename = os.path.join(sample_folder, sample + '_final_results.html')
    hidden_filename = os.path.join(sample_folder, '.' + sample + '_final_individual_results.tab')

    def complete_to_rating(row):
        if row.fraction_covered >= 0.8 and row.fraction_covered <= 1.2:
            return 'likely'
        elif row.fraction_covered > 1.2 or (row.fraction_covered < 0.8 and row.fraction_covered > 0.5):
            return 'unlikely'
        else:
            return 'unprobable'

    def mapping_to_rating(row):
        if row.percentage == 'X':
            return 'neutral'
        elif row.percentage >= 80:
            return 'likely'
        elif row.percentage < 80 and row.percentage > 60:
            return 'unlikely'
        else:
            return 'unprobable'


    def apply_img_tag(row):
        return '<div class=summary>' + '\n' + \
    '<div class=numeric-values>' + '\n' + \
    '<div class=\"percentage ' + row.perc_rating + '\">' + 'MAPPING %<br>' + str(row.percentage) + '</div>' + '\n' + \
    '<div class=\"complete ' + row.complete_rating + '\">' + 'ALIGN FR<br>' + str(row.fraction_covered) + '</div>' + '\n' + \
    '</div>' + '\n' + \
    '<a href=' + row.images + ' target=\"_blank\">' + '\n' + \
    '<img src=' + row.images + ' alt=' + "\"" + row.id + "\"" + '>' + '\n' + \
    '</a>' + '\n' + \
    '</div>'

    def italic_species(row):
        return '<i>' + row.species + '</i>'

    df['perc_rating'] = df.apply(lambda x: mapping_to_rating(x), axis=1)

    df['complete_rating'] = df.apply(lambda x: complete_to_rating(x), axis=1)

    df['images'] = df.apply(lambda x: apply_img_tag(x), axis=1)

    df['species'] = df.apply(lambda x: italic_species(x), axis=1)

    df.drop(['percentage', 'fraction_covered', 'perc_rating', 'complete_rating'], axis = 1, inplace = True)

    df.rename(columns={'images':sample}, inplace=True)

    df.to_csv(hidden_filename, sep='\t', index=False)

    table = tabulate(df, headers='keys', tablefmt='html', showindex=False)
    table = html.unescape(table)
    table = table.replace("style=\"text-align: right;\"", "")

    final_html = html_template.replace('TABLESUMMARY', table)
    with open(html_filename, 'w+') as f:
        f.write(final_html)

def summary_to_html_group(group_folder, html_template):
    group = group_folder.split("/")[-1]
    html_filename = os.path.join(group_folder, group + '_final_results.html')
    individual_files = []
    for root, _, files in os.walk(group_folder):
        for name in files:
            if name.endswith("final_individual_results.tab") and name.startswith("."):
                individual_files.append(os.path.join(root, name))
    individual_dfs = []
    sample_list_column = []
    for file in individual_files:
        df = pd.read_csv(file, sep='\t')
        del df['contig_name']
        individual_dfs.append(df)
        sample_list_column.append(df.columns.tolist()[-1])

    dfm = individual_dfs[0]
    for df_ in individual_dfs[1:]:
        dfm = dfm.merge(df_, on=['id','length', 'species', 'description'], how='outer')

    coun_df = dfm.drop(['length', 'species', 'description'], axis = 1).groupby('id').count().sum(axis=1).reset_index(name='N')

    dfm = dfm.merge(coun_df, on='id', how='outer')

    columns_reorder = ['id','length', 'species', 'description', 'N'] + sample_list_column

    dfm = dfm[columns_reorder]

    dfm.fillna('-', inplace=True)

    dfm = dfm.sort_values(by=['N','length'], ascending=[False,False]).reset_index(drop=True)

    table = tabulate(dfm, headers='keys', tablefmt='html', showindex=False)
    table = html.unescape(table)

    final_html = html_template.replace('TABLESUMMARY', table)

    with open(html_filename, 'w+') as f:
        f.write(final_html)

    return dfm

def summary_to_tab_group(group_folder):
    group = group_folder.split("/")[-1]
    tab_filename = os.path.join(group_folder, group + '_final_results.tab')
    individual_files = []
    for root, _, files in os.walk(group_folder):
        for name in files:
            if name.endswith("final_results.tab"):
                individual_files.append(os.path.join(root, name))
    individual_dfs = []

    for file in individual_files:
        sample = file.split('/')[-1].replace('_final_results.tab', '')
        df = pd.read_csv(file, sep='\t')
        df.drop(['contig_name', 'images'], axis = 1, inplace=True)
        df.rename(columns={'fraction_covered':'Fr_cov_' + sample, 'percentage':'Map%_' + sample}, inplace=True)
        individual_dfs.append(df)

    dfm = individual_dfs[0]
    for df_ in individual_dfs[1:]:
        dfm = dfm.merge(df_, on=['id','length', 'species', 'description'], how='outer')

    count_df = dfm.filter(regex='Fr_cov.*|^id$', axis=1).groupby('id').count().sum(axis=1).reset_index(name='N')

    dfm = dfm.merge(count_df, on='id', how='outer')

    columns_reorder = dfm.columns.tolist()[0:4]
    columns_reorder.append(dfm.columns.tolist()[-1])
    columns_reorder = columns_reorder + dfm.columns.tolist()[4:-1]

    dfm = dfm[columns_reorder]

    dfm = dfm.sort_values(by=['N','length'], ascending=[False,False]).reset_index(drop=True)

    dfm.to_csv(tab_filename, sep='\t', index=False)

    return


def main():

    def get_arguments():

        parser = argparse.ArgumentParser(prog = 'summary_report_pid.py', description= 'Creates a summary report in tsv and hml from plasmidID execution')

        parser.add_argument('-i', '--input', dest="input_folder", metavar="input_folder", type=str, required=True, help='REQUIRED.Input pID folder')
        parser.add_argument('-g', '--group', required=False,  action='store_false', help='Creates a group report instead of individual (Default True)')


        arguments = parser.parse_args()

        return arguments

    args = get_arguments()

    input_folder = os.path.abspath(args.input_folder)
    #output_dir = input_folder

    #LOGGING
    #Create log file with date and time
    #right_now = str(datetime.date.today())
    #right_now_full = "_".join(right_now.split(" "))

    #log_filename = 'logs/summary_pid' + "_" + right_now_full + ".log"
    #log_full_path = os.path.join(output_dir, log_filename)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    #formatter = logging.Formatter('%(asctime)s:%(message)s')

    #file_handler = logging.FileHandler(log_full_path)
    #file_handler.setLevel(logging.DEBUG)
    #file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    #stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    #logger.addHandler(file_handler)

    #####################START PIPELINE################

    logger.info(args)
    #CALCULATE MASH DISTANCE
    logger.info('Creating summary')

    if args.group == True:
        summary_to_html_group(input_folder, html_template)
        summary_to_tab_group(input_folder)
    else:
        percentage_file, complete_file, representative_file = extract_files(input_folder)
        check_file_exists(complete_file)
        check_file_exists(representative_file)
        percentage_df = percentage_to_df(percentage_file)
        len_description_df = len_description_to_df(representative_file)
        summary_df = complete_report_df(complete_file, len_description_df, percentage_df)
        final_individual_dataframe = include_images(input_folder, summary_df)
        summary_to_html(input_folder, final_individual_dataframe, html_template)

    logger.info('DONE')

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logger.exception(e)
        raise
