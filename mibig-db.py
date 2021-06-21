
# mibig-db.py
# reads a .csv file containing URL of MIBiG hit table for the hydroxylase of a cluster (the query for the MIBiG hits)
# for each query, parses through MIBIG hit table and grabs info for each hit (amino acid sequence, compound name, SMILES, organism name, the original query)
# merges all hits from every MIBiG hit table and removes duplicates
# outputs the dataframe containing all this information as a .csv file

# %% codecell
import numpy as np
import pandas as pd
import json
import requests
import lxml
from unicodedata import normalize
import matplotlib.pyplot as plt
import statistics as stats

# %% codecell


def return_mibighit_df(url):
    mibighit_df = pd.read_html(url)[0]

    return mibighit_df
# %% codecell


def append_translation_from_mibighit_row(df):
    translation_list = []

    for index, row in df.iterrows():
        query_mibig_bgc = row['MIBiG Cluster']
        query_protein_id = row['MIBiG Protein']

        hits_bgc_jsondetailed_url = "https://mibig.secondarymetabolites.org/repository/{}/{}.1.json".format(
            query_mibig_bgc, query_mibig_bgc)
        r = requests.get(hits_bgc_jsondetailed_url).json()

        protein_trans = []

        for feature in r['records'][0]['features']:
            if 'CDS' in str(feature):
                if query_protein_id in str(feature):
                    protein_trans.append(feature['qualifiers']['translation'][0])

        if len(protein_trans) == 1:
            translation_list.append(protein_trans[0])
        else:
            translation_list.append('N/A')

    df['aa_sequence'] = translation_list
# %% codecell


def append_compoundname_from_mibighit_row(df):
    compound_list = []

    for index, row in df.iterrows():
        query_mibig_bgc = row['MIBiG Cluster']

        hits_bgc_jsonannotation_url = "https://mibig.secondarymetabolites.org/repository/{}/{}.json".format(
            query_mibig_bgc, query_mibig_bgc)
        r = requests.get(hits_bgc_jsonannotation_url).json()

        if len(r['cluster']['compounds']) > 1:
            name = r['cluster']['compounds'][0]['compound']
            compound_list.append('Multiple including {}'.format(name))
        else:
            name = r['cluster']['compounds'][0]['compound']
            compound_list.append(name)

    df['compound_name'] = compound_list

# %% codecell


def append_compoundsmiles_from_mibighit_row(df):
    smiles_list = []

    for index, row in df.iterrows():
        query_mibig_bgc = row['MIBiG Cluster']

        hits_bgc_jsonannotation_url = "https://mibig.secondarymetabolites.org/repository/{}/{}.json".format(
            query_mibig_bgc, query_mibig_bgc)
        r = requests.get(hits_bgc_jsonannotation_url).json()

        try:
            smiles = r['cluster']['compounds'][0]['chem_struct']
            smiles_list.append(smiles)
        except:
            smiles_list.append('N/A')
    df['compound_smiles'] = smiles_list
# %% codecell


def append_organismname_from_mibighit_row(df):
    organism_list = []

    for index, row in df.iterrows():
        query_mibig_bgc = row['MIBiG Cluster']

        hits_bgc_jsonannotation_url = "https://mibig.secondarymetabolites.org/repository/{}/{}.json".format(
            query_mibig_bgc, query_mibig_bgc)
        r = requests.get(hits_bgc_jsonannotation_url).json()

        try:
            name = r['cluster']['organism_name']
            organism_list.append(name)
        except:
            organism_list.append('N/A')
    df['organism_name'] = organism_list
# %% codecell


def gather_info_from_mibighit_url(url, note):
    df = return_mibighit_df(url)
    append_compoundname_from_mibighit_row(df)
    append_translation_from_mibighit_row(df)
    append_compoundsmiles_from_mibighit_row(df)
    append_organismname_from_mibighit_row(df)
    df.drop(['% ID', '% Coverage', 'BLAST Score', 'E-value'], axis=1)
    df['note'] = note

    return df


# %% codecell
def merge_dfs_drop_duplicates(df_list):
    merged_df = pd.concat(df_list).reset_index()
    merged_df = merged_df.drop_duplicates(subset='MIBiG Protein', keep="first")
    merged_df = merged_df.drop(['index', '% ID', '% Coverage', 'BLAST Score', 'E-value'], axis=1)

    return merged_df


# %% codecell
query_df = pd.read_csv('raw-data/mibig_hit_queries.csv')
query_return_list = []

for index, row in query_df.iterrows():
    query_return_list.append(
        gather_info_from_mibighit_url(row['mibig_hits_url'], row['query_compound'])
    )

query_return_list_noduplicates = merge_dfs_drop_duplicates(query_return_list)

query_return_list_noduplicates.to_csv(
    'processed-data/mibig_hydroxylase_db.csv', index=False)
