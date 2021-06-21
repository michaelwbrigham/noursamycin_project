from math import sqrt
from sklearn import decomposition
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
import os
import tempfile
import subprocess


def check_aaseq_nan(df, i_name, i_seq):
    if df[i_seq].isnull().values.any():
        dropped_protein_name = df[i_name][df[i_seq].isnull()]
        dropped_index = df[df[i_seq].isnull()].index
        dropped_list = list(zip(dropped_index, dropped_protein_name))
        print('Dropped following entries due to no aa seq: {}\t{} entries'.format(
            dropped_list, len(dropped_list)))

    df = df[df[i_seq].notna()]

    return df


def write_fasta_from_df(df, path, file_name, i_seq, i_name, i_cat, i_aa):
    hyd_records = []

    for index, row in df.iterrows():
        record = SeqRecord(Seq(str(row[i_seq])))
        record.id = str(row[i_name] + ' | ' + row[i_cat] + ' | ' + row[i_aa])
        record.description = ''
        hyd_records.append(record)

    SeqIO.write(hyd_records, path+"/"+file_name, "fasta")


def create_muscle_aligment_bashcmd(path, i_file_name):
    file_path = path+"/"+i_file_name
    cmd = MuscleCommandline(input=file_path, out=file_path.replace('.faa', '_musclealign.faa'))
    subprocess.run(str(cmd).split(), stdout=subprocess.DEVNULL)


def aaseq_df_to_align_onehot_df(query_df, i_name='user_name', i_seq='aa_sequence', i_cat='sub_type', i_aa='substrate', path='cladogram-output'):
    print("Column for seq name: {}\nColumn for seq: {}".format(i_name, i_seq))
    print("Number of input entries: {}".format(len(query_df)))
    query_df = check_aaseq_nan(query_df, i_name, i_seq)
    if not os.path.exists(path):
        try:
            os.mkdir(path)
        except OSError:
            print("Creation of the directory {} failed".format(path))
    write_fasta_from_df(query_df, path, "seq_fasta.faa", i_seq, i_name, i_cat, i_aa)
    print('Extracted {} aa sequences from df'.format(len(query_df)))
    create_muscle_aligment_bashcmd(path, "seq_fasta.faa")
    print('Muscle aligned aa sequences with default settings')


def readfile_apply_clad(path):
    query_df = pd.read_excel(path)
    aaseq_df_to_align_onehot_df(query_df)


readfile_apply_clad('raw-data/prunded_mibig_hydroxylase_db_lit.xls')
