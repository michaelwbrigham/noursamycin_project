from Bio import AlignIO

import pandas as pd

alignment = AlignIO.read("cladogram-output/unaligned_seq_fasta_musclealign.faa", "fasta")
print("Alignment length %i" % alignment.get_alignment_length())


asn_o = ''
for record in alignment:
    # print("%s - %s" % (record.seq, record.id))

    if 'AsnO' in str(record.id):
        asn_o = record.seq
        print(asn_o)

# [125, 144, 146, 155, 157, 158, 241, 287, 305]

# %%


def get_key_alignmentpositions(key_asno_positions, asno_alignment):
    counter = 0
    position_c = 0
    key_residues = []
    key_alignment_positions = []

    for position in asno_alignment:
        position_c = position_c + 1

        if position != '-':
            counter = counter + 1

            if counter in key_asno_positions:
                key_residues.append(position)
                key_alignment_positions.append(position_c)

    return key_alignment_positions


key_alignment_positions = get_key_alignmentpositions(
    [125, 144, 146, 155, 157, 158, 241, 287, 305, 301], asn_o)

print(key_alignment_positions)
# %%


def return_keypositions(key_positions, query_id, query_alignment):
    counter = 0
    position_c = 0
    key_residues = []
    key_alignment_positions = []

    for position in query_alignment:
        position_c = position_c+1

        if position_c in key_positions:
            key_residues.append(position)
            key_alignment_positions.append(query_alignment.index(position))

    key_residues_febinding = ''.join([key_residues[3], key_residues[4], key_residues[7]])
    key_residues_freecarboxylbinding = ''.join([key_residues[2], key_residues[8]])
    key_residues_aminobinding = key_residues[0]
    key_residues_sidechain = ''.join([key_residues[1], key_residues[5], key_residues[6]])
    key_residues_akg = key_residues[9]
    return [query_id, key_residues_febinding, key_residues_freecarboxylbinding, key_residues_aminobinding, key_residues_sidechain, key_residues_akg]
    # print('FE: {} - F.COOH: {} - F.NH2: {} - F.SC: {} | {}'.format(key_residues_febinding, key_residues_freecarboxylicbinding, key_residues_aminebinding, key_residues_sidechain, query_id))


cols = ['user_name', 'Fe binding residues', 'Carboxyl binding residues',
        'Amino binding residue', 'Side chain binding residues', 'aKD binding residue']
lst = []

for record in alignment:
    lst.append(return_keypositions(key_alignment_positions, record.id, record.seq))

df1 = pd.DataFrame(lst, columns=cols)

df1

query_df = pd.read_excel('raw-data/prunded_mibig_hydroxylase_db_lit.xls')

df = pd.merge(df1, query_df, on=['user_name'])

df

df[df['Fe binding residues'] == 'HEH']
