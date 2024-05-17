import pandas as pd
import numpy as np
from Bio import SeqIO, Seq
import subprocess
import time
from multiprocess import Pool
import re
import os
from function import return_position_in_repeat_alignment
import warnings
warnings.simplefilter('ignore')

family_list = ['LTR7_HERVH']

s = time.time()
print('start searching binding sites of provirus')

### load data
KZFP_columns = ['KZFP summits chr', 'KZFP summits start', 'KZFP summits end', 'KZFP peak name', 'KZFP score', 'KZFP strand', 'KZFP signal value',
                'KZFP p value', 'KZFP q value', 'KZFP peak chr', 'KZFP peak start', 'KZFP peak end', 'KZFP peak length', 'KZFP accession',
                'KZFP gene symbol', 'KZFP experiment', 'KZFP dataset origin', 'KZFP assay', 'KZFP run type', 'KZFP Cell type', 'KZFP Cell type class', 'KZFP control']

# provirus metadata
Dfam_RM_family = pd.read_csv('../data/TE/provirus_metadata_with_Liftover.csv')
Dfam_RM_family.index = Dfam_RM_family['repeat adjusted subfamily name']
Dfam_RM_family_dict = Dfam_RM_family.to_dict()

# provirus annotation
Dfam_RM = pd.read_csv('../data/TE/provirus_annotation.csv')
Dfam_RM.index = Dfam_RM['repeat name']
Dfam_RM

# KZFP ChIP-seq metadata
KZFP_ChIP_metadata = pd.read_csv('../data/ChIP-seq/KZFP_ChIP-seq_metadata.csv')

# overlap between TE and KZFP peaks
Dfam_RM_overlap_KZFP = pd.read_table('../data/overlap/provirus_KZFP_overlap.bed', header=None)
Dfam_RM_overlap_KZFP.columns = Dfam_RM.columns.tolist() + KZFP_columns + ['overlap length']

# KZFP targets
KZFP_target = pd.read_csv('../data/targets/provirus_targets_for_analysis.csv')


# hg38 genome
hg38_dict = dict()
with open('../data/UCSC/hg38.fa') as f:

    for record in SeqIO.parse(f, 'fasta'):

        hg38_dict[record.id] = record


print(time.time()-s, 'Loading done')


# 全てのマップファイルからアライメントにおけるポジションを記録する
map_dict = dict()
s = time.time()
for i, family in enumerate(family_list):

    with open('../data/MSA/{}.fasta.map'.format(family)) as f:
        
        fr = f.read().split('\n')

        for line in fr:

            split = line.split(', ')

            if '>' in line:

                repeat_name = line.split(' ')[0][1:]
                map_dict[repeat_name] = dict()
        
            elif len(split)== 3 and '#' not in line:

                split = line.split(', ')
                map_dict[repeat_name][int(split[1])] = split[2]


### process data

# extract full-length provirus
Dfam_RM_fil = Dfam_RM[(Dfam_RM['5-LTR name'].isna()==False) & (Dfam_RM['3-LTR name'].isna()==False) & (Dfam_RM['repeat adjusted subfamily name'].isna()==False)]
Dfam_RM_overlap_KZFP_fil = Dfam_RM_overlap_KZFP[Dfam_RM_overlap_KZFP['repeat name'].isin(Dfam_RM_fil['repeat name'])].copy()

# add peak summit position in repeat
def return_position_in_repeat(summit, r_s, r_e, strand):

    if strand == '+':

        return summit - r_s
    
    else:

        return r_e - summit

position_in_repeat_list = list()
for summit, r_s, r_e, strand in Dfam_RM_overlap_KZFP_fil[['KZFP summits start', 'repeat start', 'repeat end', 'repeat strand']].values:

    position = return_position_in_repeat(summit, r_s, r_e, strand)
    position_in_repeat_list.append(position)

Dfam_RM_overlap_KZFP_fil['summit start in repeat'] = position_in_repeat_list

## record KZFP summit position
KZFP_summits_in_repeat_dict = dict()
for i, (family, repeat_name, peak_name, strand, position) in enumerate(Dfam_RM_overlap_KZFP_fil[['repeat family name', 'repeat name', 'KZFP peak name', 'repeat strand', 'summit start in repeat']].values):

    KZFP_summits_in_repeat_dict[(family, peak_name)] = [repeat_name, strand, position]


# summarize KZFP targets
target_dict = dict()
KZFP_target_nodup = KZFP_target[KZFP_target['repeat family name'].isin(family_list)]
KZFP_target_nodup = KZFP_target_nodup[KZFP_target_nodup[['KZFP gene symbol', 'repeat family name']].duplicated()==False]
for KZFP, family in KZFP_target_nodup[['KZFP gene symbol', 'repeat family name']].values:

    if KZFP not in target_dict:

        target_dict[KZFP] = [family]

    else:

        target_dict[KZFP] += [family]


# KRAB-ZFP
position_in_repeat_alignment_list = list()
for repeat_name, position in Dfam_RM_overlap_KZFP_fil[['repeat name', 'summit start in repeat']].values:

    position = return_position_in_repeat_alignment(repeat_name, position, map_dict)
    position_in_repeat_alignment_list.append(position)

Dfam_RM_overlap_KZFP_fil['summit position in repeat alignment'] = position_in_repeat_alignment_list


print(time.time()-s, 'processing done')


### output KZFP peak sequence


# obtain overlap data
Dfam_KZFP_familly_dict = dict()
for i, (KZFP, family_list) in enumerate(target_dict.items()):

    KZFP_df = Dfam_RM_overlap_KZFP_fil[Dfam_RM_overlap_KZFP_fil['KZFP gene symbol']==KZFP]
    
    for family in family_list:

        df = KZFP_df[KZFP_df['repeat family name']==family]
        Dfam_KZFP_familly_dict[(family, KZFP)] = df


# obtain peak sequence
Dfam_KZFP_familly_sequence_dict = dict()
for i, (key, df) in enumerate(Dfam_KZFP_familly_dict.items()):

    family, KZFP = key
    df = Dfam_KZFP_familly_dict[(family, KZFP)]
    df_nodup = df[df['KZFP peak name'].duplicated()==False]

    record_list = list()
    for chr, start, strand, name in df_nodup[['KZFP summits chr', 'KZFP summits start', 'repeat strand', 'KZFP peak name']].values:

        record = hg38_dict[chr][start-50:start+50]
        record.id = name
        record.description = name
        record.name = name

        if strand == '-':

            record.seq = record.seq.reverse_complement()

        record_list.append(record)
    
    Dfam_KZFP_familly_sequence_dict[(family, KZFP)] = record_list


# output peak sequence
for i, (key, df) in enumerate(Dfam_KZFP_familly_sequence_dict.items()):

    family, KZFP = key

    with open('../data/motif/peak_sequences/{}_{}_sequence.fasta'.format(family, KZFP), 'w') as f:

        for record in df:

            SeqIO.write(record, f, 'fasta')


print(time.time()-s, 'obtaining peak sequences done')



### scan KZFP motif using fimo

motif_PATH_set = set(os.listdir('../data/motif/raw_motif'))
s = time.time()

for i, (KZFP, exp) in enumerate(KZFP_ChIP_metadata[['KZFP gene symbol', 'KZFP experiment']].value_counts().index[2:]):

    if KZFP in target_dict.keys():
        target_list = target_dict[KZFP]
    
    else:
        continue

    motif = '../data/motif/raw_motif/{}_{}_motif_sig_500_231120.meme'.format(KZFP, exp)

    for k, family in enumerate(target_list):

        fasta = '../data/motif/peak_sequences/{}_{}_sequence.fasta'.format(family, KZFP)
        oc = '../data/motif/fimo/{}_{}_{}'.format(family, KZFP, exp)
        option = '--thresh 0.05 --no-qvalue'

        cmd = 'fimo --oc {} {} {}'.format(oc, motif, fasta)
        cmd = 'fimo --oc {} {} {} {}'.format(oc, option, motif, fasta)

        output_str = subprocess.run(cmd, shell=True, capture_output=True, text=True)

print(time.time()-s, 'scanning done')



### process fimo output files

KZFP_experiment_metadata = KZFP_ChIP_metadata[KZFP_ChIP_metadata[['KZFP gene symbol', 'KZFP experiment']].duplicated()==False]
KZFP_experiment_metadata_dict = dict()

for KZFP in pd.unique(KZFP_experiment_metadata['KZFP gene symbol']):

    df = KZFP_experiment_metadata[KZFP_experiment_metadata['KZFP gene symbol']==KZFP]
    KZFP_experiment_metadata_dict[KZFP] = df


# obtain fimo output files
fimo_dict = dict()
for i, (family, KZFP) in enumerate(KZFP_target[['repeat family name', 'KZFP gene symbol']].value_counts().index):

    exp_list = KZFP_experiment_metadata_dict[KZFP]
    for exp in exp_list['KZFP experiment']:

        try:
            fimo = pd.read_table('../data/motif/fimo/{}_{}_{}/fimo.tsv'.format(family, KZFP, exp))
            fimo['motif center'] = abs((fimo['stop'] + fimo['start']) // 2 - 50)
            fimo['log10 p-value'] = -fimo['p-value'].apply(np.log10)
            fimo = fimo[fimo['sequence_name'].isna()==False]
 
            fimo_dict[(family, KZFP, exp)] = fimo
        
        except:
            continue


# add motif position in repeat alignment
def return_motif_position_in_repeat(position, start_or_end):

    return position + start_or_end - 50


fimo_modified_dict = dict()
for i, (family, KZFP) in enumerate(KZFP_target[['repeat family name', 'KZFP gene symbol']].value_counts().index):

    exp_list = KZFP_experiment_metadata_dict[KZFP]
    for exp in exp_list['KZFP experiment']:

        try:
            fimo = fimo_dict[(family, KZFP, exp)]
        
        except:
            continue

        motif_position_list = [[], []]
        motif_alignment_position_list = [[], []]

        motif_position_list = [[], [], [], []]
        
        for k, (peak_name, start, end) in enumerate(fimo[['sequence_name', 'start', 'stop']].values):

            try:
                
                repeat_name, strand, position = KZFP_summits_in_repeat_dict[(family, peak_name)]

            except:

                for i in range(4):

                    motif_position_list[i].append(np.nan)

                continue

            start_in_r = position + start - 50
            end_in_r = position + end - 50

            motif_position_list[0].append(start_in_r)
            motif_position_list[1].append(end_in_r)

            try:

                motif_position_list[2].append(return_position_in_repeat_alignment(repeat_name, start_in_r, map_dict))

            except:

                motif_position_list[2].append(np.nan)
                #print('start', start_in_r, end_in_r)

            try:
                
                motif_position_list[3].append(return_position_in_repeat_alignment(repeat_name, end_in_r, map_dict))
            
            except:
                motif_position_list[3].append(np.nan)
                #print('end', start_in_r, end_in_r)


        fimo['motif start in repeat'] = motif_position_list[0]
        fimo['motif end in repeat'] = motif_position_list[1]
        fimo['motif start in repeat alignment'] = motif_position_list[2]
        fimo['motif end in repeat alignment'] = motif_position_list[3]

        fimo_modified_dict[(family, KZFP, exp)] = fimo.copy()


# output fimo results
for i, (key, df) in enumerate(fimo_modified_dict.items()):

    name = '{}_{}_{}_fimo_results.csv'.format(key[0], key[1], key[2])
    df.to_csv('../data/motif/fimo/{}'.format(name), index=False)

print(time.time()-s, 'processing fimo outputs done')


### summarize fimo results

motif_metadata_df = pd.DataFrame(fimo_modified_dict.keys(), columns=['repeat family name', 'KZFP gene symbol', 'KZFP experiment'])

# obtain fimo results
fimo_modified_dict = dict()
for i, key in enumerate(motif_metadata_df[['repeat family name', 'KZFP gene symbol', 'KZFP experiment']].values):

    name = '{}_{}_{}_fimo_results.csv'.format(key[0], key[1], key[2])

    df = pd.read_csv('../data/motif/fimo/{}'.format(name))
    df['log10 p-value'] = -df['p-value'].apply(np.log10)
    df_fil = df[df['log10 p-value']>4]
        
    fimo_modified_dict[tuple(key)] = df_fil


# obtain the number of input sequence (for calculating discovery rates)
seuqnece_num_dict = dict()
for i, key in enumerate(motif_metadata_df[['repeat family name', 'KZFP gene symbol', 'KZFP experiment']].values):

    family, KZFP, exp = key
    count = 0
    with open('../data/motif/peak_sequences/{}_{}_sequence.fasta'.format(family, KZFP)) as f:

        for record in SeqIO.parse(f, 'fasta'):
            count += 1
    
    seuqnece_num_dict[tuple(key)] = count


# calculate discovery rate
fimo_groupby_dict = dict()
fimo_groupby_fil_dict = dict()
for i, (key, fimo) in enumerate(fimo_modified_dict.items()):
    
    groupby_count = fimo[['motif start in repeat alignment', 'motif end in repeat alignment', 'strand']].value_counts().rename('overlap count').to_frame()
    groupby_pvalue = fimo[['motif start in repeat alignment', 'motif end in repeat alignment', 'strand', 'log10 p-value']].groupby(by=['motif start in repeat alignment', 'motif end in repeat alignment', 'strand']).median()
    groupby_center = fimo[['motif start in repeat alignment', 'motif end in repeat alignment', 'strand', 'motif center']].groupby(by=['motif start in repeat alignment', 'motif end in repeat alignment', 'strand']).mean()

    try:
        
        groupby = pd.concat([groupby_count, groupby_pvalue, groupby_center], axis=1).sort_values(by='overlap count', ascending=False)
        groupby['relative overlap count'] = groupby['overlap count'] / groupby['overlap count'].iloc[0]
        groupby['discovery rate'] = groupby['overlap count'] / seuqnece_num_dict[tuple(key)] * 100
        fimo_groupby_dict[key] = groupby

    except:
        pass


# summarize results
data_list = list()
for i, (family, KZFP, exp) in enumerate(motif_metadata_df[['repeat family name', 'KZFP gene symbol', 'KZFP experiment']].values):

    try:
        groupby = fimo_groupby_dict[(family, KZFP, exp)]
        groupby_fil = groupby[(groupby['discovery rate']>50)].sort_values(by='log10 p-value', ascending=False)
    
    except:
        #print(family, KZFP, exp, 'no data')
        continue

    if len(groupby_fil) == 0:
        #print(family, KZFP, exp, 'no candiate')
        continue

    overlap, pvalue, distance, relative_count, proportion = groupby_fil.iloc[0]
    start, end, strand = groupby_fil.index[0]
    center = (int(start) + int(end)) // 2

    data_list.append([family, KZFP, exp, start, end, center, strand, overlap, proportion, distance, pvalue])

motif_metadata_modified_df = pd.DataFrame(data_list, columns=['repeat family name', 'KZFP gene symbol', 'KZFP experiment', 'motif start', 'motif end', 'motif center', 'motif strand', 'peak overlap motif', 'discovery proportion', 'distance to peak summit', 'log10 p-value'])


# add E-value of motif
regex_Evalue = re.compile(r'E= (.+)')
regex_Evalue_float = re.compile(r'E= (.+)e-(.+)')
E_value_list = [[], []]
for i, (KZFP, exp) in enumerate(motif_metadata_modified_df[['KZFP gene symbol', 'KZFP experiment']].values):

    with open('../data/motif/raw_motif/{}_{}_motif_sig_500_231120.meme'.format(KZFP, exp)) as f:

        fr = f.read().split('\n')

        for line in fr:

            if 'letter-probability' in line:

                eff = regex_Evalue_float.search(line).group(1)
                log = regex_Evalue_float.search(line).group(2)
                E_value_list[0].append(float(eff))
                E_value_list[1].append(int(log))


motif_metadata_modified_df['E-value significant figure'] = E_value_list[0]
motif_metadata_modified_df['E-value log10'] = E_value_list[1]

motif_metadata_modified_df.to_csv('../data/motif/provirus_KZFP_binding_sites.csv', index=False)

print(time.time()-s, 'finish searching binding sites of provirus')