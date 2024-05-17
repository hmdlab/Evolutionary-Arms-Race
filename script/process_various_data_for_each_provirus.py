import pandas as pd
import numpy as np
from Bio import SeqIO
import time
from Bio import Phylo
import baltic as bt
from function import return_position_in_repeat_alignment, return_yaxis
import warnings
warnings.simplefilter('ignore')

family_list = ['LTR7_HERVH', 'THE1_THE1-int']

s = time.time()
print('start processing provirus family data')

### load data
TRIM28_columns = ['TRIM28 summits chr', 'TRIM28 summits start', 'TRIM28 summits end', 'TRIM28 peak name', 'TRIM28 score', 'TRIM28 strand', 'TRIM28 signal value',
                  'TRIM28 p value', 'TRIM28 q value', 'TRIM28 peak chr', 'TRIM28 peak start', 'TRIM28 peak end', 'TRIM28 peak length', 'TRIM28 accession',
                  'TRIM28 experiment', 'TRIM28 assay', 'TRIM28 run type', 'TRIM28 Cell type', 'TRIM28 Cell type class', 'TRIM28 control']

KZFP_columns = ['KZFP summits chr', 'KZFP summits start', 'KZFP summits end', 'KZFP peak name', 'KZFP score', 'KZFP strand', 'KZFP signal value',
                'KZFP p value', 'KZFP q value', 'KZFP peak chr', 'KZFP peak start', 'KZFP peak end', 'KZFP peak length', 'KZFP accession',
                'KZFP gene symbol', 'KZFP experiment', 'KZFP dataset origin', 'KZFP assay', 'KZFP run type', 'KZFP Cell type', 'KZFP Cell type class', 'KZFP control']

TF_columns = ['TF summits chr', 'TF summits start', 'TF summits end', 'TF peak name', 'TF score', 'TF strand', 'TF signal value',
              'TF p value', 'TF q value', 'TF peak chr', 'TF peak start', 'TF peak end', 'TF peak length', 'TF accession',
              'TF experiment', 'TF gene symbol', 'TF assay', 'TF run type', 'TF Cell type', 'TF Cell type class', 'TF control']

# repeat familyのデータを取得する
Dfam_RM_family = pd.read_csv('../data/TE/provirus_metadata_with_Liftover.csv')
Dfam_RM_family.index = Dfam_RM_family['repeat adjusted subfamily name']

# repeatのデータを取得する
Dfam_RM = pd.read_csv('../data/TE/provirus_annotation.csv')
Dfam_RM.index = Dfam_RM['repeat name']

# KZFP peaks
Dfam_RM_overlap_KZFP = pd.read_table('../data/overlap/provirus_KZFP_overlap.bed', header=None)
Dfam_RM_overlap_KZFP.columns = Dfam_RM.columns.tolist() + KZFP_columns + ['overlap length']

# TRIM28 peaks
Dfam_RM_overlap_TRIM28 = pd.read_table('../data/overlap/provirus_TRIM28_overlap.bed', header=None)
Dfam_RM_overlap_TRIM28.columns = Dfam_RM.columns.tolist() + TRIM28_columns + ['overlap length']

# TF peaks
Dfam_RM_overlap_TF = pd.read_table('../data/overlap/provirus_TF_overlap.bed', header=None)
Dfam_RM_overlap_TF.columns = Dfam_RM.columns.tolist() + TF_columns + ['overlap length']

# Liftover
LiftOver_df_summary = pd.read_csv('../data/LiftOver/provirus_Liftover_summary.csv', index_col=0)
LiftOver_df = pd.read_csv('../data/LiftOver/provirus_Liftover.csv', index_col=0)
LiftOver_df.index = LiftOver_df['repeat name']

# map files
map_dict = dict()
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


# tree and branch length
tree_dict = dict()
branch_length_dict = dict()
for family in family_list:

    tree = Phylo.read('../data/phylogenetic tree/{}.contree'.format(family), format='newick')
    tree_dict[family] = tree

    for node in tree.get_terminals():

        branch_length_dict[node.name] = node

print(time.time()-s, 'Loading done')


# process data

Dfam_RM_fil = Dfam_RM[(Dfam_RM['5-LTR name'].isna()==False) & (Dfam_RM['3-LTR name'].isna()==False) & (Dfam_RM['repeat adjusted subfamily name'].isna()==False)]
Dfam_RM_overlap_TRIM28_fil = Dfam_RM_overlap_TRIM28[Dfam_RM_overlap_TRIM28['repeat name'].isin(Dfam_RM_fil['repeat name'])].copy()
Dfam_RM_overlap_KZFP_fil = Dfam_RM_overlap_KZFP[Dfam_RM_overlap_KZFP['repeat name'].isin(Dfam_RM_fil['repeat name'])].copy()
Dfam_RM_overlap_TF_fil = Dfam_RM_overlap_TF[Dfam_RM_overlap_TF['repeat name'].isin(Dfam_RM_fil['repeat name'])].copy()
LiftOver_df_fil = LiftOver_df[LiftOver_df['repeat name'].isin(Dfam_RM_fil['repeat name'])].copy()


# add repeat subfamily name and repeat name
Dfam_RM_fil_dict = Dfam_RM_fil.to_dict()

Dfam_RM_overlap_TRIM28_fil['repeat adjusted subfamily name'] = Dfam_RM_overlap_TRIM28_fil['repeat name'].apply(lambda x:Dfam_RM_fil_dict['repeat adjusted subfamily name'][x])
Dfam_RM_overlap_KZFP_fil['repeat adjusted subfamily name'] = Dfam_RM_overlap_KZFP_fil['repeat name'].apply(lambda x:Dfam_RM_fil_dict['repeat adjusted subfamily name'][x])
Dfam_RM_overlap_TF_fil['repeat adjusted subfamily name'] = Dfam_RM_overlap_TF_fil['repeat name'].apply(lambda x:Dfam_RM_fil_dict['repeat adjusted subfamily name'][x])
LiftOver_df_fil['repeat adjusted subfamily name'] = LiftOver_df_fil['repeat name'].apply(lambda x:Dfam_RM_fil_dict['repeat adjusted subfamily name'][x])
LiftOver_df_fil['repeat adjusted name'] = LiftOver_df_fil['repeat name'].apply(lambda x:Dfam_RM_fil_dict['repeat adjusted name'][x])

# add repeat summit position in repeat alignment
def return_position_in_repeat(summit, r_s, r_e, strand):

    if strand == '+':
        return summit - r_s
    
    else:
        return r_e - summit

def return_region_in_repeat(summit, LTR5_s, LTR5_e, LTR3_s, LTR3_e, strand):
            
    if LTR5_s <= summit and summit <= LTR5_e : 
        return '5-LTR'
            
    elif LTR3_s <= summit and summit <= LTR3_e:
        return '3-LTR'

    else:

        if (LTR5_e <= summit and summit <= LTR3_s) or (LTR3_e <= summit and summit <= LTR5_s):
            return 'Int'

        else:
            print(summit, LTR5_s, LTR5_e, LTR3_s, LTR3_e, strand)

# TRIM28
position_in_repeat_list = list()
region_in_repeat_list = list()
position_in_repeat_alignment_list = list()
for summit, r_s, r_e, strand, repeat_name, LTR5_s, LTR5_e, LTR3_s, LTR3_e in Dfam_RM_overlap_TRIM28_fil[['TRIM28 summits start', 'repeat start', 'repeat end', 'repeat strand', 'repeat name', '5-LTR start', '5-LTR end', '3-LTR start', '3-LTR end']].values:

    position = return_position_in_repeat(summit, r_s, r_e, strand)
    region = return_region_in_repeat(summit, LTR5_s, LTR5_e, LTR3_s, LTR3_e, strand)
    position_MSA = return_position_in_repeat_alignment(repeat_name, position, map_dict)

    position_in_repeat_list.append(position)
    region_in_repeat_list.append(region)
    position_in_repeat_alignment_list.append(position_MSA)

Dfam_RM_overlap_TRIM28_fil['summit start in repeat'] = position_in_repeat_list
Dfam_RM_overlap_TRIM28_fil['summit region in repeat'] = region_in_repeat_list
Dfam_RM_overlap_TRIM28_fil['summit position in repeat alignment'] = position_in_repeat_alignment_list


# KRAB-ZFP
position_in_repeat_list = list()
region_in_repeat_list = list()
position_in_repeat_alignment_list = list()
for summit, r_s, r_e, strand, repeat_name, LTR5_s, LTR5_e, LTR3_s, LTR3_e in Dfam_RM_overlap_KZFP_fil[['KZFP summits start', 'repeat start', 'repeat end', 'repeat strand', 'repeat name', '5-LTR start', '5-LTR end', '3-LTR start', '3-LTR end']].values:

    position = return_position_in_repeat(summit, r_s, r_e, strand)
    region = return_region_in_repeat(summit, LTR5_s, LTR5_e, LTR3_s, LTR3_e, strand)
    position_MSA = return_position_in_repeat_alignment(repeat_name, position, map_dict)

    position_in_repeat_list.append(position)
    region_in_repeat_list.append(region)
    position_in_repeat_alignment_list.append(position_MSA)

Dfam_RM_overlap_KZFP_fil['summit start in repeat'] = position_in_repeat_list
Dfam_RM_overlap_KZFP_fil['summit region in repeat'] = region_in_repeat_list
Dfam_RM_overlap_KZFP_fil['summit position in repeat alignment'] = position_in_repeat_alignment_list

# TF
position_in_repeat_list = list()
region_in_repeat_list = list()
position_in_repeat_alignment_list = list()
for summit, r_s, r_e, strand, repeat_name, LTR5_s, LTR5_e, LTR3_s, LTR3_e in Dfam_RM_overlap_TF_fil[['TF summits start', 'repeat start', 'repeat end', 'repeat strand', 'repeat name', '5-LTR start', '5-LTR end', '3-LTR start', '3-LTR end']].values:

    position = return_position_in_repeat(summit, r_s, r_e, strand)
    region = return_region_in_repeat(summit, LTR5_s, LTR5_e, LTR3_s, LTR3_e, strand)
    position_MSA = return_position_in_repeat_alignment(repeat_name, position, map_dict)

    position_in_repeat_list.append(position)
    region_in_repeat_list.append(region)
    position_in_repeat_alignment_list.append(position_MSA)

Dfam_RM_overlap_TF_fil['summit start in repeat'] = position_in_repeat_list
Dfam_RM_overlap_TF_fil['summit region in repeat'] = region_in_repeat_list
Dfam_RM_overlap_TF_fil['summit position in repeat alignment'] = position_in_repeat_alignment_list

# Liftover
branch_list = ['Vertebrata', 'Mammalia', 'Theria', 'Eutheria', 'Boreoeutheria', 'Euarchontoglires', 'Primatomorpha', 
               'Primates', 'Haplorrhini', 'Simiiformes', 'Catarrhini', 'Hominoidea', 'Hominidae', 'Homininae', 'Hominini', 'Homo sapiens']

LiftOver_crosstab = pd.crosstab(Dfam_RM_fil['repeat adjusted subfamily name'], LiftOver_df_fil['branch'], normalize='index') * 100

for branch in branch_list:
    
    if branch not in LiftOver_crosstab.columns:

        LiftOver_crosstab[branch] = 0

# add branch length for root
LiftOver_df_fil['branch length'] = LiftOver_df_fil['repeat adjusted name'].apply(lambda x:branch_length_dict[x].branch_length if x in branch_length_dict.keys() else np.nan)
Dfam_RM_fil['branch length'] = Dfam_RM_fil['repeat adjusted name'].apply(lambda x:branch_length_dict[x].branch_length if x in branch_length_dict.keys() else np.nan)

branch_length_groupby = LiftOver_df_fil[['repeat adjusted subfamily name', 'branch length']].groupby(by='repeat adjusted subfamily name').mean()
Dfam_RM_Liftover_metadata = pd.concat([Dfam_RM_family, LiftOver_crosstab, branch_length_groupby['branch length']], axis=1)


print(time.time()-s, 'processing done')


### process and output data for each TE families

# peaks
for i, family in enumerate(family_list):

    KZFP = Dfam_RM_overlap_KZFP_fil[Dfam_RM_overlap_KZFP_fil['repeat family name']==family]
    TRIM28 = Dfam_RM_overlap_TRIM28_fil[Dfam_RM_overlap_TRIM28_fil['repeat family name']==family]
    TF = Dfam_RM_overlap_TF_fil[Dfam_RM_overlap_TF_fil['repeat family name']==family]

    KZFP.to_csv('../data/overlap/{}_KZFP.csv'.format(family), index=False)
    TRIM28.to_csv('../data/overlap/{}_TRIM28.csv'.format(family), index=False)
    TF.to_csv('../data/overlap/{}_TF.csv'.format(family), index=False)


# Liftover
for i, family in enumerate(family_list):

    subfamily_list = Dfam_RM_family[Dfam_RM_family['repeat family name']==family]
    Dfam = Dfam_RM_fil[(Dfam_RM_fil['repeat adjusted subfamily name'].isin(subfamily_list.index)) & (Dfam_RM_fil['repeat name'].isin(map_dict.keys()))]

    Lift = LiftOver_crosstab.loc[subfamily_list.index]
    summary = LiftOver_df_summary.loc[Dfam['repeat name']]

    Lift.to_csv('../data/Liftover/{}_Liftover_rate.csv'.format(family))
    summary.to_csv('../data/Liftover/{}_Liftover.csv'.format(family))



# tree
for family in family_list:

    Lift = Dfam_RM_Liftover_metadata[Dfam_RM_Liftover_metadata['repeat family name']==family].sort_values(by=branch_list, ascending=False)
    branch, subfamily = Lift.iloc[0][['branch', 'repeat adjusted subfamily name']]

    if family == 'THE1_THE1-int':
        subfamily = 'THE1_THE1-int_1'
        
    Dfam = LiftOver_df_fil[(LiftOver_df_fil['repeat adjusted subfamily name']==subfamily) & (LiftOver_df_fil['branch']==branch)].sort_values(by='branch length')

    outgroup = Dfam.iloc[-1]['repeat adjusted name']
    if family == 'LTR7_HERVH':

        outgroup = Dfam.iloc[-2]['repeat adjusted name']

    tree = tree_dict[family]
    tree.root_with_outgroup(outgroup)

    Phylo.write(tree, '../data/phylogenetic tree/{}_reroot.contree'.format(family), format='newick')


# TE annotation
for family in family_list:

    Dfam = Dfam_RM_fil[Dfam_RM_fil['repeat family name']==family]
    tree = bt.loadNewick('../data/phylogenetic tree/{}_reroot.contree'.format(family))
    order = return_yaxis(tree).sort_values('y')
    order_fil = order[order['repeat adjusted name'].isin(Dfam_RM['repeat adjusted name'])]

    # 調整する
    Dfam.index = Dfam['repeat adjusted name']
    Dfam = Dfam.loc[order_fil['repeat adjusted name']]
    Dfam['branch length'] = Dfam['repeat adjusted name'].apply(lambda x:branch_length_dict[x].branch_length)

    # 出力する
    Dfam.to_csv('../data/TE/{}.annotation.csv'.format(family), index=False)


# MSA
MSA_dict = dict()

for family in family_list:

    seq_list = list()
    name_list = list()

    with open('../data/MSA/{}_alignment_to_consensus_merge.fasta'.format(family)) as f:

        for record in SeqIO.parse(f, 'fasta'):

            seq_list.append(list(str(record.seq)))
            name_list.append(record.id)

    MSA = pd.DataFrame(seq_list, index=name_list)
    MSA_dict[family] = MSA

    MSA.to_csv('../data/MSA/{}_MSA.csv'.format(family))


# conserved level in each position 
for i, family in enumerate(family_list):

    MSA = MSA_dict[family]
    df = Dfam_RM[Dfam_RM['repeat family name']==family]
    df = df[df.index.isin(MSA.index)]

    for subfamily in pd.unique(df['repeat adjusted subfamily name']):

        subfamily_df = df[df['repeat adjusted subfamily name']==subfamily]
        subfamily_MSA = MSA.loc[subfamily_df['repeat name']]

        conserved_list = [pd.DataFrame(index=['A', 'C', 'G', 'T', '-'])]
        for col in subfamily_MSA.columns:

            count = subfamily_MSA[col].apply(lambda x:x.upper()).value_counts(normalize=True).rename(col)
            conserved_list.append(count)
        
        conserved_df = pd.concat(conserved_list, axis=1).fillna(0)
        conserved_df.to_csv('../data/MSA/conserved/{}_conserved.csv'.format(subfamily))


print(time.time()-s, 'All data outputed')
