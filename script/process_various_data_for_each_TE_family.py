import pandas as pd
import numpy as np
from Bio import SeqIO
import time
from Bio import Phylo
import baltic as bt
from function import return_position_in_repeat_alignment, return_yaxis
import warnings
warnings.simplefilter('ignore')

s = time.time()
print('start processing TE family data')

family_list = ['LTR5', 'L1P_5end']

### load data

TRIM28_columns = ['TRIM28 summits chr', 'TRIM28 summits start', 'TRIM28 summits end', 'TRIM28 peak name', 'TRIM28 score', 'TRIM28 strand', 'TRIM28 signal value',
                  'TRIM28 p value', 'TRIM28 q value', 'TRIM28 peak chr', 'TRIM28 peak start', 'TRIM28 peak end', 'TRIM28 peak length', 'TRIM28 accession',
                  'TRIM28 experiment', 'TRIM28 assay', 'TRIM28 run type', 'TRIM28 Cell type', 'TRIM28 Cell type class', 'TRIM28 control']

KZFP_columns = ['KZFP summits chr', 'KZFP summits start', 'KZFP summits end', 'KZFP peak name', 'KZFP score', 'KZFP strand', 'KZFP signal value',
                'KZFP p value', 'KZFP q value', 'KZFP peak chr', 'KZFP peak start', 'KZFP peak end', 'KZFP peak length', 'KZFP accession',
                'KZFP gene symbol', 'KZFP experiment', 'KZFP dataset origin', 'KZFP assay', 'KZFP run type', 'KZFP Cell type', 'KZFP Cell type class', 'KZFP control']

# TE metadata
Dfam_RM_family = pd.read_csv('../data/TE/TE_metadata_with_Liftover.csv')
Dfam_RM_family.index = Dfam_RM_family['repeat subfamily name']

# TE annotation
Dfam_RM = pd.read_csv('../data/TE/TE_annotation.csv')
Dfam_RM.index = Dfam_RM['repeat name']

# overlap with KZFP peaks
Dfam_RM_overlap_KZFP = pd.read_table('../data/overlap/TE_KZFP_overlap.bed', header=None)
Dfam_RM_overlap_KZFP.columns = Dfam_RM.columns.tolist() + KZFP_columns + ['overlap length']

# overlap with TRIM28 peaks
Dfam_RM_overlap_TRIM28 = pd.read_table('../data/overlap/TE_TRIM28_overlap.bed', header=None)
Dfam_RM_overlap_TRIM28.columns = Dfam_RM.columns.tolist() + TRIM28_columns + ['overlap length']

# Liftover
LiftOver_df = pd.read_csv('../data/Liftover/TE_Liftover.csv')
LiftOver_df.index = LiftOver_df['repeat name']

LiftOver_summary = pd.read_csv('../data/Liftover/TE_Liftover_summary.csv', index_col=0)

# screening for results
screening_result = pd.read_csv('../data/screening/screening_TE_evasion_candidates.csv')

# map files
map_dict = dict()
s = time.time()
for i, family in enumerate(family_list):

    with open('../data/MSA/{}.fasta.map'.format(family)) as f:
        
        fr = f.read().split('\n')

        for line in fr:

            split = line.split(', ')

            if '>' in line:

                repeat_name = line[1:]
                map_dict[repeat_name] = dict()
        
            elif len(split)== 3 and '#' not in line:

                split = line.split(', ')
                map_dict[repeat_name][int(split[1])] = split[2]



# tree and branch length
tree_dict = dict()
branch_length_dict = dict()
for family in family_list:

    try:

        tree = Phylo.read('../data/phylogenetic tree/{}.contree'.format(family), format='newick')
        tree_dict[family] = tree
    
    except:

        #del tree
        trees = Phylo.parse('../data/phylogenetic tree/{}.treefile'.format(family), format='newick')
        print(family)

        for i, tree in enumerate(trees):

            if len(tree.get_terminals())>=2:
                tree_dict[family] = tree
                print('done')
                break
    
    for node in tree.get_terminals():

        branch_length_dict[node.name] = node


print(time.time()-s, 'Loading done')



### process all data

# filter 
# remove TE subfamilies that have low copies (<10) and TE family without subfamily classification
Dfam_RM_family_fil = Dfam_RM_family[Dfam_RM_family['repeat copy count after filtering']>=10]
family_num = Dfam_RM_family_fil['repeat family name'].value_counts()
Dfam_RM_family_fil = Dfam_RM_family_fil[Dfam_RM_family_fil['repeat family name'].isin(family_num[family_num>=2].index)]

# extract TE 
TE_class_list = ['DNA', 'ERV/LTR', 'ERV/Int', 'LINE', 'SINE', 'Retroposon']
Dfam_RM_family_fil = Dfam_RM_family_fil[Dfam_RM_family_fil['repeat class'].isin(TE_class_list)]

# extract full-length TE
Dfam_RM_fil = Dfam_RM.copy()
Dfam_RM_fil = Dfam_RM_fil[Dfam_RM['repeat filtering']]

Dfam_RM_overlap_TRIM28_fil = Dfam_RM_overlap_TRIM28[Dfam_RM_overlap_TRIM28['repeat name'].isin(set(Dfam_RM_fil['repeat name']))].copy()
Dfam_RM_overlap_KZFP_fil = Dfam_RM_overlap_KZFP[Dfam_RM_overlap_KZFP['repeat name'].isin(set(Dfam_RM_fil['repeat name']))].copy()
LiftOver_df_fil = LiftOver_df.loc[Dfam_RM_fil.index]
LiftOver_summary_fil = LiftOver_summary.loc[Dfam_RM_fil.index]

# add evolutionary age
branch_dict = {'Vertebrata':530, 'Tetrapoda':400, 'Amniota':312.1, 'Mammalia':177, 'Theria':159, 'Eutheria':105.0, 'Boreoeutheria':96, 'Euarchontoglires': 90, 'Primatomorpha':76, 
               'Primates':74, 'Haplorrhini':63, 'Simiiformes':43.2, 'Catarrhini':29.4, 'Hominoidea':20.2, 'Hominidae':15.7, 'Homininae':9.1, 'Hominini':6.7, 'Homo sapiens':0, None: np.nan}

Dfam_RM_fil['branch'] = LiftOver_df_fil['branch']
Dfam_RM_fil['Age'] = Dfam_RM_fil['branch'].apply(lambda x:branch_dict[x])

# add repeat family name
Dfam_RM_family_dict = Dfam_RM_family.to_dict()
Dfam_RM_fil_dict = Dfam_RM_fil.to_dict()

Dfam_RM_overlap_TRIM28_fil['repeat family name'] = Dfam_RM_overlap_TRIM28_fil['repeat name'].apply(lambda x:Dfam_RM_fil_dict['repeat family name'][x])
Dfam_RM_overlap_KZFP_fil['repeat family name'] = Dfam_RM_overlap_KZFP_fil['repeat name'].apply(lambda x:Dfam_RM_fil_dict['repeat family name'][x])

# add peak summit position in repeat alignment
def return_position_in_repeat(summit, r_s, r_e, strand):

    if strand == '+':

        return summit - r_s
    
    else:

        return r_e - summit

# TRIM28
position_in_repeat_list = list()
position_in_repeat_alignment_list = list()
for summit, r_s, r_e, strand, repeat_name in Dfam_RM_overlap_TRIM28_fil[['TRIM28 summits start', 'repeat start', 'repeat end', 'repeat strand', 'repeat name']].values:
    position = return_position_in_repeat(summit, r_s, r_e, strand)
    position_in_MSA = return_position_in_repeat_alignment(repeat_name, position, map_dict)
    position_in_repeat_list.append(position)
    position_in_repeat_alignment_list.append(position_in_MSA)

Dfam_RM_overlap_TRIM28_fil['summit start in repeat'] = position_in_repeat_list
Dfam_RM_overlap_TRIM28_fil['summit position in repeat alignment'] = position_in_repeat_alignment_list

# KZFP
position_in_repeat_list = list()
position_in_repeat_alignment_list = list()
for summit, r_s, r_e, strand, repeat_name in Dfam_RM_overlap_KZFP_fil[['KZFP summits start', 'repeat start', 'repeat end', 'repeat strand', 'repeat name']].values:
    position = return_position_in_repeat(summit, r_s, r_e, strand)
    position_in_MSA = return_position_in_repeat_alignment(repeat_name, position, map_dict)
    position_in_repeat_list.append(position)
    position_in_repeat_alignment_list.append(position_in_MSA)

Dfam_RM_overlap_KZFP_fil['summit start in repeat'] = position_in_repeat_list
Dfam_RM_overlap_KZFP_fil['summit position in repeat alignment'] = position_in_repeat_alignment_list

# Liftover
branch_list = ['Vertebrata', 'Mammalia', 'Theria', 'Eutheria', 'Boreoeutheria', 'Euarchontoglires', 'Primatomorpha', 
               'Primates', 'Haplorrhini', 'Simiiformes', 'Catarrhini', 'Hominoidea', 'Hominidae', 'Homininae', 'Hominini', 'Homo sapiens']

LiftOver_crosstab = pd.crosstab(Dfam_RM_fil['repeat subfamily name'], LiftOver_df_fil['branch'], normalize='index')[branch_list] * 100

# add branch length for inference of root
LiftOver_df_fil['branch length'] = LiftOver_df_fil['repeat name'].apply(lambda x:branch_length_dict[x].branch_length if x in branch_length_dict.keys() else np.nan)
Dfam_RM_fil['branch length'] = Dfam_RM_fil['repeat name'].apply(lambda x:branch_length_dict[x].branch_length if x in branch_length_dict.keys() else np.nan)

branch_length_groupby = Dfam_RM_fil[['repeat subfamily name', 'branch length']].groupby(by='repeat subfamily name').mean()
Dfam_RM_Liftover_metadata = pd.concat([Dfam_RM_family_fil, LiftOver_crosstab, branch_length_groupby['branch length']], axis=1)


print(time.time()-s, 'Processing done')


### process and output data for each TE families

# screening result
for family in family_list:
    df = screening_result[(screening_result['repeat family name']==family)]
    df.to_csv('../data/screening/{}_screening.csv'.format(family), index=False)

# peak data
for i, family in enumerate(family_list):

    KZFP = Dfam_RM_overlap_KZFP_fil[Dfam_RM_overlap_KZFP_fil['repeat family name']==family]
    TRIM28 = Dfam_RM_overlap_TRIM28_fil[Dfam_RM_overlap_TRIM28_fil['repeat family name']==family]

    KZFP.to_csv('../data/overlap/{}_KZFP.csv'.format(family), index=False)
    TRIM28.to_csv('../data/overlap/{}_TRIM28.csv'.format(family), index=False)

# Liftover
for i, family in enumerate(family_list):

    subfamily_list = Dfam_RM_family_fil[Dfam_RM_family_fil['repeat family name']==family]
    Dfam = Dfam_RM_fil[(Dfam_RM_fil['repeat subfamily name'].isin(subfamily_list.index)) & (Dfam_RM_fil['repeat name'].isin(map_dict.keys()))]

    Lift = LiftOver_crosstab.loc[subfamily_list.index]
    summary = LiftOver_summary_fil.loc[Dfam['repeat name']]

    Lift.to_csv('../data/Liftover/{}_Liftover_rate.csv'.format(family))
    summary.to_csv('../data/Liftover/{}_Liftover.csv'.format(family))

# rerooted tree
for family in family_list:

    Lift = Dfam_RM_Liftover_metadata[Dfam_RM_Liftover_metadata['repeat family name']==family].sort_values(by=branch_list, ascending=False)
    branch, subfamily = Lift.iloc[0][['branch', 'repeat subfamily name']]

    # because the oldest subfamily of L1P_5end is L1PA17_5end.
    if family == 'L1P_5end':
        subfamily = 'L1PA17_5end'

    Dfam = Dfam_RM_fil[(Dfam_RM_fil['repeat subfamily name']==subfamily) & (LiftOver_df_fil['branch']==branch)].sort_values(by='branch length')
    Dfam = Dfam[Dfam['branch length'].isna()==False]
    outgroup = Dfam.iloc[-1]['repeat name']
    print(family, subfamily, len(Dfam), outgroup)

    tree = tree_dict[family]
    tree.root_with_outgroup(outgroup)    

    tree_dict[family] = tree
    Phylo.write(tree, '../data/phylogenetic tree/{}_reroot.contree'.format(family), format='newick')


# TE annotation with the position in tree
for family in family_list:

    Dfam = Dfam_RM_fil[Dfam_RM_fil['repeat family name']==family]
    tree = bt.loadNewick('../data/phylogenetic tree/{}_reroot.contree'.format(family))
    order = return_yaxis(tree).sort_values('y')

    Dfam = Dfam.loc[order['repeat adjusted name']]
    Dfam['branch length'] = Dfam['repeat name'].apply(lambda x:branch_length_dict[x].branch_length)

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

    for subfamily in pd.unique(df['repeat subfamily name']):

        subfamily_df = df[df['repeat subfamily name']==subfamily]
        subfamily_MSA = MSA.loc[subfamily_df['repeat name']]

        conserved_list = [pd.DataFrame(index=['A', 'C', 'G', 'T', '-'])]
        for col in subfamily_MSA.columns:

            count = subfamily_MSA[col].apply(lambda x:x.upper()).value_counts(normalize=True).rename(col)
            conserved_list.append(count)
        
        conserved_df = pd.concat(conserved_list, axis=1).fillna(0)
        conserved_df.to_csv('../data/MSA/conserved/{}_conserved.csv'.format(subfamily))


print(time.time()-s, 'All data outputed')
