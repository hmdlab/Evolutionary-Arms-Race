import pandas as pd
import numpy as np
from Bio import SeqIO
import time
from Bio import Phylo
import baltic as bt
from function import return_position_in_repeat_alignment, return_yaxis
import warnings
warnings.simplefilter('ignore')

family_list = ['MER11']
family = 'MER11'

s = time.time()
print('start processing MER11 data')

### load data
repeat_columns = ['repeat chr', 'repeat start', 'repeat end', 'repeat name', 'repeat score', 'repeat strand', 'repeat subfamily name', 'repeat family name', 
                  'repeat class', 'repeat classification', 'repeat region', 'repeat length', 'repeat substitution', 'repeat deletion', 'repeat insertion', 
                  'repeat consensus length', 'repeat alignment length', 'repeat no alignment length', 'repeat alignment ratio', 'repeat no alignment ratio',
                  'repeat substitution rate', 'repeat substitution rate filtering', 'repeat length filtering', 'repeat zerocoverage filtering', 'repeat filtering', 'used for phylogenetic tree']

TRIM28_columns = ['TRIM28 summits chr', 'TRIM28 summits start', 'TRIM28 summits end', 'TRIM28 peak name', 'TRIM28 score', 'TRIM28 strand', 'TRIM28 signal value',
                  'TRIM28 p value', 'TRIM28 q value', 'TRIM28 peak chr', 'TRIM28 peak start', 'TRIM28 peak end', 'TRIM28 peak length', 'TRIM28 accession',
                  'TRIM28 experiment', 'TRIM28 assay', 'TRIM28 run type', 'TRIM28 Cell type', 'TRIM28 Cell type class', 'TRIM28 control']

KZFP_columns = ['KZFP summits chr', 'KZFP summits start', 'KZFP summits end', 'KZFP peak name', 'KZFP score', 'KZFP strand', 'KZFP signal value',
                'KZFP p value', 'KZFP q value', 'KZFP peak chr', 'KZFP peak start', 'KZFP peak end', 'KZFP peak length', 'KZFP accession',
                'KZFP gene symbol', 'KZFP experiment', 'KZFP dataset origin', 'KZFP assay', 'KZFP run type', 'KZFP Cell type', 'KZFP Cell type class', 'KZFP control']

# TE metadata
Dfam_RM_family = pd.read_csv('../data/TE/TE_metadata_with_Liftover.csv')
Dfam_RM_family.index = Dfam_RM_family['repeat subfamily name']

# MER11 annotation
Dfam_RM_MER11 = pd.read_csv('../data/TE/MER11_annotation.csv')
Dfam_RM_MER11.index = Dfam_RM_MER11['repeat name']

# overlap with KZFP peaks
Dfam_RM_overlap_KZFP = pd.read_table('../data/overlap/TE_KZFP_overlap.bed', header=None)
Dfam_RM_overlap_KZFP.columns = repeat_columns + KZFP_columns + ['overlap length']

# overlap with TRIM28 peaks
Dfam_RM_overlap_TRIM28 = pd.read_table('../data/overlap/TE_TRIM28_overlap.bed', header=None)
Dfam_RM_overlap_TRIM28.columns = repeat_columns + TRIM28_columns + ['overlap length']

# Liftover
LiftOver_df = pd.read_csv('../data/Liftover/TE_Liftover.csv')
LiftOver_df.index = LiftOver_df['repeat name']

LiftOver_summary = pd.read_csv('../data/Liftover/TE_Liftover_summary.csv', index_col=0)

# screening for results
screening_result = pd.read_csv('../data/screening/screening_TE_evasion_candidates.csv')

# map files
map_dict = dict()
for i, family in enumerate(family_list):

    with open('../data/MSA/{}_nodup.fasta.map'.format(family)) as f:
        
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
MER11_tree = Phylo.read('../data/phylogenetic tree/MER11_nodup.contree', format='newick')

MER11_branch_length_dict = dict()
for node in MER11_tree.get_terminals():

    name = node.name
    branch_length = node.branch_length

    MER11_branch_length_dict[name] = node


print(time.time()-s, 'Loading done')



### process data

# extract MER11 copies
Dfam_RM_overlap_KZFP_fil = Dfam_RM_overlap_KZFP[Dfam_RM_overlap_KZFP['repeat name'].isin(Dfam_RM_MER11.index)]
Dfam_RM_overlap_TRIM28_fil = Dfam_RM_overlap_TRIM28[Dfam_RM_overlap_TRIM28['repeat name'].isin(Dfam_RM_MER11.index)]
LiftOver_df_fil = LiftOver_df.loc[Dfam_RM_MER11.index]
LiftOver_summary_fil = LiftOver_summary.loc[Dfam_RM_MER11.index]
Dfam_RM_MER11['branch'] = LiftOver_df_fil['branch']

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

LiftOver_crosstab = pd.crosstab(Dfam_RM_MER11['repeat adjusted subfamily name'], LiftOver_df_fil['branch'], normalize='index') * 100


# add branch length
Dfam_RM_MER11['branch length'] = Dfam_RM_MER11['repeat name'].apply(lambda x:MER11_branch_length_dict[x].branch_length if x in MER11_branch_length_dict.keys() else np.nan)
LiftOver_df_fil['branch length'] = LiftOver_df_fil['repeat name'].apply(lambda x:MER11_branch_length_dict[x].branch_length if x in MER11_branch_length_dict.keys() else np.nan)


print(time.time()-s, 'processing done')


### output for MER11

# screening
for family in family_list:
    df = screening_result[(screening_result['repeat family name']==family)]
    df.to_csv('../data/screening/{}_screening.csv'.format(family), index=False)

# peaks
Dfam_RM_overlap_KZFP_fil.to_csv('../data/overlap/{}_KZFP.csv'.format(family), index=False)
Dfam_RM_overlap_TRIM28_fil.to_csv('../data/overlap/{}_TRIM28.csv'.format(family), index=False)

# Liftover
for i, family in enumerate(family_list):

    Lift = LiftOver_crosstab
    summary = LiftOver_summary_fil.loc[Dfam_RM_MER11['repeat name']]

    Lift.to_csv('../data/Liftover/{}_Liftover_rate.csv'.format(family))
    summary.to_csv('../data/Liftover/{}_Liftover.csv'.format(family))

# tree
data = Dfam_RM_MER11[(Dfam_RM_MER11['repeat adjusted subfamily name']=='MER11_2') & (Dfam_RM_MER11['branch']=='Catarrhini')].sort_values(by='branch length', ascending=False)
outgroup = data.index[0]

MER11_tree.root_with_outgroup(outgroup_targets=outgroup)
Phylo.write(trees=MER11_tree, file='../data/phylogenetic tree/MER11_reroot.contree', format='newick')


# TE annotation
for family in family_list:

    Dfam = Dfam_RM_MER11
    tree = bt.loadNewick('../data/phylogenetic tree/{}_reroot.contree'.format(family))
    order = return_yaxis(tree).sort_values('y')
    order_fil = order[order['repeat adjusted name'].isin(Dfam.index)]

    # 調整する
    Dfam = Dfam.loc[order_fil['repeat adjusted name']]
    Dfam['branch length'] = Dfam['repeat name'].apply(lambda x:MER11_branch_length_dict[x].branch_length)

    # 出力する
    Dfam.to_csv('../data/TE/{}.annotation.csv'.format(family), index=False)
    
# MSA
MSA_dict = dict()

for family in family_list:

    seq_list = list()
    name_list = list()

    with open('../data/MSA/{}_nodup_alignment_to_consensus_merge.fasta'.format(family)) as f:

        for record in SeqIO.parse(f, 'fasta'):

            seq_list.append(list(str(record.seq)))
            name_list.append(record.id)

    MSA = pd.DataFrame(seq_list, index=name_list)
    MSA_dict[family] = MSA

    MSA.to_csv('../data/MSA/{}_MSA.csv'.format(family))


# conserved level in each position 
for i, family in enumerate(family_list):

    MSA = MSA_dict[family]
    df = Dfam_RM_MER11
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
