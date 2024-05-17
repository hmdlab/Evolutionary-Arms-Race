import pandas as pd
import numpy as np
import time
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import warnings
warnings.simplefilter('ignore')

s = time.time()
print('start screening for evasions')


### load data

# TE metadata
Dfam_RM_family = pd.read_csv('../data/TE/TE_metadata_with_Liftover.csv')
Dfam_RM_family.index = Dfam_RM_family['repeat subfamily name']
Dfam_RM_family_dict = Dfam_RM_family.to_dict()

# TE annotation
Dfam_RM = pd.read_csv('../data/TE/TE_annotation.csv')
Dfam_RM.index = Dfam_RM['repeat name']

KZFP_columns = ['KZFP summits chr', 'KZFP summits start', 'KZFP summits end', 'KZFP peak name', 'KZFP score', 'KZFP strand', 'KZFP signal value',
                'KZFP p value', 'KZFP q value', 'KZFP peak chr', 'KZFP peak start', 'KZFP peak end', 'KZFP peak length', 'KZFP accession',
                'KZFP gene symbol', 'KZFP experiment', 'KZFP dataset origin', 'KZFP assay', 'KZFP run type', 'KZFP Cell type', 'KZFP Cell type class', 'KZFP control']


# KZFP metadata
KZFP_dataset_df = pd.read_csv('../data/ChIP-seq/KZFP_metadata.csv', index_col=0)
KZFP_dataset_dict = KZFP_dataset_df.to_dict()

# KZFP ChIP-seq metadata
KZFP_metadata = pd.read_csv('../data/ChIP-seq/KZFP_ChIP-seq_metadata.csv')

# overlap between TE and KZFP
Dfam_RM_overlap_KZFP = pd.read_table('../data/overlap/TE_KZFP_overlap.bed', header=None)
Dfam_RM_overlap_KZFP.columns = Dfam_RM.columns.tolist() + KZFP_columns + ['overlap length']

# KZFP targets
KZFP_target = pd.read_csv('../data/targets/TE_targets_for_analysis.csv')
print('KZFP targets: N={}'.format(len(KZFP_target)))


print(time.time()-s, 'Loading done')


### process data

# remove TE subfamilies that have low copies (<10) and TE family without subfamily classification
Dfam_RM_family_fil = Dfam_RM_family[Dfam_RM_family['repeat copy count after filtering']>=10]
family_num = Dfam_RM_family_fil['repeat family name'].value_counts()
Dfam_RM_family_fil = Dfam_RM_family_fil[Dfam_RM_family_fil['repeat family name'].isin(family_num[family_num>=2].index)]

# extract TE 
TE_class_list = ['DNA', 'ERV/LTR', 'ERV/Int', 'LINE', 'SINE', 'Retroposon']
Dfam_RM_family_fil = Dfam_RM_family_fil[Dfam_RM_family_fil['repeat class'].isin(TE_class_list)]

KZFP_target_fil = KZFP_target[KZFP_target['repeat subfamily name'].isin(Dfam_RM_family_fil['repeat subfamily name'])]

# extract full-length TE
Dfam_RM_fil = Dfam_RM.copy()
Dfam_RM_fil = Dfam_RM_fil[Dfam_RM['repeat filtering']]
Dfam_RM_overlap_KZFP_fil = Dfam_RM_overlap_KZFP[Dfam_RM_overlap_KZFP['repeat name'].isin(set(Dfam_RM_fil['repeat name']))].copy()

# add repeat family name
Dfam_RM_overlap_KZFP_fil['repeat family name'] = Dfam_RM_overlap_KZFP_fil['repeat subfamily name'].apply(lambda x:Dfam_RM_family_dict['repeat family name'][x])


# create family dictionary
family_annotation_dict = dict()
for cluster, family, repeat_class, region in Dfam_RM_family[['repeat subfamily name', 'repeat family name', 'repeat class', 'repeat region']].values:

    family_annotation_dict[family] = [repeat_class, region]

print(time.time()-s, 'processing done')



### count data

# obtain overlap data for each TE-KZFP association
Dfam_RM_fil_family_dict = dict()
Dfam_RM_family_overlap_KZFP_dict = dict()
for i, (family, KZFP) in enumerate(KZFP_target_fil[['repeat family name', 'KZFP gene symbol']].value_counts().index):

    # record TE family data
    df = Dfam_RM_fil[Dfam_RM_fil['repeat family name']==family]
    Dfam_RM_fil_family_dict[family] = df

    # record the overlap of TE family and KRAB-ZFP
    df = Dfam_RM_overlap_KZFP_fil[(Dfam_RM_overlap_KZFP_fil['repeat family name']==family) & (Dfam_RM_overlap_KZFP_fil['KZFP gene symbol']==KZFP)]
    Dfam_RM_family_overlap_KZFP_dict[(family, KZFP)] = df


# calculate binding rate
binding_rate_dict = dict()
ChIP_count = KZFP_metadata['KZFP gene symbol'].value_counts()
for i, (family, KZFP) in enumerate(KZFP_target_fil[['repeat family name', 'KZFP gene symbol']].value_counts().index):

    # obtain data
    Dfam = Dfam_RM_fil_family_dict[family]
    KZFP_overlap = Dfam_RM_family_overlap_KZFP_dict[(family, KZFP)]

    # calculate binding rate
    crosstab = pd.crosstab(KZFP_overlap['repeat name'], KZFP_overlap['KZFP accession']) >= 1
    count = crosstab.sum(axis=1).rename('KZFP binding') / ChIP_count[KZFP]
    
    data = pd.concat([Dfam, count], axis=1)
    data['KZFP binding'] = data['KZFP binding'].fillna(0)

    binding_rate_dict[(family, KZFP)] = data


# record binding rate of TE subfamilies
binding_rate_subfamily_dict = dict()
for i, (family, KZFP) in enumerate(KZFP_target_fil[['repeat family name', 'KZFP gene symbol']].value_counts().index):

    data = binding_rate_dict[(family, KZFP)]
    groupby = data[['repeat subfamily name', 'KZFP binding']].groupby(by=['repeat subfamily name']).mean()

    for subfamily, binding_rate in zip(groupby.index, groupby['KZFP binding']):

        binding_rate_subfamily_dict[(subfamily, KZFP)] = binding_rate


print(time.time()-s, 'counting done')



### screening
tukeyhsd_list = list()
s = time.time()

for i, (family, KZFP) in enumerate(KZFP_target_fil[['repeat family name', 'KZFP gene symbol']].value_counts().index):

    data = binding_rate_dict[(family, KZFP)]

    # turkeyの検定
    binding_list = data['KZFP binding']
    sample_list = data['repeat subfamily name']
    
    tukeyhsd = pairwise_tukeyhsd(binding_list, sample_list)
    result_table = tukeyhsd._results_table
    df = pd.DataFrame(result_table.data[1:], columns=['group1', 'group2', 'meandiff', 'p-adj', 'lower', 'upper', 'reject'])

    df['repeat family name'] = family
    df['KZFP gene symbol'] = KZFP

    # 記録する
    tukeyhsd_list.append(df)

print(time.time()-s, 'screening done')



### write additional information
branch_dict = {'Vertebrata':530, 'Dipnotetrapodomorpha':413, 'Tetrapoda':352, 'Amniota':312, 'Mammalia':179.2, 'Theria':159, 'Eutheria':105, 'Boreoeutheria':96, 'Euarchontoglires': 90, 'Primatomorpha':76, 
               'Primates':74, 'Haplorrhini':63, 'Simiiformes':43.2, 'Catarrhini':29.4, 'Hominoidea':20.2, 'Hominidae':15.8, 'Homininae':9.1, 'Hominini':6.7, 'Homo sapiens':0, np.nan:None}

tukeyhsd_df = pd.concat(tukeyhsd_list, axis=0)
tukeyhsd_df.columns = ['subfamily1', 'subfamily2', 'difference of binding rate', 'adjusted p value', 'lower', 'upper', 'reject', 'repeat family name', 'KZFP gene symbol']

# group1: high binding rate, group2: low binding rate
group_list = [[], []]
for group1, group2, meandiff in tukeyhsd_df[['subfamily1', 'subfamily2', 'difference of binding rate']].values:

    if meandiff > 0:

        group1, group2 = group2, group1

    group_list[0].append(group1)
    group_list[1].append(group2)
    
tukeyhsd_df['subfamily1'] = group_list[0]
tukeyhsd_df['subfamily2'] = group_list[1]
tukeyhsd_df['difference of binding rate'] = tukeyhsd_df['difference of binding rate'].apply(abs)

# サブファミリーの年齢を記載する
tukeyhsd_df['subfamily1 branch'] = tukeyhsd_df['subfamily1'].apply(lambda x:Dfam_RM_family_dict['branch'][x])
tukeyhsd_df['subfamily2 branch'] = tukeyhsd_df['subfamily2'].apply(lambda x:Dfam_RM_family_dict['branch'][x])

tukeyhsd_df['subfamily1 age'] = tukeyhsd_df['subfamily1 branch'].apply(lambda x:branch_dict[x])
tukeyhsd_df['subfamily2 age'] = tukeyhsd_df['subfamily2 branch'].apply(lambda x:branch_dict[x])
tukeyhsd_df['KZFP age'] = tukeyhsd_df['KZFP gene symbol'].apply(lambda x:KZFP_dataset_dict['Age adjusted'][x])

# classを出す
tukeyhsd_df['repeat class'] = tukeyhsd_df['repeat family name'].apply(lambda x:family_annotation_dict[x][0])
tukeyhsd_df['repeat region'] = tukeyhsd_df['repeat family name'].apply(lambda x:family_annotation_dict[x][1])


# decrease in younger subfamily?
def return_diff_age(meandiff, age1, age2):

    if age1 == age2:

        return 'same'
    
    elif age1 > age2:

        return 'decrease'
    
    else:

        return 'increase'
        
value_list = list()
for meandiff, age1, age2 in tukeyhsd_df[['difference of binding rate', 'subfamily1 age', 'subfamily2 age']].values:
    
    value = return_diff_age(meandiff, age1, age2)
    value_list.append(value)

tukeyhsd_df['decrease in younger subfamily?'] = value_list


# サブファミリーごとの結合率を記録する
value_list = [[], []]
for group1, group2, KZFP in tukeyhsd_df[['subfamily1', 'subfamily2', 'KZFP gene symbol']].values:

    value1 = binding_rate_subfamily_dict[(group1, KZFP)]
    value2 = binding_rate_subfamily_dict[(group2, KZFP)]

    value_list[0].append(value1)
    value_list[1].append(value2)

tukeyhsd_df['subfamily1 binding rate'] = value_list[0]
tukeyhsd_df['subfamily2 binding rate'] = value_list[1]


print(time.time()-s, 'Writing additional information done')



### output
rename_columns = ['repeat family name', 'KZFP gene symbol', 'subfamily1', 'subfamily2', 'difference of binding rate', 'decrease in younger subfamily?', 'adjusted p value', 'subfamily1 binding rate', 'subfamily2 binding rate',
                  'repeat class', 'repeat region', 'subfamily1 branch', 'subfamily2 branch', 'subfamily1 age', 'subfamily2 age', 'KZFP age']

condition1 = tukeyhsd_df['KZFP age']>=tukeyhsd_df['subfamily2 age']
condition2 = tukeyhsd_df['adjusted p value']<0.05
condition3 = tukeyhsd_df['difference of binding rate']>=0.1
condition4 = tukeyhsd_df['decrease in younger subfamily?']=='decrease'

decrease_df = tukeyhsd_df[condition1 & condition2 & condition3 & condition4].sort_values(by='difference of binding rate', ascending=False)
decrease_df = decrease_df[decrease_df[['repeat family name', 'KZFP gene symbol']].duplicated()==False]

tukeyhsd_df[rename_columns].to_csv('../data/screening/screening_TE_evasion_raw_data.csv', index=False)
decrease_df[rename_columns].to_csv('../data/screening/screening_TE_evasion_candidates.csv', index=False)

print('candidate for evasion event: N={}'.format(len(decrease_df)))
print(time.time()-s, 'finish screening for evasions')
