import pandas as pd
import numpy as np
from scipy.stats import binomtest
from statsmodels.stats.multitest import fdrcorrection
from Bio import SeqIO
import time

s = time.time()
print('start enrichment analysis of proviruses')

### load data

# provirus metadata
Dfam_RM_family = pd.read_csv('../data/TE/provirus_metadata_with_Liftover.csv')
Dfam_RM_family.index = Dfam_RM_family['repeat adjusted subfamily name']
Dfam_RM_family_dict = Dfam_RM_family.to_dict()


# provirus annotation
Dfam_RM = pd.read_csv('../data/TE/provirus_annotation.csv')
Dfam_RM.index = Dfam_RM['repeat name']
Dfam_RM

KZFP_columns = ['KZFP summits chr', 'KZFP summits start', 'KZFP summits end', 'KZFP peak name', 'KZFP score', 'KZFP strand', 'KZFP signal value',
                'KZFP p value', 'KZFP q value', 'KZFP peak chr', 'KZFP peak start', 'KZFP peak end', 'KZFP peak length', 'KZFP accession',
                'KZFP gene symbol', 'KZFP experiment', 'KZFP dataset origin', 'KZFP assay', 'KZFP run type', 'KZFP Cell type', 'KZFP Cell type class', 'KZFP control']

# KZFP metadata
KZFP_dataset_df = pd.read_csv('../data/ChIP-seq/KZFP_metadata.csv', index_col=0)
KZFP_dataset_dict = KZFP_dataset_df.to_dict()

# KZFP ChIP-seq metadata
KZFP_metadata = pd.read_csv('../data/ChIP-seq/KZFP_ChIP-seq_metadata.csv')

# KZFP peaks
KZFP_peak_df = pd.read_table('../data/ChIP-seq/KZFP_peaks.bed', header=None)
KZFP_peak_df.columns = KZFP_columns

# overlap between provirus and KZFP
Dfam_RM_overlap_KZFP = pd.read_table('../data/overlap/provirus_KZFP_overlap.bed', header=None)
Dfam_RM_overlap_KZFP.columns = Dfam_RM.columns.tolist() + KZFP_columns + ['overlap length']


# hg38 genome
hg38_dict = dict()
length = 0

with open('../data/UCSC/hg38.fa') as f:

    for record in SeqIO.parse(f, 'fasta'):

        hg38_dict[record.id] = record


# unmapped regions
zero_coverage_df = pd.read_csv('../data/ChIP-seq/405samples_zerocoverage.csv')
zero_coverage_df.head()


print('Loading done')

### processed data

# filter non-full-length provirus
Dfam_RM_fil = Dfam_RM[(Dfam_RM['5-LTR name'].isna()==False) & (Dfam_RM['3-LTR name'].isna()==False) & (Dfam_RM['repeat adjusted subfamily name'].isna()==False)]
Dfam_RM_overlap_KZFP_fil = Dfam_RM_overlap_KZFP[Dfam_RM_overlap_KZFP['repeat name'].isin(Dfam_RM_fil['repeat name'])].copy()

# write subfamily name
Dfam_RM_fil_dict = Dfam_RM_fil.to_dict()
Dfam_RM_overlap_KZFP_fil['repeat adjusted subfamily name'] = Dfam_RM_overlap_KZFP_fil['repeat name'].apply(lambda x:Dfam_RM_fil_dict['repeat adjusted subfamily name'][x])


# return peak summit position in repeat
def return_position_in_repeat(summit, r_s, r_e, strand):

    if strand == '+':

        return summit - r_s
    
    else:

        return r_e - summit


# return peak summit position in repeat region (5-LTR, internal region, 3-LTR)
def return_region_in_repeat(summit, LTR5_s, LTR5_e, LTR3_s, LTR3_e, strand):
            
    if LTR5_s <= summit and summit <= LTR5_e :
                
        return '5-LTR'
            
    elif LTR3_s <= summit and summit <= LTR3_e:

        return '3-LTR'

    else:

        if (LTR5_e <= summit and summit <= LTR3_s) or (LTR3_e <= summit and summit <= LTR5_s):

            return 'Int'

        else:

            return None


position_in_repeat_list = list()
for summit, r_s, r_e, strand in Dfam_RM_overlap_KZFP_fil[['KZFP summits start', 'repeat start', 'repeat end', 'repeat strand']].values:

    position = return_position_in_repeat(summit, r_s, r_e, strand)
    position_in_repeat_list.append(position)

Dfam_RM_overlap_KZFP_fil['summit start in repeat'] = position_in_repeat_list

position_in_repeat_list = list()
for summit, LTR5_s, LTR5_e, LTR3_s, LTR3_e, strand in Dfam_RM_overlap_KZFP_fil[['KZFP summits start', '5-LTR start', '5-LTR end', '3-LTR start', '3-LTR end', 'repeat strand']].values:

    position = return_region_in_repeat(summit, LTR5_s, LTR5_e, LTR3_s, LTR3_e, strand)
    position_in_repeat_list.append(position)

Dfam_RM_overlap_KZFP_fil['summit region in repeat'] = position_in_repeat_list 


# fileter chromosomes
chrom_list = ['chr{}'.format(i) for i in range(1, 23)] + ['chrX']

Dfam_RM_fil_chrom = Dfam_RM_fil[Dfam_RM_fil['repeat chr'].isin(chrom_list)]
KZFP_peak_df_chrom = KZFP_peak_df[KZFP_peak_df['KZFP summits chr'].isin(chrom_list)]


# calculate effective genome length
chrom_length = 0
for chrom, record in hg38_dict.items():

    if chrom in chrom_list:

        chrom_length += len(record.seq)

zero_coverage_df_fil_chrom = zero_coverage_df[(zero_coverage_df['N=405']) & (zero_coverage_df['chr'].isin(chrom_list))]
remove_length = (zero_coverage_df_fil_chrom['end'] - zero_coverage_df_fil_chrom['start']).sum()
effective_length = chrom_length - remove_length


# obtain provirus subfamily annotation and length in genome
Dfam_cluster_chrom_dict = dict()
Dfam_cluster_chrom_length_dict = dict()

for i, subfamily in enumerate(Dfam_RM_family['repeat adjusted subfamily name']):

    # clusterの情報
    df = Dfam_RM_fil_chrom[Dfam_RM_fil_chrom['repeat adjusted subfamily name']==subfamily]
    Dfam_cluster_chrom_dict[subfamily] = df

    # 長さの情報
    LTR5_length = (df['5-LTR end'] - df['5-LTR start']).sum()
    LTR3_length = (df['3-LTR end'] - df['3-LTR start']).sum()
    Int_length = (df['repeat end'] - df['repeat start']).sum() - LTR5_length - LTR3_length
    
    Dfam_cluster_chrom_length_dict[subfamily] = [LTR5_length, LTR3_length, Int_length]

print(time.time()-s, 'Processing done')


# count peak
KZFP_peak_num = KZFP_peak_df_chrom['KZFP gene symbol'].value_counts()

# count repeat
repeat_count = Dfam_RM_fil_chrom['repeat adjusted subfamily name'].value_counts()

# count overlap
Dfam_RM_overlap_KZFP_fil_chrom = Dfam_RM_overlap_KZFP_fil[Dfam_RM_overlap_KZFP_fil['repeat name'].isin(Dfam_RM_fil_chrom['repeat name'])]

crosstab_dict = dict()
for region in ['5-LTR', '3-LTR', 'Int']:
    data = Dfam_RM_overlap_KZFP_fil_chrom[Dfam_RM_overlap_KZFP_fil_chrom['summit region in repeat']==region]
    crosstab = pd.crosstab(data['repeat adjusted subfamily name'], data['KZFP gene symbol'])
    crosstab_dict[region] = crosstab

repeat_overlap_count_dict = dict()
for region in ['5-LTR', '3-LTR', 'Int']:
    data = Dfam_RM_overlap_KZFP_fil_chrom[Dfam_RM_overlap_KZFP_fil_chrom['summit region in repeat']==region]
    Dfam_RM_overlap_KZFP_fil_nodup = data[(data[['repeat name', 'KZFP gene symbol']].duplicated()==False)]
    repeat_overlap_count = Dfam_RM_overlap_KZFP_fil_nodup[['repeat adjusted subfamily name', 'KZFP gene symbol']].value_counts()
    repeat_overlap_count_dict[region] = repeat_overlap_count


print(time.time()-s, 'Counting done')



### binomial test
ChIP_count = KZFP_metadata['KZFP gene symbol'].value_counts()
binom_result_list = list()
for i, cluster in enumerate(pd.unique(Dfam_RM_family['repeat adjusted subfamily name'])):

    for i, region in enumerate(['5-LTR', '3-LTR', 'Int']):

        for KZFP in crosstab.columns:

            try:

                df = crosstab_dict[region].loc[cluster]
                overlap_for_proportion = repeat_overlap_count_dict[region][(cluster, KZFP)]
            
            except:

                binom_result_list.append([KZFP, cluster, region, 0, 1, 0])

                continue

            # obtain count
            overlap = df[KZFP]
            peak_num = KZFP_peak_num[KZFP]
            length = Dfam_cluster_chrom_length_dict[cluster][i]
            excepted = length / effective_length

            # bionomial test
            overlap = min(overlap, peak_num)
            p = binomtest(overlap, n=peak_num, p=excepted, alternative='greater').pvalue
            ratio = (overlap/peak_num) / excepted

            binom_result_list.append([KZFP, cluster, region, ratio, p, overlap, overlap_for_proportion])
    
# convert to DataFrame
binom_result_df = pd.DataFrame(binom_result_list, columns=['KZFP gene symbol', 'repeat adjusted subfamily name', 'target region', 'ratio', 'p value', 'overlap peak count to all copies', 'copy count overlaped with peaks'])

# adjusted p-value
binom_result_df['p value adjusted'] = np.where(binom_result_df['p value']==0, binom_result_df[binom_result_df['p value']!=0]['p value'].min(), binom_result_df['p value'])
binom_result_df['log10 p value'] = binom_result_df['p value adjusted'].apply(np.log10).apply(abs)

print(time.time()-s, 'Binomial test done')



### define primary and secondary targets

# definition
def return_primary(p, v):

    if p < 0.05:

        if v >= 0.9:

            return 'Primary'
        
        elif v < 0.9:

            return 'Secondary'
        
    else:

        return False
    
binom_result_df = binom_result_df.sort_values(by=['KZFP gene symbol', 'log10 p value', 'ratio'], ascending=False)
binom_result_target_df = list()

for KZFP in pd.unique(binom_result_df['KZFP gene symbol']):

    df = binom_result_df[(binom_result_df['KZFP gene symbol']==KZFP)].copy()

    # calculate q-value
    df['q value'] = fdrcorrection(df['p value adjusted'])[1]
    df['log10 q value'] = df['q value'].apply(np.log10).apply(abs)

    # normalize q-value
    primary = df['log10 q value'] / df['log10 q value'].iloc[0]
        
    target = list()
    for qvalue, value in zip(df['q value'], primary):

        target.append(return_primary(qvalue, value))
        
    df['normalized score'] = primary
    df['rank'] = target

    binom_result_target_df.append(df)


binom_result_target_df = pd.concat(binom_result_target_df, axis=0)
binom_result_target_df.head()


### Write additional information

binom_result_target_df['copy number'] = binom_result_target_df['repeat adjusted subfamily name'].apply(lambda x:repeat_count[x])
binom_result_target_df['proportion of copy overlaped with peaks'] = binom_result_target_df['copy count overlaped with peaks'] / binom_result_target_df['copy number']
binom_result_target_df['binding rate'] = binom_result_target_df['overlap peak count to all copies'] / binom_result_target_df['copy number']
binom_result_target_df['binding rate'] = [binding / ChIP_count[KZFP] for binding, KZFP in binom_result_target_df[['binding rate', 'KZFP gene symbol']].values]
binom_result_target_df['repeat family name'] = binom_result_target_df['repeat adjusted subfamily name'].apply(lambda x:Dfam_RM_family_dict['repeat family name'][x])
binom_result_target_df['repeat class'] = binom_result_target_df['repeat adjusted subfamily name'].apply(lambda x:Dfam_RM_family_dict['repeat class'][x])
binom_result_target_df['repeat classification'] = binom_result_target_df['repeat adjusted subfamily name'].apply(lambda x:Dfam_RM_family_dict['repeat classification'][x])


# TEファミリーとKZFPの出現年代を記録していく

branch_dict = {'Vertebrata':530, 'Dipnotetrapodomorpha':413, 'Tetrapoda':352, 'Amniota':312, 'Mammalia':179.2, 'Theria':159, 'Eutheria':105, 'Boreoeutheria':96, 'Euarchontoglires': 90, 'Primatomorpha':76, 
               'Primates':74, 'Haplorrhini':63, 'Simiiformes':43.2, 'Catarrhini':29.4, 'Hominoidea':20.2, 'Hominidae':15.8, 'Homininae':9.1, 'Hominini':6.7, 'Homo sapiens':0, np.nan:None}

Age_dict = {352.0: 'Tetrapoda', 312.1: 'Amniota', 312.0: 'Amniota', 163.7:'Mammalia', 159.0: 'Theria', 105.0: 'Eutheria', 96.0: 'Boreoeutheria', 97.5: 'Boreoeutheria', 90.0: 'Euarchontoglires', 90.9: 'Euarchontoglires', 76: 'Primatomorpha', 75.9: 'Primatomorpha', 
            74.0: 'Primates', 63: 'Haplorrhini', 43.2:'Simiiformes', 43.1:'Simiiformes', 29.4: 'Catarrhini', 29.1: 'Catarrhini', 20.2:'Hominoidea', 15.8: 'Hominidae', 9.1: 'Homininae', 6.7: 'Hominini', 0: 'Homo sapiens', 'nan': np.nan}

def return_Age(x, dict):

    if x in dict.keys():

        return dict[x]
    
    else:

        return None

binom_result_target_df['emergence era of TE subfamily'] = binom_result_target_df['repeat adjusted subfamily name'].apply(lambda x:Dfam_RM_family_dict['branch'][x])
binom_result_target_df['evolutionary age of TE subfamily'] = binom_result_target_df['emergence era of TE subfamily'].apply(lambda x:return_Age(x, branch_dict))
binom_result_target_df['evolutionary age of KZFP for analysis'] = binom_result_target_df['KZFP gene symbol'].apply(lambda x:return_Age(x, KZFP_dataset_dict['Age adjusted']))
binom_result_target_df['emergence era of KZFP for analysis'] = binom_result_target_df['evolutionary age of KZFP for analysis'].apply(lambda x:return_Age(x, Age_dict))
binom_result_target_df['evolutionary age of KZFP in Imbeault et al.'] = binom_result_target_df['KZFP gene symbol'].apply(lambda x:return_Age(x, KZFP_dataset_dict['Age']))
binom_result_target_df['evolutionary age of KZFP in Tribolet-Hardy et al.'] = binom_result_target_df['KZFP gene symbol'].apply(lambda x:return_Age(x, KZFP_dataset_dict['Age (Genome research)']))

# primateかどうかを記録する
binom_result_target_df['emegence era of TE subfamily is in primate'] = binom_result_target_df['evolutionary age of TE subfamily'] <= 74
binom_result_target_df['emegence era of KZFP is in primate'] = binom_result_target_df['evolutionary age of KZFP for analysis'] <= 74

# KZFPのデータセットにもbranchを挿入する
KZFP_dataset_df['branch'] = KZFP_dataset_df['Age'].apply(lambda x:return_Age(x, Age_dict))


print(time.time()-s, 'Writing additional information done')



### output
rename_columns = ['KZFP gene symbol', 'repeat adjusted subfamily name', 'target region', 'ratio', 'p value', 'q value', 'log10 q value', 'normalized score', 'rank',
                  'overlap peak count to all copies', 'copy count overlaped with peaks', 'copy number', 'proportion of copy overlaped with peaks', 'binding rate', 
                  'repeat family name', 'repeat class', 'repeat classification', 'emergence era of TE subfamily', 'evolutionary age of TE subfamily', 
                  'evolutionary age of KZFP for analysis', 'emergence era of KZFP for analysis', 'evolutionary age of KZFP in Imbeault et al.', 'evolutionary age of KZFP in Tribolet-Hardy et al.']


condition1 = (binom_result_target_df['log10 q value']>10)
condition2 = (binom_result_target_df['ratio']>2)
condition3 = (binom_result_target_df['overlap peak count to all copies']>=5)
condition4 = (binom_result_target_df['proportion of copy overlaped with peaks']>=0.1)

binom_result_target_df['criteria'] = condition1 & condition2 & condition3 & condition4

binom_result_target_df[rename_columns].to_csv('../data/targets/provirus_targets_raw.csv', index=False)
binom_result_target_df[binom_result_target_df['criteria']][rename_columns].to_csv('../data/targets/provirus_targets_for_analysis.csv', index=False)


print(time.time()-s, 'finish enrichment analysis of proviruses')
