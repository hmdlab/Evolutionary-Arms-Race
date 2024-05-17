import pandas as pd
import numpy as np
from scipy.stats import binomtest
from statsmodels.stats.multitest import fdrcorrection
from Bio import SeqIO
import warnings
import time
warnings.simplefilter('ignore')

s = time.time()
print('start enrichment analysis of TEs')

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

# KZFP peaks
KZFP_peak_df = pd.read_table('../data/ChIP-seq/KZFP_peaks.bed', header=None)
KZFP_peak_df.columns = KZFP_columns

# overlap between TE and KZFP peaks
Dfam_RM_overlap_KZFP = pd.read_table('../data/overlap/TE_KZFP_overlap.bed', header=None)
Dfam_RM_overlap_KZFP.columns = Dfam_RM.columns.tolist() + KZFP_columns + ['overlap length']


# hg38 genome
hg38_dict = dict()
length = 0

with open('../data/UCSC/hg38.fa') as f:

    for record in SeqIO.parse(f, 'fasta'):

        hg38_dict[record.id] = record


# unmapped region
zero_coverage_df = pd.read_csv('../data/ChIP-seq/405samples_zerocoverage.csv')
zero_coverage_df.head()

print(time.time()-s, 'Loading done')



### processed data

# filter chromosomes
chrom_list = ['chr{}'.format(i) for i in range(1, 23)] + ['chrX']

Dfam_RM_chrom = Dfam_RM[Dfam_RM['repeat chr'].isin(chrom_list)]
KZFP_peak_df_chrom = KZFP_peak_df[KZFP_peak_df['KZFP summits chr'].isin(chrom_list)]
Dfam_RM_overlap_KZFP_chrom = Dfam_RM_overlap_KZFP[Dfam_RM_overlap_KZFP['repeat chr'].isin(chrom_list)]


# estimate effective genome length
chrom_length = 0
for chrom, record in hg38_dict.items():

    if chrom in chrom_list:

        chrom_length += len(record.seq)


zero_coverage_df_fil_chrom = zero_coverage_df[(zero_coverage_df['N=405']) & (zero_coverage_df['chr'].isin(chrom_list))]
remove_length = (zero_coverage_df_fil_chrom['end'] - zero_coverage_df_fil_chrom['start']).sum()

effective_length = chrom_length - remove_length

print(time.time()-s, 'Processing done')


### count overlap

# record TE subfamily annotation and length in genome
Dfam_cluster_chrom_dict = dict()
Dfam_cluster_chrom_length_dict = dict()

for i, subfamily in enumerate(Dfam_RM_family['repeat subfamily name']):

    df = Dfam_RM_chrom[Dfam_RM_chrom['repeat subfamily name']==subfamily]
    Dfam_cluster_chrom_dict[subfamily] = df.copy()

    length = (df['repeat end'] - df['repeat start']).sum()
    Dfam_cluster_chrom_length_dict[subfamily] = length


# count KZFP peaks
KZFP_peak_num = KZFP_peak_df_chrom['KZFP gene symbol'].value_counts()

# count repeat copy
repeat_count = Dfam_RM_chrom['repeat subfamily name'].value_counts()

# count overlap
Dfam_RM_overlap_KZFP_chrom_TE = Dfam_RM_overlap_KZFP_chrom[Dfam_RM_overlap_KZFP_chrom['repeat name'].isin(Dfam_RM_chrom['repeat name'])]
crosstab = pd.crosstab(Dfam_RM_overlap_KZFP_chrom_TE['repeat subfamily name'], Dfam_RM_overlap_KZFP_chrom_TE['KZFP gene symbol'])

# count repeat overlap with KRAB-ZFPs
Dfam_RM_overlap_KZFP_fil_nodup = Dfam_RM_overlap_KZFP_chrom_TE[(Dfam_RM_overlap_KZFP_chrom_TE[['repeat name', 'KZFP gene symbol']].duplicated()==False)]
repeat_overlap_count = Dfam_RM_overlap_KZFP_fil_nodup[['repeat subfamily name', 'KZFP gene symbol']].value_counts()

print(time.time()-s, 'Counting done')



### binomial test
binom_result_list = list()
for i, cluster in enumerate(pd.unique(Dfam_RM_family['repeat subfamily name'])):

    for KZFP in crosstab.columns:

        try:

            df = crosstab.loc[cluster]
            overlap_for_proportion = repeat_overlap_count[(cluster, KZFP)]
        
        except:
            binom_result_list.append([KZFP, cluster, 0, 1, 0])
            continue

        # obtain count
        overlap = df[KZFP]
        peak_num = KZFP_peak_num[KZFP]
        length = Dfam_cluster_chrom_length_dict[cluster]
        excepted = length / effective_length

        # binomial test
        overlap = min(overlap, peak_num)
        p = binomtest(overlap, n=peak_num, p=excepted, alternative='greater').pvalue
        ratio = (overlap/peak_num) / excepted

        binom_result_list.append([KZFP, cluster, ratio, p, overlap, overlap_for_proportion])

# convert to DataFrame
binom_result_df = pd.DataFrame(binom_result_list, columns=['KZFP gene symbol', 'repeat subfamily name', 'ratio', 'p value', 'overlap peak count to all copies', 'copy count overlaped with peaks'])

# adjust p value
binom_result_df['p value adjusted'] = np.where(binom_result_df['p value']==0, binom_result_df[binom_result_df['p value']!=0]['p value'].min(), binom_result_df['p value'])
binom_result_df['log10 p value'] = binom_result_df['p value adjusted'].apply(np.log10).apply(abs)

print(time.time()-s, 'Binomial test done')


### calculate overlap to full-length TE

# count repeat
Dfam_RM_chrom_intact = Dfam_RM_chrom[Dfam_RM_chrom['repeat filtering']]
repeat_count_intact = Dfam_RM_chrom_intact['repeat subfamily name'].value_counts()

# count overlap
Dfam_RM_overlap_KZFP_chrom_intact = Dfam_RM_overlap_KZFP_chrom[Dfam_RM_overlap_KZFP_chrom['repeat filtering']]
crosstab_intact = pd.crosstab(Dfam_RM_overlap_KZFP_chrom_intact['repeat subfamily name'], Dfam_RM_overlap_KZFP_chrom_intact['KZFP gene symbol'])

Dfam_RM_overlap_KZFP_nodup_intact = Dfam_RM_overlap_KZFP_chrom_intact[(Dfam_RM_overlap_KZFP_chrom_intact[['repeat name', 'KZFP gene symbol']].duplicated()==False)]
repeat_overlap_count_intact = Dfam_RM_overlap_KZFP_nodup_intact[['repeat subfamily name', 'KZFP gene symbol']].value_counts()


intact_repeat_overlap_list = list()
for i, cluster in enumerate(pd.unique(Dfam_RM_family['repeat subfamily name'])):

    for KZFP in crosstab.columns:

        # overlapが観測されなかったTEファミリーは飛ばす
        try:

            df = crosstab_intact.loc[cluster]
            overlap_for_proportion = repeat_overlap_count_intact[(cluster, KZFP)]
        
        except:

            intact_repeat_overlap_list.append([KZFP, cluster, 0, 0])

            continue

        overlap = df[KZFP]

        intact_repeat_overlap_list.append([KZFP, cluster, overlap, overlap_for_proportion])
    

intact_repeat_overlap_df = pd.DataFrame(intact_repeat_overlap_list, columns=['KZFP gene symbol', 'repeat subfamily name', 'overlap peak count to full-length copies', 'full-length copy count overlaped with peaks'])



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

# calculate  
binom_result_df = binom_result_df.sort_values(by=['KZFP gene symbol', 'log10 p value', 'ratio'], ascending=[False, False, False])
binom_result_target_df = list()
for KZFP in pd.unique(binom_result_df['KZFP gene symbol']):

    df = binom_result_df[binom_result_df['KZFP gene symbol']==KZFP].copy()

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



### write additional infomation
binom_result_target_df['copy number'] = binom_result_target_df['repeat subfamily name'].apply(lambda x:repeat_count[x])
binom_result_target_df['proportion of copy overlaped with peaks'] = binom_result_target_df['copy count overlaped with peaks'] / binom_result_target_df['copy number']
binom_result_target_df['repeat family name'] = binom_result_target_df['repeat subfamily name'].apply(lambda x:Dfam_RM_family_dict['repeat family name'][x])
binom_result_target_df['repeat class'] = binom_result_target_df['repeat subfamily name'].apply(lambda x:Dfam_RM_family_dict['repeat class'][x])
binom_result_target_df['repeat classification'] = binom_result_target_df['repeat subfamily name'].apply(lambda x:Dfam_RM_family_dict['repeat classification'][x])


# write additional infomation of full-length TEs
binom_result_target_df['full-length copy number'] = binom_result_target_df['repeat subfamily name'].apply(lambda x:repeat_count_intact[x] if x in repeat_count_intact.keys() else np.nan)
binom_result_target_df['overlap peak count to full-length copies'] = intact_repeat_overlap_df['overlap peak count to full-length copies']
binom_result_target_df['full-length copy count overlaped with peaks'] = intact_repeat_overlap_df['full-length copy count overlaped with peaks']
binom_result_target_df['proportion of full-length copy overlaped with peaks'] = binom_result_target_df['full-length copy count overlaped with peaks'] / binom_result_target_df['full-length copy number']



### write the evolutionary age of TE and KRAB-ZFPs
branch_dict = {'Vertebrata':530, 'Dipnotetrapodomorpha':413, 'Tetrapoda':352, 'Amniota':312, 'Mammalia':179.2, 'Theria':159, 'Eutheria':105, 'Boreoeutheria':96, 'Euarchontoglires': 90, 'Primatomorpha':76, 
               'Primates':74, 'Haplorrhini':63, 'Simiiformes':43.2, 'Catarrhini':29.4, 'Hominoidea':20.2, 'Hominidae':15.8, 'Homininae':9.1, 'Hominini':6.7, 'Homo sapiens':0, np.nan:None}

Age_dict = {352.0: 'Tetrapoda', 312.1: 'Amniota', 312.0: 'Amniota', 163.7:'Mammalia', 159.0: 'Theria', 105.0: 'Eutheria', 96.0: 'Boreoeutheria', 97.5: 'Boreoeutheria', 90.0: 'Euarchontoglires', 90.9: 'Euarchontoglires', 76: 'Primatomorpha', 75.9: 'Primatomorpha', 
            74.0: 'Primates', 63: 'Haplorrhini', 43.2:'Simiiformes', 43.1:'Simiiformes', 29.4: 'Catarrhini', 29.1: 'Catarrhini', 20.2:'Hominoidea', 15.8: 'Hominidae', 9.1: 'Homininae', 6.7: 'Hominini', 0: 'Homo sapiens', 'nan': np.nan}

def return_Age(x, dict):

    if x in dict.keys():

        return dict[x]
    
    else:

        return None
    

binom_result_target_df['emergence era of TE subfamily'] = binom_result_target_df['repeat subfamily name'].apply(lambda x:Dfam_RM_family_dict['branch'][x])
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

rename_columns = ['KZFP gene symbol', 'repeat subfamily name', 'ratio', 'p value', 'q value', 'log10 q value', 'normalized score', 'rank', 'criteria',
                  'overlap peak count to all copies', 'copy count overlaped with peaks', 'copy number', 'proportion of copy overlaped with peaks', 'full-length copy number', 'overlap peak count to full-length copies', 'full-length copy count overlaped with peaks', 'proportion of full-length copy overlaped with peaks',
                  'repeat family name', 'repeat class', 'repeat classification', 'emergence era of TE subfamily', 'evolutionary age of TE subfamily', 'emegence era of TE subfamily is in primate',
                  'emergence era of KZFP for analysis', 'evolutionary age of KZFP for analysis', 'evolutionary age of KZFP in Imbeault et al.', 'evolutionary age of KZFP in Tribolet-Hardy et al.', 'emegence era of KZFP is in primate']


# 処理を含める
TE_class_list = ['DNA', 'ERV/LTR', 'LINE', 'SINE', 'Retroposon']

# KZFPを制限する
KZFP_dataset_df_fil = KZFP_dataset_df[KZFP_dataset_df['classification (Genome research)']=='protein_coding']

# 制限する
condition1 = binom_result_target_df['q value'] < 0.05
condition2 = (binom_result_target_df['repeat class'].isin(TE_class_list))
condition3 = binom_result_target_df['KZFP gene symbol'].isin(KZFP_dataset_df_fil.index)
Primary = binom_result_target_df['rank']=='Primary'
Secondary = (binom_result_target_df['rank']=='Secondary') & (binom_result_target_df['log10 q value']>10) & (binom_result_target_df['ratio']>2) & (binom_result_target_df['full-length copy count overlaped with peaks']>=5) & (binom_result_target_df['proportion of full-length copy overlaped with peaks']>=0.1)
#condition3 = (binom_result_target_df['ratio']>2) & (binom_result_target_df['q value']<0.05) & (binom_result_target_df['repeat class'].isin(TE_class_list))

binom_result_target_df['criteria'] = (Primary | Secondary) & (condition1) & (condition2) & (condition3)

print(binom_result_target_df[binom_result_target_df['criteria']]['rank'].value_counts())

binom_result_target_df[rename_columns].to_csv('../data/targets/TE_targets_raw.csv', index=False)
binom_result_target_df[binom_result_target_df['criteria']][rename_columns].to_csv('../data/targets/TE_targets_for_analysis.csv', index=False)


print(time.time()-s, 'finish enrichment analysis of TEs')