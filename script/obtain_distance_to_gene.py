import pandas as pd
import re
import warnings
import time
warnings.simplefilter('ignore')

s = time.time()
print('start obtaining distance')

### process gtf

# load gtf
gtf = pd.read_table('../data/gencode/gencode.v40.annotation.gtf', skiprows=5, header=None)

# process gtf
gtf_gene = gtf[gtf[2]=='gene']
gtf_gene.columns = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 8]

gene_feature_list = [[], [], []]
for row in gtf_gene[8]:

    value_list = row.split('; ')
    
    for i, value in enumerate(value_list[:3]):

        v_list = value.split(' ')
        gene_feature_list[i].append(v_list[1][1:-1])


gtf_gene['gene_id version'] = gene_feature_list[0]
gtf_gene['gene_type'] = gene_feature_list[1]
gtf_gene['gene_name'] = gene_feature_list[2]

regex_gene_id = re.compile(r'(.+)\.')
gene_id_list = list()
for row in gtf_gene['gene_id version']:
    id = regex_gene_id.search(row).group(1)
    gene_id_list.append(id)

gtf_gene['gene_id'] = gene_id_list

# obtain TSS
gtf_gene_TSS = gtf_gene.copy()
TSS_list = list()
for start, end, strand in gtf_gene_TSS[['start', 'end', 'strand']].values:

    if strand == '+':

        TSS_list.append(start)
    
    else:

        TSS_list.append(end)

gtf_gene_TSS['start'] = TSS_list
gtf_gene_TSS['end'] = gtf_gene_TSS['start'] + 1

# output TSS annotation
columns = ['chr', 'start', 'end', 'score', 'strand', 'frame', 'source', 'feature', 'gene_id', 'gene_id version', 'gene_type', 'gene_name']
gtf_gene_TSS[columns].to_csv('../data/gencode/gencode.v40.annotation_TSS.csv', index=False)



### search nearby genes

# load annotation data
family = 'LTR7_HERVH'
Dfam_RM = pd.read_csv('../data/TE/{}.annotation.csv'.format(family))
gtf_gene_TSS = pd.read_csv('../data/gencode/gencode.v40.annotation_TSS.csv')

# calucurate distance
result_list = list()
for i, repeat in enumerate(Dfam_RM.values):

    chr, start, end = repeat[:3]

    for gene in gtf_gene_TSS.values:

        TSS_chr, TSS_start = gene[:2] 

        if chr != TSS_chr:

            continue

        if start <= TSS_start and TSS_start < end:

            distance = 0
        
        else:

            distance = min(abs(TSS_start-start), abs(TSS_start-end))
        

        if distance <= 50000:

            result_list.append(list(repeat) + list(gene) + [distance])

print(len(result_list))


# output nearby genes data
result_df = pd.DataFrame(result_list, columns=Dfam_RM.columns.tolist()+gtf_gene_TSS.columns.tolist()+['distance'])
result_df.to_csv('../data/overlap/{}.annotation_window_gencode_TSS.bed'.format(family))


print(time.time(), 'finish obtaining distance')
