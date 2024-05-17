import pandas as pd
import subprocess

Dfam = pd.read_csv('../data/TE/LTR7_HERVH.annotation.csv')
Dfam.to_csv('../data/TE/LTR7_HERVH.annotation.bed', index=False, header=None, sep='\t')

a = '../data/TE/LTR7_HERVH.annotation.bed'
b = '../data/roadmap/E003_15_coreMarks_hg38lift_stateno.bed'
o = '../data/overlap/LTR7_HERVH_E003_15models.bed'


cmd = 'bedtools intersect -wo -a {} -b {} > {}'.format(a, b, o)

subprocess.run(cmd, shell=True)

print('finish obtaining overlap between TE and peaks')