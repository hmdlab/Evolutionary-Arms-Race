import subprocess
import pandas as pd

# obtain overlap between TE and peaks.
TE_list = ['TE', 'provirus']
TE_file_list = ['TE_annotation', 'provirus_annotation']

peak_list = ['KZFP', 'TRIM28', 'TF']
peak_file_list = ['KZFP_peaks.bed', 'TRIM28_peaks.bed', 'TF_peaks.bed']

for TE, TE_file in zip(TE_list, TE_file_list):

    csv = pd.read_csv('../data/TE/{}.csv'.format(TE_file))
    bed = csv.to_csv('../data/TE/{}.bed'.format(TE_file), index=False, header=None, sep='\t')

    for peak, peak_file in zip(peak_list, peak_file_list):

        a = '../data/TE/{}.bed'.format(TE_file)
        b = '../data/ChIP-seq/{}'.format(peak_file)
        o = '../data/overlap/{}_{}_overlap.bed'.format(TE, peak)
        cmd = 'bedtools intersect -wo -a {} -b {} > {}'.format(a, b, o)

        subprocess.run(cmd, shell=True)


print('finish obtaining overlap between TE and peaks')