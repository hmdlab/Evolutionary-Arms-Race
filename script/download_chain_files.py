import pandas as pd
import subprocess

chain_files_df = pd.read_csv('../Supplementary Table/Supplementary Table 2_Liftover species.csv')

for i, genome in enumerate(chain_files_df['genome name']):

    genome_name = genome[0].upper()+genome[1:]

    path = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38To{}.over.chain.gz'.format(genome_name)
    cmd = 'wget -P ../data/UCSC/chain_files {}'.format(path)
    subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)

    unzip_cmd = 'unpigz -f ../data/UCSC/chain_files/{}'.format('hg38To{}.over.chain.gz'.format(genome_name))
    subprocess.run(unzip_cmd, shell=True, stdout=subprocess.DEVNULL)

    print(i, genome)