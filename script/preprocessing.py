import os
import subprocess
import time

s = time.time()

if __name__ == "__main__":

    os.makedirs('../data/gencode', exist_ok=True)
    os.makedirs('../data/UCSC', exist_ok=True)
    os.makedirs('../data/roadmap', exist_ok=True)
    os.makedirs('../data/overlap', exist_ok=True)
    os.makedirs('../data/targets', exist_ok=True)
    os.makedirs('../data/screening', exist_ok=True)
    os.makedirs('../data/MSA/conserved', exist_ok=True)
    os.makedirs('../data/motif/fimo', exist_ok=True)
    os.makedirs('../data/motif/peak_sequences', exist_ok=True)

    # download public data
    subprocess.run('bash download_genome.sh', shell=True)

    # obtain overlap
    subprocess.run('python obtain_overlap_between_TE_and_peaks.py', shell=True)

    # enrichment analysis
    subprocess.run('python enrichment_analysis_TE.py', shell=True)
    subprocess.run('python enrichment_analysis_provirus.py', shell=True)

    # screening
    subprocess.run('python screening_for_evasion.py', shell=True)

    # search binding site
    subprocess.run('python search_binding_sites_TE.py', shell=True)
    subprocess.run('python search_binding_sites_provirus.py', shell=True)

    # output each TE family data
    subprocess.run('python process_various_data_for_each_TE_family.py', shell=True)
    subprocess.run('python process_various_data_for_MER11.py', shell=True)
    subprocess.run('python process_various_data_for_each_provirus.py', shell=True)

    # obtain overlap between LTR7_HERVH and roadmap
    subprocess.run('python obtain_overlap_between_LTR7_and_roadmap.py', shell=True)

    # obtain distance between LTR7_HERVH and TSS
    subprocess.run('python obtain_distance_to_gene.py', shell=True)

    print(time.time()-s, 'All preprocessing done')