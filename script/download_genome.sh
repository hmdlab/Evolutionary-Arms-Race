mkdir -p ../data/UCSC
mkdir -p ../data/gencode

wget -P ../data/UCSC http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
wget -P ../data/gencode  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz
wget -P ../data/roadmap https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E003_15_coreMarks_hg38lift_stateno.bed.gz

unpigz -f ../data/UCSC/hg38.fa.gz
unpigz -f ../data/gencode/gencode.v40.annotation.gtf.gz
unpigz -f ../data/roadmap/E003_15_coreMarks_hg38lift_stateno.bed.gz