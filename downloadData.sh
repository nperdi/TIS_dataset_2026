#CREATE data/genome/folder

conda install samtools
conda install bedtools

mkdir -p data
mkdir -p data/GENOME

wget -O data/GENOME/hg38.chrom.size https://hgdownload.gi.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.chrom.sizes
grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y)[[:space:]]' data/GENOME/hg38.chrom.size > data/GENOME/hg38.chrom.size.clean

wget -O data/GENOME/hg38.fa.gz https://hgdownload.gi.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
gunzip data/GENOME/hg38.fa.gz
samtools faidx data/GENOME/hg38.fa

wget -O data/ensembl/Homo_sapiens.GRCh38.116.gtf.gz https://ftp.ensembl.org/pub/release-116/vertebrates/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz
wget -O data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz