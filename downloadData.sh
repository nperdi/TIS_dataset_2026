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
