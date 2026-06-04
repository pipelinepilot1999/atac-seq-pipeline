#!/bin/bash
# Download reference genome files for ATAC-seq pipeline
# Run on LOGIN NODE only
# Usage: bash scripts/download_genome.sh

GENOMEDIR="/N/slate/sappagan/atac_project/genome"

mkdir -p $GENOMEDIR

# hg38 genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz \
    -O $GENOMEDIR/hg38.fa.gz

# Decompress
gunzip $GENOMEDIR/hg38.fa.gz

# Chromosome sizes
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes \
    -O $GENOMEDIR/hg38.chrom.sizes

# ENCODE blacklist v2
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz \
    -O $GENOMEDIR/hg38-blacklist.v2.bed.gz

gunzip $GENOMEDIR/hg38-blacklist.v2.bed.gz

echo "Done — genome files:"
ls -lh $GENOMEDIR
