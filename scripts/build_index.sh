#!/bin/bash
# Build bowtie2 index from hg38 genome
# Run on LOGIN NODE — takes ~35 minutes
# Usage: bash scripts/build_index.sh

GENOMEDIR="/N/slate/sappagan/atac_project/genome"

mkdir -p $GENOMEDIR/bowtie2_index

bowtie2-build \
    --threads 8 \
    $GENOMEDIR/hg38.fa \
    $GENOMEDIR/bowtie2_index/hg38

echo "Done — index files:"
ls -lh $GENOMEDIR/bowtie2_index/
