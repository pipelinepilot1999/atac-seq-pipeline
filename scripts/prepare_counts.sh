#!/bin/bash
# ============================================================
# Prepare consensus peaks and count matrix for DESeq2
# Run after Snakemake pipeline completes
# Usage: bash scripts/prepare_counts.sh
# ============================================================

PROJDIR="/N/slate/sappagan/atac_project"

# Step 1: Build consensus peak set
cat ${PROJDIR}/results/peaks/*_peaks_filtered.bed | \
    cut -f1-3 | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - > ${PROJDIR}/results/peaks/consensus_peaks.bed

echo "Consensus peaks: $(wc -l < ${PROJDIR}/results/peaks/consensus_peaks.bed)"

# Step 2: Count reads per peak per sample
bedtools multicov \
    -bams ${PROJDIR}/results/dedup/MCF7_rep1.dedup.bam \
          ${PROJDIR}/results/dedup/MCF7_rep2.dedup.bam \
          ${PROJDIR}/results/dedup/Normal_donor_A.dedup.bam \
          ${PROJDIR}/results/dedup/Normal_donor_B.dedup.bam \
    -bed ${PROJDIR}/results/peaks/consensus_peaks.bed \
    > ${PROJDIR}/results/peaks/consensus_counts.tsv

# Step 3: Add header
echo -e "chr\tstart\tend\tMCF7_rep1\tMCF7_rep2\tNormal_donor_A\tNormal_donor_B" | \
    cat - ${PROJDIR}/results/peaks/consensus_counts.tsv \
    > ${PROJDIR}/results/peaks/consensus_counts_header.tsv

echo "Count matrix: $(wc -l < ${PROJDIR}/results/peaks/consensus_counts_header.tsv) rows"
echo "Done"
