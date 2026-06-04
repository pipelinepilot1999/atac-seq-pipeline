#!/bin/bash
# ============================================================
# HOMER motif analysis on differential peaks
# Run after differential_analysis.R
# Requires: module load homer/4.11.1
# Usage: bash scripts/run_homer.sh
# ============================================================

module load homer/4.11.1

PROJDIR="/N/slate/sappagan/atac_project"
PREPARSED="${PROJDIR}/genome/homer_preparsed"

mkdir -p ${PROJDIR}/results/motifs/cancer \
         ${PROJDIR}/results/motifs/normal \
         ${PREPARSED}

findMotifsGenome.pl \
    ${PROJDIR}/results/cancer_open_peaks.bed \
    hg38 \
    ${PROJDIR}/results/motifs/cancer/ \
    -size 200 -p 4 \
    -preparsedDir ${PREPARSED}

findMotifsGenome.pl \
    ${PROJDIR}/results/normal_open_peaks.bed \
    hg38 \
    ${PROJDIR}/results/motifs/normal/ \
    -size 200 -p 4 \
    -preparsedDir ${PREPARSED}

echo "HOMER complete"
