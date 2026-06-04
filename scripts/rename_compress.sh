#!/bin/bash

OUTDIR="/N/slate/sappagan/atac_project/raw_data"

# Rename
mv $OUTDIR/SRR14103347_1.fastq $OUTDIR/MCF7_rep1_1.fastq
mv $OUTDIR/SRR14103347_2.fastq $OUTDIR/MCF7_rep1_2.fastq
mv $OUTDIR/SRR14103348_1.fastq $OUTDIR/MCF7_rep2_1.fastq
mv $OUTDIR/SRR14103348_2.fastq $OUTDIR/MCF7_rep2_2.fastq
mv $OUTDIR/SRR10388801_1.fastq $OUTDIR/Normal_donor_A_1.fastq
mv $OUTDIR/SRR10388801_2.fastq $OUTDIR/Normal_donor_A_2.fastq
mv $OUTDIR/SRR14305402_1.fastq $OUTDIR/Normal_donor_B_1.fastq
mv $OUTDIR/SRR14305402_2.fastq $OUTDIR/Normal_donor_B_2.fastq

echo "Renaming done"

# Compress
gzip $OUTDIR/*.fastq &
wait

echo "Done — final files:"
ls -lh $OUTDIR
