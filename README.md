# Reproducible ATAC-seq Workflow for Differential Chromatin Accessibility in Breast Cancer (MCF7 vs Normal)

## Overview

A modular Snakemake pipeline for processing bulk ATAC-seq data from raw FASTQ files through peak calling, followed by differential chromatin accessibility analysis comparing MCF7 breast cancer cells to normal breast epithelium. Built on IU HPC (Quartz) with SLURM job scheduling.


## Biological Question

Which regions of the genome are differentially accessible between MCF7 breast cancer cells and normal breast epithelial tissue, and what transcription factors and biological pathways are associated with these changes?

## Key Results

| Metric | Value |
|--------|-------|
| Consensus peaks | 139,921 |
| Significantly differential peaks (padj < 0.05, \|log2FC\| > 1) | 52,429 |
| More open in cancer (MCF7) | ~22,000 |
| More open in normal | ~30,000 |

### Top Motifs

- **Cancer-open regions**: FOXA1, CTCF, BORIS — consistent with ER+ breast cancer biology. FOXA1 is the known pioneer factor that opens chromatin for estrogen receptor binding in MCF7.
- **Normal-open regions**: AP-1 family (Fra1, Fosl2, Fos), NF1 — tissue maintenance and homeostasis transcription factors.

### GO Enrichment

- **Cancer-open**: cell growth regulation, GTPase signaling, neuron projection development (hijacked developmental pathways)
- **Normal-open**: cell junction assembly, extracellular matrix organization, synapse organization (tissue integrity programs lost in cancer)

## Dataset

| Sample | SRA Accession | Condition | Source |
|--------|--------------|-----------|--------|
| MCF7_rep1 | SRR14103347 | Cancer | ENCODE ENCSR422SUG |
| MCF7_rep2 | SRR14103348 | Cancer | ENCODE ENCSR422SUG |
| Normal_donor_A | SRR10388801 | Normal | ENCODE ENCSR846ZBX |
| Normal_donor_B | SRR14305402 | Normal | ENCODE ENCSR654UYP |

All samples are paired-end 101bp Illumina reads aligned to the hg38 human genome assembly.

### Known Limitations

- MCF7 is a cell line (homogeneous) while normal samples are primary tissue (heterogeneous cell types). This introduces variance asymmetry that DESeq2's statistical model may not fully account for.
- Normal samples are from different donors (not true biological replicates of the same individual), contributing to higher inter-sample variability.
- These limitations are documented as design constraints, not errors.

## Pipeline Architecture

```
Raw FASTQ
  │
  ├── FastQC (quality control)
  │
  ├── Trim Galore (adapter removal)
  │     Nextera adapter: CTGTCTCTTATACACATCT
  │
  ├── Bowtie2 (alignment to hg38)
  │     with read group tags
  │
  ├── samtools (filtering)
  │     remove chrM, MAPQ < 30, improper pairs
  │
  ├── Picard MarkDuplicates (deduplication)
  │
  ├── MACS2 (peak calling)
  │     --shift -75 --extsize 150 (ENCODE standard)
  │
  ├── bedtools (blacklist filtering)
  │     ENCODE blacklist v2
  │
  ├── deeptools bamCoverage (bigwig generation)
  │
  └── QC metrics summary
```

### Tn5 Shift Correction

Tn5 shift correction is applied within MACS2 using `--shift -75 --extsize 150` (ENCODE ATAC-seq standard), rather than a separate base-pair-resolution +4/-5 shift step, since this project performs peak calling and differential accessibility analysis — not TF footprinting, which would require single-base precision.

## Project Structure

```
atac_project/
├── workflow/
│   ├── Snakefile                  # Master pipeline
│   └── rules/
│       ├── qc.smk                 # FastQC, Trim Galore
│       ├── alignment.smk          # Bowtie2, samtools filter, Picard dedup
│       ├── peaks.smk              # MACS2 peak calling, blacklist filter
│       └── coverage.smk           # BigWig generation, QC metrics
├── scripts/
│   ├── download_fastq.sh          # Download raw FASTQ from SRA
│   ├── rename_compress.sh         # Rename and compress FASTQ files
│   ├── download_genome.sh         # Download hg38 reference + blacklist
│   ├── build_index.sh             # Build bowtie2 index
│   ├── qc_metrics.py              # Parse Picard metrics + count peaks
│   ├── prepare_counts.sh          # Build consensus peaks + count matrix
│   └── run_homer.sh               # HOMER motif analysis
├── analysis/
│   └── differential_analysis.R    # DESeq2 + ChIPseeker + GO enrichment
├── config/
│   └── config.yaml                # All pipeline parameters
├── envs/
│   └── atac-pipeline.yml          # Conda environment (full reproducibility)
├── results/                       # All pipeline and analysis outputs
├── raw_data/                      # Raw FASTQ files (not in repo)
├── genome/                        # Reference genome files (not in repo)
└── logs/                          # SLURM and tool logs
```

## Reproducibility

### Environment Setup

```bash
# Create the conda environment
conda env create -f envs/atac-pipeline.yml

# Activate
conda activate atac-pipeline
```

### Tool Versions

| Tool | Version | Purpose |
|------|---------|---------|
| FastQC | 0.12.1 | Read quality control |
| Trim Galore | 2.2.0 | Adapter and quality trimming |
| Bowtie2 | 2.5.5 | Read alignment |
| samtools | 1.23.1 | BAM filtering and indexing |
| Picard | 3.4.0 | PCR duplicate removal |
| MACS2 | 2.2.9.1 | Peak calling |
| deeptools | 3.5.6 | BigWig generation |
| bedtools | 2.31.1 | Blacklist filtering, consensus peaks |
| Snakemake | 7.32.4 | Pipeline orchestration |
| HOMER | 4.11.1 | Motif enrichment analysis |
| DESeq2 | (R/Bioconductor) | Differential accessibility |
| ChIPseeker | (R/Bioconductor) | Peak annotation |
| clusterProfiler | (R/Bioconductor) | GO enrichment |

### Running the Full Pipeline

```bash
# Step 1: Download and prepare data (login node)
bash scripts/download_fastq.sh
bash scripts/rename_compress.sh
bash scripts/download_genome.sh
bash scripts/build_index.sh

# Step 2: Run Snakemake pipeline (submits to SLURM)
snakemake \
    --snakefile workflow/Snakefile \
    --configfile config/config.yaml \
    --jobs 8 \
    --cluster "sbatch --partition=general --account=r00086 \
               --mem={resources.mem_mb}M --time={resources.runtime} \
               --cpus-per-task={threads} \
               --output=logs/slurm_%j.out \
               --error=logs/slurm_%j.err" \
    --default-resources mem_mb=8000 runtime=60

# Step 3: Prepare count matrix
bash scripts/prepare_counts.sh

# Step 4: Differential analysis
Rscript analysis/differential_analysis.R

# Step 5: Motif analysis
bash scripts/run_homer.sh
```

## QC Summary

| Sample | Duplication Rate | Peak Count |
|--------|-----------------|------------|
| MCF7_rep1 | 18.3% | 90,084 |
| MCF7_rep2 | 16.7% | 99,492 |
| Normal_donor_A | 13.3% | 56,197 |
| Normal_donor_B | 12.7% | 44,829 |

All samples pass ENCODE quality thresholds:
- Duplication rates < 40%
- Peak counts within 50K-200K range
- 55% adapter content confirms short ATAC fragments (Tn5 tagmentation read-through)
- Per base sequence content FAIL expected (Tn5 insertion bias)

## Output Files

| File | Description |
|------|-------------|
| `results/volcano_plot.pdf` | Differential accessibility volcano plot |
| `results/deseq2_results.csv` | Full DESeq2 results for all 139,921 peaks |
| `results/significant_peaks.csv` | 52,429 significant differential peaks |
| `results/GO_cancer_dotplot.pdf` | GO enrichment for cancer-open regions |
| `results/GO_normal_dotplot.pdf` | GO enrichment for normal-open regions |
| `results/motifs/cancer/` | HOMER motif results for cancer-open peaks |
| `results/motifs/normal/` | HOMER motif results for normal-open peaks |
| `results/qc/qc_summary.tsv` | Per-sample QC metrics |
| `results/bigwigs/*.bw` | Genome-wide coverage tracks for visualization |

## References

- Buenrostro, J.D. et al. (2013). Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins and nucleosome position. *Nature Methods*, 10(12), 1213-1218.
- ENCODE Project Consortium. ENCODE ATAC-seq pipeline. https://github.com/ENCODE-DCC/atac-seq-pipeline
- Ewels, P.A. et al. nf-core/atacseq. https://github.com/nf-core/atacseq
- Amemiya, H.M., Kundaje, A., & Boyle, A.P. (2019). The ENCODE Blacklist: Identification of Problematic Regions of the Genome. *Scientific Reports*, 9, 9354.

## Author

Sathvik Sai Appagana — MS Bioinformatics


## License

MIT License
