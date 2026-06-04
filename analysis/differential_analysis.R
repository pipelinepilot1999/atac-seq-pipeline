# ============================================================
# Differential Chromatin Accessibility Analysis
# MCF7 (cancer) vs Normal breast epithelium
# Run after Snakemake pipeline completes
# Usage: Rscript analysis/differential_analysis.R
# ============================================================

library(DESeq2)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GenomicRanges)

setwd("/N/slate/sappagan/atac_project")

# ── Step 1: Read count matrix ────────────────────────────────
counts <- read.delim("results/peaks/consensus_counts_header.tsv")
peak_info <- counts[, 1:3]
count_matrix <- as.matrix(counts[, 4:7])

# ── Step 2: Sample metadata ─────────────────────────────────
col_data <- data.frame(
    condition = c("cancer", "cancer", "normal", "normal"),
    row.names = c("MCF7_rep1", "MCF7_rep2", "Normal_donor_A", "Normal_donor_B")
)

# ── Step 3: DESeq2 ──────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = col_data,
    design = ~ condition
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "cancer", "normal"))

# ── Step 4: Prepare results ─────────────────────────────────
res_df <- as.data.frame(res)
res_df$chr   <- peak_info$chr
res_df$start <- peak_info$start
res_df$end   <- peak_info$end

res_df$color <- ifelse(
    res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
    ifelse(res_df$log2FoldChange > 0, "red", "blue"),
    "grey"
)

# ── Step 5: Volcano plot ────────────────────────────────────
pdf("results/volcano_plot.pdf", width = 8, height = 6)
plot(
    res_df$log2FoldChange,
    -log10(res_df$padj),
    col = res_df$color,
    pch = 20, cex = 0.5,
    xlab = "log2 Fold Change (Cancer vs Normal)",
    ylab = "-log10(adjusted p-value)",
    main = "Differential Chromatin Accessibility\nMCF7 vs Normal Breast"
)
abline(h = -log10(0.05), lty = 2)
abline(v = c(-1, 1), lty = 2)
legend("topright",
    legend = c("More open in Cancer", "More open in Normal", "Not significant"),
    col = c("red", "blue", "grey"), pch = 20
)
dev.off()

# ── Step 6: Save results ────────────────────────────────────
write.csv(res_df, "results/deseq2_results.csv", row.names = FALSE)
sig_peaks <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]
write.csv(sig_peaks, "results/significant_peaks.csv", row.names = FALSE)
cat("Significant peaks:", nrow(sig_peaks), "\n")

# ── Step 7: Peak annotation ─────────────────────────────────
sig_gr <- GRanges(
    seqnames = sig_peaks$chr,
    ranges = IRanges(start = sig_peaks$start, end = sig_peaks$end)
)
sig_gr$log2FoldChange <- sig_peaks$log2FoldChange
sig_gr$padj <- sig_peaks$padj

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- annotatePeak(sig_gr, TxDb = txdb, annoDb = "org.Hs.eg.db", tssRegion = c(-3000, 3000))
anno_df <- as.data.frame(peakAnno)

pdf("results/peak_annotation_pie.pdf", width = 8, height = 8)
plotAnnoPie(peakAnno)
dev.off()

# ── Step 8: GO enrichment ───────────────────────────────────
cancer_open <- anno_df[anno_df$log2FoldChange > 1, ]
normal_open <- anno_df[anno_df$log2FoldChange < -1, ]

cancer_go <- enrichGO(
    gene = unique(cancer_open$SYMBOL), OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05
)
normal_go <- enrichGO(
    gene = unique(normal_open$SYMBOL), OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05
)

write.csv(as.data.frame(cancer_go), "results/GO_cancer_open.csv", row.names = FALSE)
write.csv(as.data.frame(normal_go), "results/GO_normal_open.csv", row.names = FALSE)

pdf("results/GO_cancer_dotplot.pdf", width = 10, height = 8)
dotplot(cancer_go, showCategory = 15, title = "GO: Regions More Open in Cancer")
dev.off()

pdf("results/GO_normal_dotplot.pdf", width = 10, height = 8)
dotplot(normal_go, showCategory = 15, title = "GO: Regions More Open in Normal")
dev.off()

# ── Step 9: BED files for HOMER ─────────────────────────────
cancer_bed <- sig_peaks[sig_peaks$log2FoldChange > 1, c("chr", "start", "end")]
write.table(cancer_bed, "results/cancer_open_peaks.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

normal_bed <- sig_peaks[sig_peaks$log2FoldChange < -1, c("chr", "start", "end")]
write.table(normal_bed, "results/normal_open_peaks.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Done. Cancer peaks:", nrow(cancer_bed), "Normal peaks:", nrow(normal_bed), "\n")
