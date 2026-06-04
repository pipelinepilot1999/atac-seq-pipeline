rule bigwig:
    input:
        bam = "results/dedup/{sample}.dedup.bam"
    output:
        bw = "results/bigwigs/{sample}.bw"
    threads: 4
    resources:
        mem_mb    = 16000,
        runtime   = 60,
        partition = "general",
        account   = "r00086"
    shell:
        """
        bamCoverage \
            -b {input.bam} \
            -o {output.bw} \
            --normalizeUsing RPGC \
            --effectiveGenomeSize 2913022398 \
            --extendReads \
            --binSize 10 \
            -p {threads}
        """

rule qc_metrics:
    input:
        metrics = expand("results/dedup/{sample}.metrics.txt", sample=SAMPLES),
        peaks   = expand("results/peaks/{sample}_peaks_filtered.bed", sample=SAMPLES)
    output:
        summary = "results/qc/qc_summary.tsv"
    threads: 1
    resources:
        mem_mb    = 4000,
        runtime   = 30,
        partition = "general",
        account   = "r00086"
    script:
        "../../scripts/qc_metrics.py"
