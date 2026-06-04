rule call_peaks:
    input:
        bam = "results/dedup/{sample}.dedup.bam"
    output:
        peaks = "results/peaks/{sample}_peaks.narrowPeak"
    params:
        genome_size = config["peaks"]["genome_size"],
        qvalue      = config["peaks"]["qvalue"],
        name        = "{sample}"
    threads: 4
    resources:
        mem_mb    = 16000,
        runtime   = 60,
        partition = "general",
        account   = "r00086"
    shell:
        """
        macs2 callpeak \
            -t {input.bam} \
            -f BAMPE \
            -n {params.name} \
            -g {params.genome_size} \
            -q {params.qvalue} \
            --nomodel \
            --shift -75 \
            --extsize 150 \
            --keep-dup all \
            --outdir results/peaks/
        """

rule blacklist_filter:
    input:
        peaks     = "results/peaks/{sample}_peaks.narrowPeak",
        blacklist = config["genome"]["blacklist"]
    output:
        peaks = "results/peaks/{sample}_peaks_filtered.bed"
    threads: 1
    resources:
        mem_mb    = 8000,
        runtime   = 60,
        partition = "general",
        account   = "r00086"
    shell:
        """
        bedtools subtract \
            -a {input.peaks} \
            -b {input.blacklist} \
            > {output.peaks}
        """
