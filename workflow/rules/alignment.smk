rule align:
    input:
        r1 = "results/trimmed/{sample}_1_val_1.fq.gz",
        r2 = "results/trimmed/{sample}_2_val_2.fq.gz"
    output:
        bam = "results/aligned/{sample}.sorted.bam"
    params:
        index = config["genome"]["index"]
    threads: 8
    resources:
        mem_mb = 16000,
        runtime = 240,
        partition = "general",
        account = "r00086"
    shell:
        """
        bowtie2 -x {params.index} -1 {input.r1} -2 {input.r2} --threads {threads} --rg-id {wildcards.sample} --rg "SM:{wildcards.sample}" --rg "PL:ILLUMINA" --rg "LB:ATAC" | samtools sort -@ {threads} -o {output.bam}
        """

rule filter:
    input:
        bam = "results/aligned/{sample}.sorted.bam"
    output:
        bam = "results/filtered/{sample}.filtered.bam"
    params:
        mapq = config["alignment"]["mapq"]
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 120,
        partition = "general",
        account = "r00086"
    shell:
        """
        samtools view -b -q {params.mapq} -f 2 -@ {threads} {input.bam} | samtools view -b -@ {threads} -o {output.bam} -e 'rname != "chrM"'
        """

rule dedup:
    input:
        bam = "results/filtered/{sample}.filtered.bam"
    output:
        bam     = "results/dedup/{sample}.dedup.bam",
        metrics = "results/dedup/{sample}.metrics.txt"
    threads: 4
    resources:
        mem_mb    = 32000,
        runtime   = 120,
        partition = "general",
        account   = "r00086"
    shell:
        """
        samtools index {input.bam}
        picard MarkDuplicates \
            I={input.bam} \
            O={output.bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=true \
            VALIDATION_STRINGENCY=LENIENT
        samtools index {output.bam}
        """
