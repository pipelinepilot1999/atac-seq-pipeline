rule fastqc:
    input:
        r1 = lambda wc: config["samples"][wc.sample]["R1"],
        r2 = lambda wc: config["samples"][wc.sample]["R2"]
    output:
        r1_html = "results/qc/fastqc/{sample}_1_fastqc.html",
        r2_html = "results/qc/fastqc/{sample}_2_fastqc.html"
    threads: 4
    resources:
        mem_mb    = 8000,
        runtime   = 60,
        partition = "general",
        account   = "r00086"
    shell:
        """
        fastqc {input.r1} {input.r2} \
            --outdir results/qc/fastqc/ \
            --threads {threads}
        """

rule trim:
    input:
        r1 = lambda wc: config["samples"][wc.sample]["R1"],
        r2 = lambda wc: config["samples"][wc.sample]["R2"]
    output:
        r1 = "results/trimmed/{sample}_1_val_1.fq.gz",
        r2 = "results/trimmed/{sample}_2_val_2.fq.gz"
    params:
        adapter = config["trimming"]["adapter"]
    threads: 8
    resources:
        mem_mb    = 8000,
        runtime   = 120,
        partition = "general",
        account   = "r00086"
    shell:
        """
        trim_galore --paired \
            --adapter {params.adapter} \
            --cores {threads} \
            --output_dir results/trimmed/ \
            {input.r1} {input.r2}
        """
