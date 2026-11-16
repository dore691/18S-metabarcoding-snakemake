# 1. VARIABLES

configfile: "config/config.yaml"

SAMPLES = config["samples"]
READS = config["reads"]

# 2. OBJECTIFS DU PIPELINE

rule all:
    input:
        expand("results/fastqc/{sample}_{read}_001_fastqc.html",
               sample=SAMPLES, read=READS),
        expand("results/trimmed/{sample}_{read}_trimmed_001.fastq.gz",
               sample=SAMPLES, read=READS),
        expand("results/fastqc_trimmed/{sample}_{read}_trimmed_001_fastqc.html",
               sample=SAMPLES, read=READS),

        "results/trimmed/results_dada2/seqtab.nochim.rds",
        "results/trimmed/results_dada2/ASV_length_distribution.pdf",
        "results/trimmed/results_dada2/track.txt",
        "results/trimmed/results_dada2/retention_track.pdf"

# 3. CONTROLE QUALITE

rule fastqc:
    input:
        "data/{sample}_{read}_001.fastq.gz"
    output:
        html="results/fastqc/{sample}_{read}_001_fastqc.html",
        zip="results/fastqc/{sample}_{read}_001_fastqc.zip"
    threads: 1
    conda: "envs/fastqc.yaml"
    shell:
        "fastqc --threads {threads} --outdir results/fastqc {input}"

# 4. RETRAIT DES AMORCES

rule cutadapt:
    input:
        r1="data/{sample}_R1_001.fastq.gz",
        r2="data/{sample}_R2_001.fastq.gz"
    output:
        r1="results/trimmed/{sample}_R1_trimmed_001.fastq.gz",
        r2="results/trimmed/{sample}_R2_trimmed_001.fastq.gz"
    params:
        primer_f=config["primers"]["forward"],
        primer_r=config["primers"]["reverse"]
    log:
        "logs/{sample}_cutadapt.log"
    threads: 1
    conda:
        "envs/cutadapt.yaml"
    shell:
        """
        cutadapt \
            -j {threads} \
            -g {params.primer_f} \
            -G {params.primer_r} \
            --discard-untrimmed \
            -o {output.r1} \
            -p {output.r2} \
            {input.r1} {input.r2} > {log} 2>&1
        """

# 5. CONTROLE QUALITE APRES TRIMMING

rule fastqc_trimmed:
    input:
        "results/trimmed/{sample}_{read}_trimmed_001.fastq.gz"
    output:
        html="results/fastqc_trimmed/{sample}_{read}_trimmed_001_fastqc.html",
        zip="results/fastqc_trimmed/{sample}_{read}_trimmed_001_fastqc.zip"
    threads: 1
    conda: "envs/fastqc.yaml"
    shell:
        "fastqc --threads {threads} --outdir results/fastqc_trimmed {input}"


# 6. ANALYSE DADA2 SANS TAXONOMIE
rule dada2:
    input:
        r1=expand(
            "results/trimmed/{sample}_R1_trimmed_001.fastq.gz",
            sample=SAMPLES
        ),
        r2=expand(
            "results/trimmed/{sample}_R2_trimmed_001.fastq.gz",
            sample=SAMPLES
        )

    output:
        asv="results/trimmed/results_dada2/seqtab.nochim.rds",
        length_plot="results/trimmed/results_dada2/ASV_length_distribution.pdf",
        track="results/trimmed/results_dada2/track.txt",
        retention_plot="results/trimmed/results_dada2/retention_track.pdf"

    params:
        truncLen=config["dada2"]["truncLen"],
        maxN=config["dada2"]["maxN"],
        maxEE=config["dada2"]["maxEE"],
        truncQ=config["dada2"]["truncQ"],
        rm_phix=config["dada2"]["rm_phix"],
        compress=config["dada2"]["compress"]

    threads: config["dada2"]["threads"]

    conda:
        "envs/dada2.yaml"

    log:
        "logs/dada2.log"

    script:
        "scripts/script_dada2_nochim.R"



# 7. OBJECTIF PRETRAITEMENT

rule preprocessing:
    input:
        expand(
            "results/fastqc/{sample}_{read}_001_fastqc.html",
            sample=SAMPLES, read=READS
        ),
        expand(
            "results/trimmed/{sample}_{read}_trimmed_001.fastq.gz",
            sample=SAMPLES, read=READS
        ),
        expand(
            "results/fastqc_trimmed/{sample}_{read}_trimmed_001_fastqc.html",
            sample=SAMPLES, read=READS
        )


# 8. ASSIGNATION TAXONOMIQUE PR2

rule taxonomy:
    input:
        asv="results/trimmed/results_dada2/seqtab.nochim.rds",
        database=config["taxonomy"]["database"]

    output:
        taxonomy="results/trimmed/results_dada2/taxa_PR2.rds"

    params:
        taxLevels=config["taxonomy"]["taxLevels"],
        minBoot=config["taxonomy"]["minBoot"]

    threads: config["taxonomy"]["threads"]

    conda:
        "envs/dada2.yaml"

    log:
        "logs/taxonomy_PR2.log"

    script:
        "scripts/script_taxonomie_PR2.R"


# 9. OBJECTIF PIPELINE COMPLET

rule complete:
    input:
        rules.all.input,
        "results/trimmed/results_dada2/taxa_PR2.rds"
