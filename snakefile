# Snakefile for tools used in tb-profiler

# Define the configuration file
configfile: "config.yaml"

# Load configuration values from yaml file
samples = config["samples"]
tri_in = config["tri_in"]
tri_out = config["tri_out"]
threads = config["threads"]
leading = config["leading"]
trailing = config["trailing"]
sliding = config["slidingwindow"]
minlen = config["minlen"]
map_in = config["map_in"]
map_out = config["map_out"]
map_refin = config["map_refin"]
sort_in = config["sort_in"]
sort_out = config["sort_out"]
index_in = config["index_in"]
index_out = config["index_out"]
gatk_in = config["gatk_in"]
gatk_refin = config["gatk_refin"]
gatk_out = config["gatk_out"]
norm_in = config["norm_in"]
norm_refin = config["norm_refin"]
norm_out = config["norm_out"]


# Rule all
rule all:
    input:
        expand(f"{gatk_out}/{{samples}}.vcf", samples=samples),
        expand(f"{index_out}/{{samples}}_sorted.bam.bai", samples=samples),
        expand(f"{gatk_out}/{{samples}}.vcf", samples=samples),
        expand(f"{norm_out}/{{samples}}_normalize.vcf.gz", samples=samples) 
       

# Rule for Trimmomatic
rule trimmomatic:
    input:
        r1 = f"{tri_in}/{{samples}}_1.fastq.gz",
        r2 = f"{tri_in}/{{samples}}_2.fastq.gz"
    output:
        r3 = f"{tri_out}/{{samples}}_forward_paired.fastq.gz",
        r4 = f"{tri_out}/{{samples}}_forward_unpaired.fastq.gz",
        r5 = f"{tri_out}/{{samples}}_reversed_paired.fastq.gz",
        r6 = f"{tri_out}/{{samples}}_reversed_unpaired.fastq.gz"
    log:
        "log/{samples}_trim.log"
    conda:
        "env/environment.yaml"
    threads: threads
    shell:
        """
        set -euo pipefail
        trimmomatic PE -threads {threads} {input.r1} {input.r2} \
        {output.r3} {output.r4} {output.r5} {output.r6} \
        LEADING:{leading} \
        TRAILING:{trailing} \
        SLIDINGWINDOW:{sliding} \
        MINLEN:{minlen} > {log} 2>&1
        """

# Rule for Mapping
rule mapping:
    input:
        b1 = f"{map_in}/{{samples}}_forward_paired.fastq.gz",
        b2 = f"{map_in}/{{samples}}_reversed_paired.fastq.gz",
        b3 = f"{map_refin}/sequence.fasta"
    output:
        b4 = f"{map_out}/{{samples}}.bam"
    log:
        "log/{samples}_mapping.log"
    conda:
        "env/environment.yaml"
    shell:
        """
        set -euo pipefail
        bwa mem -R '@RG\\tID:{wildcards.samples}\\tSM:{wildcards.samples}\\tPL:ILLUMINA' -M -T 50 \
        {input.b3} {input.b1} {input.b2} | samtools view -Sb - > {output.b4} 2> {log}
        """

# Rule for Sorting
rule sort:
    input:
        c1 = f"{sort_in}/{{samples}}.bam"
    output:
        c2 = f"{sort_out}/{{samples}}_sorted.bam"
    log:
        "log/{samples}_sort.log"
    conda:
        "env/environment.yaml"
    threads: threads
    shell:
        """
        set -euo pipefail
        samtools sort -@ {threads} {input.c1} -o {output.c2} 2> {log}
        """

# Rule for Indexing
rule index:
    input:
        d1 = f"{index_in}/{{samples}}_sorted.bam"
    output:
        d2 = f"{index_out}/{{samples}}_sorted.bam.bai"
    log:
        "log/{samples}_index.log"
    conda:
        "env/environment.yaml"
    shell:
        """
        set -euo pipefail
        samtools index {input.d1} 2> {log}
        """

# Rule for GATK
rule gatk:
    input:
        g1 = f"{gatk_in}/{{samples}}_sorted.bam",
        g2 = f"{gatk_refin}/sequence.fasta"
    output:
        g3 = f"{gatk_out}/{{samples}}.vcf"
    log:
        "log/{samples}_gatk.log"
    conda:
        "env/environment.yaml"
    shell:
        """
        set -euo pipefail
        gatk HaplotypeCaller \
        -A StrandBiasBySample \
        -R {input.g2} \
        -I {input.g1} \
        -O {output.g3} 2> {log}
        """

# rule for bcftools norm

rule norm:
    input:
        h1 = f"{norm_in}/{{samples}}.vcf",
        h2 = f"{norm_refin}/sequence.fasta"
    output:
        h3 = f"{norm_out}/{{samples}}_normalize.vcf.gz"
    log:
        "log/{samples}_norm.log"
    conda:
        "env/environment.yaml"
    shell:
        """
        set -euo pipefail
        bcftools norm \
        -f {input.h2} \
        -Oz -o {output.h3} \
        {input.h1} &> {log}
        tabix -p vcf {output.h3}
        """

