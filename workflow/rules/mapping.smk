import pandas as pd
import numpy as np
import math
import os
import sys


rule split_reads:
    input:
        reads=lambda wc: config["reads"][wc.sm],
    output:
        reads=temp(scatter.split("temp/reads/{{sm}}/{scatteritem}.fq.gz")),
    conda:
        "../envs/env.yml"
    resources:
        mem=2,
        hrs=8,
        load=100,  # seeting a high load here so that only a few can run at once
    params:
        sdir=SDIR,
        unzipped=scatter.split("temp/reads/{{sm}}/{scatteritem}.fq"),
    log:
        "logs/split_reads/{sm}.log",
    benchmark:
        "benchmarks/split_reads/{sm}.tbl"
    threads: 8
    priority: 10
    shell:
        """
        if [[ {input.reads} =~ \.(fasta|fasta.gz|fa|fa.gz|fastq|fastq.gz|fq|fq.gz)$ ]]; then 
            cat {input.reads} \
                | seqtk seq -F '#' \
                | rustybam fastq-split {output.reads} 
        elif [[ {input.reads} =~ \.(bam|cram|sam|sam.gz)$ ]]; then 
            samtools fasta -@ {threads} {input.reads} \
                | seqtk seq -F '#' \
                | rustybam fastq-split {output.reads} 
        fi 
        """


rule mrsfast_index:
    input:
        ref=config.get("masked_ref", rules.masked_reference.output.fasta),
    output:
        index=config.get("masked_ref", rules.masked_reference.output.fasta) + ".index",
    conda:
        "../envs/env.yml"
    log:
        "logs/mrsfast/index.{sample}.log",
    resources:
        mem=8,
        hrs=24,
    threads: 1
    shell:
        """
        mrsfast --index {input.ref}
        """

# removed these DG June 9, 2023 after 10 to 20 tries at debugging the
#        conda environment.  Even snakemake was repeatedly crashing
#        (not just the jobs).  So I just replaced the conda
#        environment with a build of mrsfast and it just worked.

#        mem=lambda wildcards, attempt, threads: 4 * attempt * threads,
#    conda:
#        "../envs/env.yml"
#        total_mem=lambda wildcards, attempt, threads: 4 * attempt * threads - 2,

rule mrsfast_alignment:
    input:
        reads="temp/reads/{sm}/{scatteritem}.fq.gz",
        index=rules.mrsfast_index.output.index,
        ref=config.get("masked_ref", rules.masked_reference.output.fasta),
    output:
        sam=temp("temp/mrsfast/{sample}/{sm}/{scatteritem}.sam.gz"),
    resources:
        total_mem=50,
        mem=50,
        hrs=2,
        load=1,
    log:
        "logs/mrsfast/{sample}/{sm}/{scatteritem}.log",
    benchmark:
        "benchmarks/{sample}/mrsfast/{sm}/{scatteritem}.tbl"
    threads: 4
    priority: 20
    shell:
        """
        echo "about to execute: extract-from-fastq36.py --in {input.reads} | mrsfast --search {input.ref} --seq /dev/stdin --disable-nohits --mem {resources.total_mem} --threads {threads} -e 2 --outcomp -o $(dirname {output.sam})/{wildcards.scatteritem} > {log} 2>&1"

        export PATH=/home/hsiehph/shared/software/packages/mrsfast/sfu-compbio-mrsfast-cf8e678:$PATH && extract-from-fastq36.py --in {input.reads} \
            | mrsfast --search {input.ref} --seq /dev/stdin \
                --disable-nohits --mem {resources.total_mem} --threads {threads} \
                -e 2 --outcomp \
                -o $(dirname {output.sam})/{wildcards.scatteritem} \
            > {log} 2>&1
        """


rule mrsfast_sort:
    input:
        sam=rules.mrsfast_alignment.output.sam,
    output:
        bam=temp("temp/mrsfast/{sample}/{sm}/{scatteritem}.bam"),
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/mrsfast/sort/{sm}/{scatteritem}_sort.log",
    benchmark:
        "benchmarks/{sample}/sort_bam/{sm}/{scatteritem}.tbl"
    resources:
        mem=4,
        hrs=24,
        load=1,
    threads: 2
    priority: 30
    shell:
        """
        zcat {input.sam} \
            | samtools view -b - \
            | samtools sort -@ {threads} \
             -T {resources.tmpdir} -m 2G \
             -o {output.bam} -
        """


rule merged_mrsfast_bam:
    input:
        bams=gather.split("temp/mrsfast/{{sample}}/{{sm}}/{scatteritem}.bam"),
    output:
        merged=temp("results/{sample}/mapping/{sm}_merged.out.gz"),
    conda:
        "../envs/env.yml"
    resources:
        mem=4,
        hrs=24,
    benchmark:
        "benchmarks/{sample}/merge_mrsfast/{sm}.tbl"
    log:
        "logs/mrsfast/{sample}/{sm}.merged.log",
    threads: 4
    priority: 40
    shell:
        """
        samtools merge -@ {threads} - {input.bams} -u \
            | samtools view -h - \
            | awk 'BEGIN {{OFS="\\t"}} {{print $1,$2,$3,$4,$5}}' \
            | pigz -p {threads} \
        > {output.merged}
        """


rule compress_mrsfast_further:
    input:
        sam=rules.merged_mrsfast_bam.output.merged,
    output:
        comp="results/{sample}/mapping/{sm}_merged_comp.out.gz",
#    conda:
#        "../envs/env.yml"
    resources:
        mem=2,
        hrs=24,
        load=25,
    benchmark:
        "benchmarks/{sample}/comp_mrsfast/{sm}.tbl"
    log:
        "logs/mrsfast/{sample}/{sm}.merged_comp.log",
    script:
        "../scripts/compress_mrsfast.py"
