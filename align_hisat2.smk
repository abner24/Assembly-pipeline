import glob

import pandas as pd
from snakemake.utils import validate

samples = pd.read_table("samples.tsv").set_index("sample", drop=False)

units = pd.read_table("units.tsv", dtype=str).set_index("sample", drop=False)

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = units.loc[wildcards.sample, ["fq1", "fq2"]].dropna()
    print(unpack(fastqs.values))
    if len(fastqs) == 2:
        return [fastqs.fq1, fastqs.fq2]
    return {"r1": fastqs.fq1}


stringtie = "/fs/plant_tools/stringtie-2.1.6/stringtie"

rule all:
    input: expand('stringtie/{sample}.gtf', sample = units.index)

rule fastp_pe:
    input:
        sample=get_fastq,
    output:
        trimmed=["trimmed/pe/{sample}_1.fq.gz", "trimmed/pe/{sample}_2.fq.gz"],
        # Unpaired reads separately
        unpaired1="trimmed/pe/{sample}.u1.fq.gz",
        unpaired2="trimmed/pe/{sample}.u2.fq.gz",
        # or in a single file
#        unpaired="trimmed/pe/{sample}.singletons.fastq",
        merged="trimmed/pe/{sample}.merged.fq.gz",
        failed="trimmed/pe/{sample}.failed.fq.gz",
        html="report/pe/{sample}.html",
        json="report/pe/{sample}.json"
    log:
        "logs/fastp/pe/{sample}.log"
    params:
        adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        extra="--merge"
    threads: 2
    wrapper:
        "v1.5.0/bio/fastp"

rule hisat2_align:
    input:
      reads=["trimmed/pe/{sample}_1.fq.gz", "trimmed/pe/{sample}_2.fq.gz"]
    output:
      "mapped/{sample}.bam"
    log:
        "logs/hisat2_align_{sample}.log"
    params:
      extra="",
      idx="/fs/abner/data/coriander/CorSat-hisat2",
    threads: 4
    wrapper:
      "v1.5.0/bio/hisat2/align"

rule sort:
    input:
        unSrt = "mapped/{sample}.bam"
    output:
        srt = "sort/{sample}.bam"
    log:
        "logs/sort/samtools_{sample}.log"
    threads: 4
    shell:
        """
        (samtools sort \
            -@ {threads} \
            -o {output.srt} \
            {input.unSrt} )2> {log}
        """

rule stringtie:
    input:
        bam = "sort/{sample}.bam"
    output:
        gtf = "stringtie/{sample}.gtf"
    params:
        st = stringtie,
        label = "{sample}"
    log:
        "logs/stringtie/stringtie_{sample}.log"
    threads: 4
    shell:
        "({params.st} "
            "-p {threads} "
            "-o {output.gtf} "
            "-l {params.label} "
            "{input.bam} )"
        "2> {log}"