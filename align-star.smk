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
    input: expand("star/pe/{sample}/Aligned.out.sam", sample = units.index)

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

rule star_index:
    input:
        fasta="/fs/papilionanthe/assembly/sup_bc/flye240222/06-repeats/vmj3.rm.fa",
    output:
        directory("vmj3"),
    message:
        "Testing STAR index"
    threads: 1
    params:
        extra="",
    log:
        "logs/star_index_vmj.log",
    wrapper:
        "v1.7.0/bio/star/index"

rule star_pe_multi:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1="trimmed/pe/{sample}_1.fq.gz",
        # paired end reads needs to be ordered so each item in the two lists match
        fq2= "trimmed/pe/{sample}_2.fq.gz",  #optional
        index = "vmj3"
    output:
        # see STAR manual for additional output files
        sam="star/pe/{sample}/Aligned.out.sam",
        log="star/pe/{sample}/Log.out",
    log:
        "logs/star/pe/{sample}.log",
    params:
        # path to STAR reference genome index
        idx="vmj3",
        # optional parameters
        extra="",
    threads: 8
    wrapper:
        "v1.7.0/bio/star/align"