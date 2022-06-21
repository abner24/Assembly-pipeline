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
    input: 
        expand("star/pe/{sample}/Aligned.out.sam", sample = units.index),
        "results/rsem.genes.results"

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

rule prepare_reference:
    input:
        # reference FASTA with either the entire genome or transcript sequences
        reference_genome="/fs/papilionanthe/assembly/sup_bc/flye240222/06-repeats/vmj3.rm.fa",
    output:
        # one of the index files created and used by RSEM (required)
        seq="index/reference.seq",
        # RSEM produces a number of other files which may optionally be specified as output; these may be provided so that snakemake is aware of them, but the wrapper doesn't do anything with this information other than to verify that the file path prefixes match that of output.seq.
        # for example,
        grp="index/reference.grp",
        ti="index/reference.ti",
    params:
        # optional additional parameters, for example,
        #extra="--gtf annotations.gtf",
        # if building the index against a reference transcript set
        extra="--gtf /fs/papilionanthe/assembly/sup_bc/flye240222/annotate/pmj.gtf",
    log:
        "logs/rsem/prepare-reference.log",
    wrapper:
        "v1.7.0/bio/rsem/prepare-reference"

rule calculate_expression:
    input:
        # input.bam or input.fq_one must be specified (and if input.fq_one, optionally input.fq_two if paired-end)
        # an aligned to transcriptome BAM
        bam="star/pe/{sample}/Aligned.out.sam",
        # one of the index files created by rsem-prepare-reference; the file suffix is stripped and passed on to rsem
        reference="index/reference.seq",
    output:
        # genes_results must end in .genes.results; this suffix is stripped and passed to rsem as an output name prefix
        # this file contains per-gene quantification data for the sample
        genes_results="results/rsem/{sample}.genes.results",
        # isoforms_results must end in .isoforms.results and otherwise have the same prefix as genes_results
        # this file contains per-transcript quantification data for the sample
        isoforms_results="results/rsem/{sample}.isoforms.results",
    params:
        # optional, specify if sequencing is paired-end
        paired_end=True,
        # additional optional parameters to pass to rsem, for example,
        extra="--seed 42",
    log:
        "logs/rsem/calculate_expression/a.log",
    wrapper:
        "v1.7.0/bio/rsem/calculate-expression"

rule rsem_generate_data_matrix:
    input: expand("results/rsem/{sample}.gene.results", sample = unit.index)
        # one or more expression files created by rsem-calculate-expression
    output:
        # a tsv containing each sample in the input as a column
        "results/rsem.genes.results",
    params:
        # optional additional parameters
        extra="",
    log:
        "logs/rsem/generate_data_matrix.log",
    wrapper:
        "v1.7.0/bio/rsem/generate-data-matrix"