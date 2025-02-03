ids = ["SAMPLE1", "SAMPLE2", "SAMPLE3"] # put in your sample IDs
read = ["1", "2"]
ref_name = "REFERENCE.fasta"

rule all:
    input:
        expand("results/fastqc/{id}_{read}_fastqc.{ext}", id=ids, ext=["html", "zip"], read=read),
        expand("results/fastp/{id}.{ext}", id=ids, ext=["html", "json"]),
        expand("results/multiqc/multiqc_report.html"),
        expand("results/map/{id}.sorted.{ext}", id=ids, ext=["bam", "bam.bai"]),
        expand("results/variants/{id}.filtered.vcf", id=ids)

rule download_data:
    output:
        expand("data/{id}_{read}.fastq", id=ids, read=read)
    log:
        "logs/download.log"
    threads: 4
    shell:
        """
        # Use fastq-dump or wget to download sequencing data
        # Fill in appropriate command
        """

rule fastqc:
    input:
        # Define your input
    output:
        expand("results/fastqc/{id}_{read}_fastqc.{ext}", id=ids, ext=["html", "zip"], read=read)
    log:
        "logs/fastqc.log"
    threads: 2
    resources:
        mem_mb=2000
    shell:
        """
        # Run FastQC for quality control
        # Fill in the command
        """

rule fastp:
    input:
        expand("data/{id}_{read}.fastq", id=ids, read=read)
    output:
        # Define your output
    log:
        "logs/fastp.log"
    threads: 4
    params:
        trim_length=50
    shell:
        """
        # Perform adapter trimming using fastp. Incorporate the params
        # Fill in the command
        """

rule multiqc:
    input:
        expand("results/fastqc/{id}_{read}_fastqc.zip", id=ids, read=read),
        expand("results/fastp/{id}.json", id=ids)
    output:
        "results/multiqc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    shell:
        """
        # Run MultiQC to summarize quality control reports
        # Fill in the command
        """

rule map_reads:
    input:
        expand("results/fastp/{id}_{read}.fastq", id=ids, read=read),
        "ref/" + ref_name
    output:
        expand("results/map/{id}.sorted.{ext}", id=ids, ext=["bam", "bam.bai"])
    log:
        "logs/mapping.log"
    threads: 8
    resources:
        mem_mb=8000
    shell:
        """
        # Align reads to reference genome
        # Fill in the command
        """

rule variant_calling:
    input:
        expand("results/map/{id}.sorted.bam", id=ids),
        "ref/" + ref_name
    output:
        expand("results/variants/{id}.vcf", id=ids)
    log:
        "logs/variant_calling.log"
    params:
        min_qual=20
    shell:
        """
        # Perform variant calling
        # Fill in the command
        """

rule filter_variants:
    input:
        expand("results/variants/{id}.vcf", id=ids)
    output:
        expand("results/variants/{id}.filtered.vcf", id=ids)
    log:
        "logs/variant_filtering.log"
    params:
        filter_criteria="QUAL > 30 && DP > 10"
    shell:
        """
        # Apply variant quality filtering (use bcftools)
        # Fill in the command
        """
