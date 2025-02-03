ids=["ERR246968",
        "ERR246969",
        "ERR246970",
        "ERR246971",
        "ERR246972",
        "ERR246973",
        "ERR246974",
        "ERR246975"]
read=["1", "2"]

rule all:
    input:
        expand("results/fastqc/{id}_{read}_fastqc.{ext}", id=ids, ext=["html", "zip"], read=read),
        expand("results/fastp/{id}.{ext}", id=ids, ext={"html", "json"}),
        expand("results/fastp/{id}_{read}.fastq", id=ids, read=read),
        f"ref/{ref_name}",
        expand(f"ref/{ref_name}.{{ext}}", ext=["amb", "ann", "bwt", "pac", "sa"]),
        expand("results/map/{id}.sorted.{ext}", id=ids, ext=["bam", "bam.bai"]),
        expand("results/variants/{id}.vcf", id=ids),
        expand("results/variants/{id}.filtered.vcf", id=ids)

rule download_files:
    output:
        expand("data/{id}_{read}.fastq", id=ids, read=read)
    shell:
        """   
        accessions=()
        data_dir="data"
        accessions_list="SRR_Acc_List.txt"

        while IFS= read -r accession_id || [[ -n "$accession_id" ]];do
        accessions+=("$accession_id")
        done < "$accessions_list"

        for fq in "${{accessions[@]}}"; do
        fasterq-dump -S "$fq" -O "$data_dir"
        done
        """
rule fastqc:
    input:
        expand("data/{id}_{read}.fastq", id=ids, read=["1", "2"])
    output:
        expand("results/fastqc/{id}_{read}_fastqc.{ext}", id=ids, ext=["html", "zip"], read=read)
    shell:
        """
        data_dir="data"
        fastqc_dir="results/fastqc"
        for fq in "${{data_dir}}"/*.fastq; do
        fastqc -f fastq "$fq" -o "$fastqc_dir"
        done
        """

rule fastp:
    input:
        expand("data/{id}_{read}.fastq", id=ids, read=read)
    output:
        qc_files=expand("results/fastp/{id}.{ext}", id=ids, ext={"html", "json"}),
        trimmed_files=expand("results/fastp/{id}_{read}.fastq", id=ids, read=read)
    shell:
        """
        data_dir="data"
        fastp_dir="results/fastp"
        for fq_1 in "${{data_dir}}"/*_1.fastq;do
        r2="${{fq_1/_1/_2}}"
        r1_out="${{fastp_dir}}/$(basename "$fq_1")"
        r2_out="${{fastp_dir}}/$(basename "$r2")"
        json="${{fastp_dir}}/$(basename "$fq_1" _1.fastq).json"
        html="${{fastp_dir}}/$(basename "$fq_1" _1.fastq).html"
        fastp -i "${{fq_1}}" -o "${{r1_out}}" -I "${{r2}}" -O "${{r2_out}}" -j "${{json}}" -h "${{html}}"
        done
        """
rule download_index_ref:
    output:
        f"ref/{ref_name}",
        expand(f"ref/{ref_name}.{{ext}}", ext=["amb", "ann", "bwt", "pac", "sa"])
    params:
        ref_name=f"{ref_name}"
    shell:
        """
        set -x
        ref_dir="ref"
        ref_path="$ref_dir/{params.ref_name}"
        wget -nc -P "$ref_dir" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz
        if file $ref_path | grep -q "BGZF"; then
        echo "file correctly zipped"
        else 
        base_name="$ref_dir/$(basename "$ref_path" .gz)"
        gunzip "$ref_path"
        bgzip "$base_name"
        fi
        bwa index "$ref_path"
        """

rule map:
    input:
        trimmed_files=expand("results/fastp/{id}_{read}.fastq", id=ids, read=read),
        ref=f"ref/{ref_name}",
    output:
        expand("results/map/{id}.sorted.{ext}", id=ids, ext=["bam", "bam.bai"]),
        expand("results/logs/bam/{id}.log", id=ids)
    params:
        bwa_mem = ["-t 8"]
    shell:
        """
        fastp_dir="results/fastp"
        map_dir="results/map"
        bam_log="results/logs/bam"
        for fq in "${{fastp_dir}}"/*_1.fastq;do
        r2="${{fq/_1/_2}}"
        _out="${{map_dir}}/$(basename "$fq" _1.fastq).sam"
        bwa mem {params.bwa_mem} "{input.ref}" "${{fq}}" "${{r2}}" > "${{_out}}"
        done

        for sam_fi in "${{map_dir}}"/*.sam;do
        unsorted_bam="${{sam_fi/sam/_unsorted.bam}}"
        log="$bam_log/$(basename "$sam_fi" .sam).log"
        samtools view -Sb "${{sam_fi}}" -o "${{unsorted_bam}}" >> "${{log}}"
        done

        for un_bam in "${{map_dir}}"/*._unsorted.bam;do
        sorted_bam="${{un_bam/._unsorted/.sorted}}"
        log="$bam_log/$(basename "$un_bam" ._unsorted.bam).log"
        samtools sort "${{un_bam}}" -o "${{sorted_bam}}" && samtools index "${{sorted_bam}}" >> "${{log}}"
        done
        """

rule variant_call:
    input:
        bams=expand("results/map/{id}.sorted.{ext}", id=ids, ext=["bam", "bam.bai"]),
        ref=f"ref/{ref_name}",
        index_files=expand(f"ref/{ref_name}.{{ext}}", ext=["amb", "ann", "bwt", "pac", "sa"])
    output:
        expand("results/variants/{id}.vcf", id=ids)
    params:
        vcf_call="-Ov -f",
        bcf_tool_call="-mv  -Ov"
    shell:
        """
        map_dir="results/map"
        variants="results/variants"
        for so_bam in "${{map_dir}}"/*.sorted.bam; do
        vcf="$variants/$(basename "$so_bam" .sorted.bam).vcf"
        bcftools mpileup {params.vcf_call} {input.ref} "${{so_bam}}" |\
        bcftools call {params.bcf_tool_call} -o "${{vcf}}"
        done
        """
rule variant_filter:
    input:
        expand("results/variants/{id}.vcf", id=ids)
    output:
        expand("results/variants/{id}.filtered.vcf", id=ids)
    params:
        vcf_filter="QUAL > 20 && DP > 10"
    shell:
        """
        set -x
        variants="results/variants"
        for vcf in "${{variants}}"/*.vcf; do
        filtered_vcf="$variants/$(basename "$vcf" .vcf).filtered.vcf"
        bcftools filter -Ov -i "{params.vcf_filter}" -o "${{filtered_vcf}}" "${{vcf}}"
        done
        """
