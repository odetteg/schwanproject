#!/bin/bash

data_dir="data"
fastqc_dir="results/fastqc"
fastp_dir="results/fastp"
trimmed_dir="results/trimmed"
ref="ref/GCF_000001405.26_GRCh38_genomic.fna.gz"
map_dir="results/map"
bam_log="results/logs/bam"
mkdir -p "$bam_log"
variants="results/variants"

accessions=()
accessions_list="SRR_Acc_List.txt"

while IFS= read -r accession_id || [[ -n "$accession_id" ]];do
accessions+=("$accession_id")
done < "$accessions_list"

echo "${accessions[@]}"

# Step 1: Obtain the dataset
echo "downloading raw fastq reads"
for fq in "${accessions[@]}"; do
fasterq-dump -S "$fq" -O "$data_dir"
done

# Step 2: Quality check with FASTQC and FASTP

mkdir -p "$fastqc_dir"
echo "Runnning FASTQC"
for fq in "${data_dir}"/*.fastq; do
fastqc -f fastq "$fq" -o "$fastqc_dir"
done

echo "Runnning FASTP"

mkdir -p "$fastp_dir"
for fq_1 in "${data_dir}"/*_1.fastq;do

r2="${fq_1/_1/_2}"
r1_out="${fastp_dir}/$(basename "$fq_1")"
r2_out="${fastp_dir}/$(basename "$r2")"
json="${fastp_dir}/$(basename "$fq_1" _1.fastq).json"
html="${fastp_dir}/$(basename "$fq_1" _1.fastq).html"
fastp -i "${fq_1}" -o "${r1_out}" -I "${r2}" -O "${r2_out}" -j "${json}" -h "${html}"
done


# Step 3: indexing index genome and mapping the reads

wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz

bwa index "$ref"


mkdir -p "$map_dir"
for fq in "${fastp_dir}"/*_1.fastq;do
r2="${fq/_1/_2}"
_out="${map_dir}/$(basename "$fq" _1.fastq).sam"
bwa mem -t 8 "${ref}" "${fq}" "${r2}" > "${_out}"
done

# Step 4: Converting sam to bam, sorting, and indexing.
echo "converting to bam, sorting, and indexing"


for sam_fi in "${map_dir}"/*.sam;do
unsorted_bam="${sam_fi/sam/_unsorted.bam}"
log="$(basename "$sam_fi" .sam).bam.log"
samtools view -Sb "${sam_fi}" -o "${unsorted_bam}" >> "${log}"
done


for un_bam in "${map_dir}"/*._unsorted.bam;do
sorted_bam="${un_bam/._unsorted/.sorted}"
log="${bam_log}/$(basename "$un_bam" .unsorted.sam).bam.log"
samtools sort "${un_bam}" -o "${sorted_bam}" && samtools index "${sorted_bam}" >> "${log}"
done


# Step 5: variant calling and filtering

echo "Running variant calling..."
mkdir -p "$variants"

for so_bam in "${map_dir}"/*.sorted.bam; do
vcf="${variants}/$(basename "$so_bam" .sorted.bam).vcf"
bcftools mpileup -Ov -f "${ref}" "${so_bam}" | bcftools call -mv  -Ov -o "${vcf}"
done


for vcf in "${variants}"/*.vcf; do
filtered_vcf="${variants}/$(basename "$vcf" .vcf).filtered.vcf"
bcftools filter -Ov -i "QUAL > 20 && DP > 10" -o "${filtered_vcf}" "${vcf}"
done

echo "variant calling completed..."
