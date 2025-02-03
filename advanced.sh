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
# In this workflow, we can refactory step 1 by simply storing the accession IDs
# In a external file then read it in our workflow and handle it as a list
# We will store all accession IDs in the file SRR.Acc_List.txt

# Working the .txt file with accession IDs
# First, set an empty list like this

accessions=()
accessions_list="SRR_Acc_List.txt"

# Next read the .txt file and add the accession IDs to our empty list using a loop

while IFS= read -r accession_id || [[ -n "$accession_id" ]];do
accessions+=("$accession_id")
done < "$accessions_list"

# Let us unpack what is happenig. We are using the while loop, and invoking the read command. The 
# Read command process the file line by line and assigns whatever is in each line to a variable "accession_id"
# we are settig the IFS (internal field seperator) to empty.
# This prevents read from trimming leading/trailing whitespace or splitting the input into multiple fields if spaces or tabs exist in the line.
# Without setting IFS=, read would split lines on whitespace by default.
# The -r option tells read not to interpret backslashes (\) as escape characters. 
# This ensures that backslashes in the input are treated as literal characters.
# The [[ -n "$accession_id" ]] part ensures that the loop processes the last line of the file correctly, 
# even if the file does not end with a newline character.

# We can then verify that all out accession ids have been added in the accession list using this command
echo "${accessions[@]}"

With that, we can then invoke a for loop to download the files
echo "downloading raw fastq reads"
for fq in "${accessions[@]}"; do
fasterq-dump -S "$fq" -O "$data_dir"
done

# Step 2: Running fastqc and fastp
# Here, instead of processing the files line by line, we will use a for loop and end up with a short relatively robust code
mkdir -p "$fastqc_dir"
echo "Runnning FASTQC"
for fq in "${data_dir}"/*.fastq; do
fastqc -f fastq "$fq" -o "$fastqc_dir"
done

echo "Runnning FASTP"

# Just like we did for the fastqc, we will use a for loop to process the raw reads. FASTP will do both quality check and filtering. 
# The report will show us what our files looked like before and after filtering
mkdir -p "$fastp_dir"
for fq_1 in "${data_dir}"/*_1.fastq;do
# We need to create some files for us to work with fastp
# fastp needs output file for both read 1 and 2, and a json file for each sample
r2="${fq_1/_1/_2}"
r1_out="${fastp_dir}/$(basename "$fq_1")"
r2_out="${fastp_dir}/$(basename "$r2")"
json="${fastp_dir}/$(basename "$fq_1" _1.fastq).json"
html="${fastp_dir}/$(basename "$fq_1" _1.fastq).html"
fastp -i "${fq_1}" -o "${r1_out}" -I "${r2}" -O "${r2_out}" -j "${json}" -h "${html}"
done

# Let us go throuhg what each command is doing. -i informs fastp about our read 1 input file, -o about our read 1 output file. 
# This applies to the -I and -O for read 2 as well. We then request fastp to give us a json file (to work well with multiqc) and a html file 
# Which will have our quality metrics
# Just a note, we are using the command basename here to extract the file name without the full path. 
# Running basename on something like  /results/fastp/ERR246968_1.fastq will give you ERR246968_1.fastq. We can then modify the path anyhow we want.z

# Step 3: indexing index genome and mapping the reads

# The third step is often to index the reference genome, then mapp our reads to the indexed genome.
# We will download the reference genome using wget then invoke the bwa command
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz
# the -nc simply tells the wget command not to download it again if the file exist

bwa index "$ref"

# We will use for loop again to process the files. This way, you will appreciate why it is often better to optimize workflows
# We will specifically work on the files in the fastp directory becuase they have been trimmed

mkdir -p "$map_dir"
for fq in "${fastp_dir}"/*_1.fastq;do
r2="${fq/_1/_2}"
_out="${map_dir}/$(basename "$fq" _1.fastq).sam"
bwa mem -t 8 "${ref}" "${fq}" "${r2}" > "${_out}"
done

# Step 4: Converting sam to bam, sorting, and indexing.
echo "converting to bam, sorting, and indexing"
# Once we have our sam files, the next thing is to convert them to bam, sort and index them. This makes the files a little bit lightweight
# and easier to search. 
# Also, many other samtool commands we will use downstream will require bam files. This conversion is therefore important.

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

# The final step in our worklow, for now will be to call the variants. Here, we will use the bcftools mpileup then call the variants with bcftool call
# Note that previous implementations used samtools mpileup before bcftools call but this has been depricated.
# The bcftools mpile up simply looks at the reference genome, marks the position of each base the "piles up each read"
echo "Running variant calling..."
mkdir -p "$variants"

for so_bam in "${map_dir}"/*.sorted.bam; do
vcf="${variants}/$(basename "$so_bam" .sorted.bam).vcf"
bcftools mpileup -Ov -f "${ref}" "${so_bam}" | bcftools call -mv  -Ov -o "${vcf}"
done


# As the final step, we will simply filter our variants. If you open any of the vcf files, you wil find a lot of information, 
# among them are columns: #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO. In the last colum, there are details such as 
# DP=1;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,1,0;MQ=60
# All these store information that we can use to filter the variants. In our case, we want to filter and keep only variants with
# A base call quality of > 20 (if you read this paper https://www.nature.com/articles/s41525-021-00227-3, you will realize that QUAL > 20)
# will give you high quality variants. We will also filter where DP is > 10

for vcf in "${variants}"/*.vcf; do
filtered_vcf="${variants}/$(basename "$vcf" .vcf).filtered.vcf"
bcftools filter -Ov -i "QUAL > 20 && DP > 10" -o "${filtered_vcf}" "${vcf}"
done

echo "variant calling completed..."

# So based on what we have done above, we have developed a "workflow" that takes up to 8 files, and processing them in five steps.
# Based on this simply workflow alone, we can process as many files as we want, using the same commands and with some level of ease. 
# You can simply run the entire workflow (see main.sh) and it will process all the files. 