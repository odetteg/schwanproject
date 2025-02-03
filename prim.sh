#!/bin/bash

# Prelimiary step: Seting up directorys
data_dir="data"
fastqc_dir="results/fastqc"
fastp_dir="results/fastp"
trimmed_dir="results/trimmed"
ref="ref/GCF_000001405.26_GRCh38_genomic.fna.gz"
bam_dir="results/map"
variants="results/variants"
# mkdir -p "$data_dir"
# Step 1: Obtain the dataset
# Our task here is simple, just to obtain the files and make sure they are in the data dir.
# We will use fastq-dump to obtain the raw reads, so I urge you to install the SRA toolkit. I will write a simple 
# command to download the files here but would recommend that you look at the advanced version

echo "downloading raw fastq files:"

# fasterq-dump -S ERR246968 -O "$data_dir" # This command is simple. We are simply calling the fasterq-dump command,
# telling it to download the file with the run ID ERR246968 and split it for us, give as both r1 and r2. 
# We the ask it to store it in the data dir. We can do the same for all other remaining SRA IDs. 
# I will be downloading and processing all eight files but urge you to may be just do this on two because of the 
# Memory issues.

fasterq-dump -S ERR246969 -O "$data_dir"
fasterq-dump -S ERR246970 -O "$data_dir"
fasterq-dump -S ERR246971 -O "$data_dir"
fasterq-dump -S ERR246972 -O "$data_dir"
fasterq-dump -S ERR246973 -O "$data_dir"
fasterq-dump -S ERR246974 -O "$data_dir"
fasterq-dump -S ERR246975 -O "$data_dir"

echo "download complete and all raw reads store in "${data_dir}" directory"

# Step 2: Quality check with FASTQC and FASTP
# Both FASTQC and FASTP requires us to supply a list of files. But for this simple workflow, we will only supply
# one file at a time

echo "Runnning FASTQC"
# We need to have the directory existing, fastqc will not create it for us.
mkdir -p "$fastqc_dir"
fastqc data/ERR246968_1.fastq data/ERR246968_2.fastq -o "$fastqc_dir"
fastqc data/ERR246969_1.fastq data/ERR246969_2.fastq -o "$fastqc_dir"
fastqc data/ERR246970_1.fastq data/ERR246970_2.fastq -o "$fastqc_dir"
fastqc data/ERR246971_1.fastq data/ERR246971_2.fastq -o "$fastqc_dir"
fastqc data/ERR246972_1.fastq data/ERR246972_2.fastq -o "$fastqc_dir"
fastqc data/ERR246973_1.fastq data/ERR246973_2.fastq -o "$fastqc_dir"
fastqc data/ERR246974_1.fastq data/ERR246974_2.fastq -o "$fastqc_dir"
fastqc data/ERR246975_1.fastq data/ERR246975_2.fastq -o "$fastqc_dir"

echo "Running FASTP"
# For fastp, the situation is slightly different. Because it does both quality check and trimming, we will need to provide both input and out files
mkdir -p "$fastp_dir"
fastp \
-i data/ERR246968_1.fastq \
-o results/fastp/ERR246968_1.fastq \
-I data/ERR246968_2.fastq \
-O results/fastp/ERR246968_2.fastq \
-j results/fastp/ERR246968.json  \
-h results/fastp/ERR246968.html 

Let us go throuhg what each command is doing. -i informs fastp about our read 1 input file, -o about our read 1 output file. 
This applies to the -I and -O for read 2 as well. We then request fastp to give us a json file (to work well with multiqc) and a html file 
Which will have our quality metrics

fastp \
-i data/ERR246969_1.fastq \
-o results/fastp/ERR246969_1.fastq \
-I data/ERR246969_2.fastq \
-O results/fastp/ERR246969_2.fastq \
-j results/fastp/ERR246969.json  \
-h results/fastp/ERR246969.html 

fastp \
-i data/ERR246970_1.fastq \
-o results/fastp/ERR246970_1.fastq \
-I data/ERR246970_2.fastq \
-O results/fastp/ERR246970_2.fastq \
-j results/fastp/ERR246970.json  \
-h results/fastp/ERR246970.html 

fastp \
-i data/ERR246971_1.fastq \
-o results/fastp/ERR246971_1.fastq \
-I data/ERR246971_2.fastq \
-O results/fastp/ERR246971_2.fastq \
-j results/fastp/ERR246971.json  \
-h results/fastp/ERR246971.html 

fastp \
-i data/ERR246972_1.fastq \
-o results/fastp/ERR246972_1.fastq \
-I data/ERR246972_2.fastq \
-O results/fastp/ERR246972_2.fastq \
-j results/fastp/ERR246972.json  \
-h results/fastp/ERR246972.html 

fastp \
-i data/ERR246973_1.fastq \
-o results/fastp/ERR246973_1.fastq \
-I data/ERR246973_2.fastq \
-O results/fastp/ERR246973_2.fastq \
-j results/fastp/ERR246973.json  \
-h results/fastp/ERR246973.html 

fastp \
-i data/ERR246974_1.fastq \
-o results/fastp/ERR246974_1.fastq \
-I data/ERR246974_2.fastq \
-O results/fastp/ERR246974_2.fastq \
-j results/fastp/ERR246974.json  \
-h results/fastp/ERR246974.html 

fastp \
-i data/ERR246975_1.fastq \
-o results/fastp/ERR246975_1.fastq \
-I data/ERR246975_2.fastq \
-O results/fastp/ERR246975_2.fastq \
-j reults/fastp/ERR246975.json  \
-h results/fastp/ERR246975.html 

# Step 3: indexing index genome and mapping the reads

# The third step is often to index the reference genome, then mapp our reads to the indexed genome.
# We will downloa the reference genome using wget then invoke the bwa command
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz
# the -nc simply tells the wget command not to download it again if the file exist

bwa index "$ref"

# To map our files, we will use bwa mem. But you can use bowtie. If you use bwa to index the ref file, just use bwa mem to map the reads
# Do the same if you used bowtie.

# Step 3.1
mkdir -p "$bam_dir"
bwa mem \
-t 8 \
"$ref" \
results/fastp/ERR246975_1.fastq \
results/fastp/ERR246975_2.fastq \
> results/bam/ERR246975.sam  


bwa mem \
-t 8 \
"$ref" \
results/fastp/ERR246974_1.fastq \
results/fastp/ERR246974_2.fastq \
> results/bam/ERR246974.sam

bwa mem \
-t 8 \
"$ref" \
results/fastp/ERR246973_1.fastq \
results/fastp/ERR246973_2.fastq \
> results/bam/ERR246973.sam

bwa mem \
-t 8 \
"$ref" \
results/fastp/ERR246972_1.fastq \
results/fastp/ERR246972_2.fastq \
> results/bam/ERR246972.sam

bwa mem \
-t 8 \
"$ref" \
results/fastp/ERR246971_1.fastq \
results/fastp/ERR246971_2.fastq \
> results/bam/ERR246971.sam

bwa mem \
-t 8 \
"$ref" \
results/fastp/ERR246970_1.fastq \
results/fastp/ERR246970_2.fastq \
> results/bam/ERR246970.sam

bwa mem \
-t 8 \
"$ref" \
results/fastp/ERR246969_1.fastq \
results/fastp/ERR246969_2.fastq \
> results/bam/ERR246969.sam

bwa mem \
-t 8 \
"$ref" \
results/fastp/ERR246968_1.fastq \
results/fastp/ERR246968_2.fastq \
> results/bam/ERR246968.sam

echo "aligning reads with bwa me finished"

# Step 4.0: Converting sam to bam, sorting, and indexing.
# Once we have our sam files, the next thing is to convert them to bam, sort and index them. This makes the files a little bit lightweight
# and easier to search. 
# Also, many other samtool commands we will use downstread will require bam files. This conversion is therefore important.
 
samtools view \
-Sb \
results/map/ERR246969.sam \
-o results/map/unsorted_ERR246969.bam

samtools view \
-Sb \
results/map/ERR246970.sam \
-o results/map/unsorted_ERR246970.bam

samtools view \
-Sb \
results/map/ERR246971.sam \
-o results/map/unsorted_ERR2471.bam

samtools view \
-Sb \
results/map/ERR246972.sam \
-o results/map/unsorted_ERR246972.bam

samtools view \
-Sb \
results/map/ERR246973.sam \
-o results/map/unsorted_ERR246973.bam

samtools view \
-Sb \
results/map/ERR246974.sam \
-o results/map/unsorted_ERR246974.bam

samtools view \
-Sb \
results/map/ERR246975.sam \
-o results/map/unsorted_ERR246975.bam

samtools view \
-Sb \
results/map/ERR246968.sam \
-o results/map/unsorted_ERR246968.bam

samtools view \
-Sb \
results/map/ERR246969.sam \
-o results/map/unsorted_ERR246969.bam

# Step 4.1: Sorting and indexing
echo "sorting and indexing"
samtools sort \
results/map/unsorted_ERR246969.bam \
-o results/map/sorted_ERR246969.bam \
&& samtools index results/map/sorted_ERR246969.bam  

samtools sort \
results/map/unsorted_ERR246970.bam \
-o results/map/sorted_ERR246970.bam \
&& samtools index results/map/sorted_ERR246970.bam  

samtools sort \
results/map/unsorted_ERR246971.bam \
-o results/map/sorted_ERR246971.bam \
&& samtools index results/map/sorted_ERR246971.bam  

samtools sort \
results/map/unsorted_ERR246972.bam \
-o results/map/sorted_ERR246972.bam \
&& samtools index results/map/sorted_ERR246972.bam  

samtools sort \
results/map/unsorted_ERR246973.bam \
-o results/map/sorted_ERR246973.bam \
&& samtools index results/map/sorted_ERR246973.bam  

samtools sort \
results/map/unsorted_ERR246974.bam \
-o results/map/sorted_ERR246974.bam \
&& samtools index results/map/sorted_ERR246974.bam  

samtools sort \
results/map/unsorted_ERR246975.bam \
-o results/map/sorted_ERR246975.bam \
&& samtools index results/map/sorted_ERR24675.bam  

# Step 5: variant calling and filtering

# The final step in our worklow, for now will be to call the variants. Here, we will use the bcftools mpileup then call the variants with bcftool call
# Note that previous implementations used samtools mpileup before bcftools call but this has been depicated.
# The bcftools mpile up simply looks at the reference genome, marks the position of each base the "piles up each read"

echo "Running variant calling..."

mkdir -p "$variants"
bcftools mpileup \
-Ov \
-f "$ref" \
results/map/ERR246975.sorted.bam \
| bcftools call -mv -Ov -o results/variants/ERR246975.vcf \

bcftools mpileup \
-Ov \
-f "$ref" \
results/map/ERR246974.sorted.bam \
| bcftools call -mv -Ov -o results/variants/ERR246974.vcf \

bcftools mpileup \
-Ov \
-f "$ref" \
results/map/ERR246973.sorted.bam \
| bcftools call -mv -Ov -o results/variants/ERR246973.vcf \

bcftools mpileup \
-Ov \
-f "$ref" \
results/map/ERR246972.sorted.bam \
| bcftools call -mv -Ov -o results/variants/ERR246972.vcf \

bcftools mpileup \
-Ov \
-f "$ref" \
results/map/ERR246971.sorted.bam \
| bcftools call -mv -Ov -o results/variants/ERR246971.vcf \



bcftools mpileup \
-Ov \
-f "$ref" \
results/map/ERR246970.sorted.bam \
| bcftools call -mv -Ov -o results/variants/ERR246970.vcf \


bcftools mpileup \
-Ov \
-f "$ref" \
results/map/ERR246969.sorted.bam \
| bcftools call -mv -Ov -o results/variants/9.vcf \

bcftools mpileup \
-Ov \
-f "$ref" \
results/map/ERR246968.sorted.bam \
| bcftools call -mv -Ov -o results/variants/ERR246968.vcf \


# As the final step, we will simpl filter our variants. If you open any of the vcf files, you wil find a lot of information, 
# among them are columns: #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO. In the last colum, there are details such as 
# DP=1;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,1,0;MQ=60
# All these store information that we can use to filter the variants. In our case, we want to filter and keep only variants with
# A base call quality of > 20 (if you read this paper https://www.nature.com/articles/s41525-021-00227-3, you will realize that QUAL > 20)
# will give you high quality variants. We will also filter where DP is > 10

bcftools filter \
-Ov \
-o results/variants/ERR246968.filtered.vcf \
-i "QUAL > 20 && DP > 10" \
results/variants/ERR246968.vcf \

echo "variant calling completed..."