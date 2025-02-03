# The schwannomatosis project
 The purpose of this repository of to analyze in part, whole genome data from the Piotrowski et al. study.

## Resources required
- Personal computer: Can be mac or windows. If you cannot do analysis in your personal computer, we will look
  at alternatives ways in the cloud/online
- Human reference genome
- Code editor. I prefer and will use visual studio. You can use sublime text and any other code editor of your choice.
- 
## Bioinformatics tools used
- SRA toolkit
- FASTQC
- FASTP
- Trimmomatic
- MULTIQC
- Burrows-Wheeler aligner
- SAMtools
- Platipus variant caller
- SeattleSeq
- Snakemake

## Analysis pipilleline
- Download raw fastq files using fastqdump
- Perform quality checks using FASTQC and FASTP
- Peroform filtering using FASTP and Trimmomatic
- Map the reads to the reference genome

## Installing SRA toolkit
- Depending on your operating system, go to the link https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit and folllow the installation instructions

## How to use this module
The tutorials are all written in the files prim.sh, advanced.sh, and prim.smk. You may want to read those as you code along. You can also read my medium blog as you code along. 
Ensure you install all tools needed before you run the entire thing. Visit my tool_installation.sh to see how to install some of the tools.
The files main.sh and main.smk contain end to end pipelines. You can refer to these as you code along and when you encounter any error. In case of any error, you can always use chatgpt (be cautious, sometimes its answers do not work) and most certainly stackexchange.
Depending on your space, you can choose to work with as few raw files as possible. I include vcf files for all 8 samples just in case you need to work with them.

## Additional materials/resorces
- [Snakemake documentation ](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html)