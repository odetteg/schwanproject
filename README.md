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