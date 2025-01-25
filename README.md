# The purpose of this repository is to analyze the data in the study by Piotrowski et al.

## Resources required
- Personal computer: Can be mac or windows. If you cannot do analysis in your personal computer, we will look
  at alternatives ways in the cloud/online
- Human reference genome
- Code editor. I prefer and will use visual studio. You can use sublime text and any other code editor of your choice.
- 
## Bioinformatics tools used
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
- Download raw fastq files using wget
- Perform quality checks using FASTQC and FASTP
- Peroform filtering using FASTP and Trimmomatic
- Map the reads to the reference genome