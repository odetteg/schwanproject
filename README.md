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
- Perform filtering using FASTP. You can also use Trimmomatic but this will not be covered at this time.
- Map the reads to the reference genome
- Call and filter variants

## Installing SRA toolkit
- Depending on your operating system, go to the link https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit and folllow the installation instructions

## How to use this module
<br>The tutorials are all written in the files prim.sh, advanced.sh, and prim.smk. You may want to read those as you code along. You can also read my medium blog as you code along.<br> 
<br>Ensure you install all tools needed before you run the entire pipeline or any code chunks. Visit my tool_installation.sh to see how to install some of the tools. This installation protocol ensures that you work within a dedicated conda environment.<br>
<br> The files main.sh and main.smk contain end to end pipelines. You can refer to these as you code along. In case of any error, you can always use chatgpt (be cautious, sometimes its answers do not work) and most certainly stackexchange. <br>
<br>Depending on your space, you can choose to work with as few raw files as possible. I include vcf files for all 8 samples just in case you need to work with them.<br>

## Additional materials/resorces
- [Snakemake documentation ](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html)

## Final note:
We have finally wrapped up our mini workshop on Snakemake! While this guide certainly does not cover everything there is to know about Snakemake, I truly hope it has sparked your interest in exploring more. The Snakemake documentation is an excellent resource, and I highly recommend checking it out to expand your knowledge further.
As with anything, practice is key - so keep experimenting and building workflows on your own. I have included a test in the test.smk file to give you a challenge and help you develop a new workflow to publish on your GitHub.
I will be doing more tutorials like this in the future, so if you're interested, feel free to connect with me on my socials. If you have any specific recommendations for what you would like to learn next, or if you have any questions about what we have covered, don't hesitate to reach out to me through my socials or email.
Good luck, and keep building!