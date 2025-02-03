# We have just gone through how we can write a workflow joining five steps using bash
# The probblem with that approach is that we fail to attain adaptability and transperency (see Molder et al.)
# That is why we will turn to snakemake, which after we will see, it makes out work easier and allows us to be more robust

# Caution: This tutorial is in no way the entirity of snakemake. I just want to pique your interest on the workflow
# This means that you will (just like I do to date) constantly consult the snakemake documentation for more advanced features
# See documentation at https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html

# Some housekeeping
# - File extensions are generally important in data analtysis and bioinformatics. In snakemake, we can use both .smk (we will use this constantly)
# Or .snakemake. We can also use snakefile but we will see why the other extensions are better because we will reverse "snakefile" for a more
# important task
# I will introduce the concept of a config file as early as now and we will see how important it is.
# As you learn, type out the code in your own files and try and run.

# Let's gooooo...
# Start by creating a file, call it whatever name you want and end it with a .smk. In my case this will be
# main.smk
# Let us rewind a bit. In our bash workflow, we ran five steps that followed each other. And we needed to put them in the exact same
# order so that one important step does not run before another prelliminary step (the pipeline would scream for errors).
# In snakemake we will not care about this order. However, instead of refering to these are steps, we call them rules.
# So the step for downloading the data is simply a rule. We can call it rule download_files. Same to running fastqc and fastp. We can call these
# rule fastqc and rule fastp respectively. So rules in snakemake are the steps or actions we run. 
# In order to create a rule, we simply involve the name "rule" which is a snakemake functionality then give the rule a name.
# In our case, the first rule would be to download our raw files, we we can simply do this:

rule download_files:
# under this rule, we can the describe whatever we want the rule to do for us. In our case, we simply want it to download for us 
# fastq reads. An important, note, when working with rules, you must know beforehand what you want the rule to do for you, how 
# you want the rule to perform the function, and what output you want after the rule has been run.
# So the basic structure of a rule is always something like this:
rule download_files:
    output:
    shell:
        """
        """
# Here we have just defined the rule, stated that we expectd some output, and told snakemake that we will be running a shell script
# Let us populate the rule to see what we now have:

rule download_files:
    output:
        "data/ERR246969_1.fastq",
        "data/ERR246969_2.fastq",
    shell:
        """
        data_dir="data"
        fasterq-dump -S ERR246969 -O "$data_dir"
        """
# Let us go through what just happened. We created a rule then described the outputs we expected from that rule and how we want these outputs to be created
# note that I have specified the outputs in "" (double quotes, you can use single quotes) and separated them using a coma.
# If you have multiple inputs, you will always separate them using a comma. The shell just contains the commands we want to run. You
# Can see that we can even define our variables there.
# Running snakemake is as easy as 

# BUT remember, here we just downloaded one file. We want to work with many files. For this reason, we will shift to whatever we did with the advanced.sh
# Before we go further, I want to introduce the concept of expand(). When we want to return many files, we can either list all files like we did in the first exame
# or use the expand() function and simply supply the natures of files. expand() will make our code cleaner.
# How expand() works. expand("the/file/{any_pattern}.file_ext", any_pattern=["this pattern", "another_patter"])
# Let us look at an example.
# expand("data/{IDs}_{read}.fastq", read=["_1", "_2"], IDs=["ERR246968", "ERR246969"]). Expand will ideally return "data/ERR246969_1.fastq", 
# "data/ERR246969_2.fastq",  "data/ERR246968_1.fastq", "data/ERR246968_2.fastq". If it is not clear yet, I hope it will be as we continue to use it.
# Let us try and invoke it in this command. The outputs I expect are fastq reads (both one and two) for these SRRs:
# ERR246968, ERR246969, ERR246970, ERR246971, ERR246972, ERR246973, ERR246974, ERR246975. We will use expand to simplify this
# Note, I will just paste the command we had previously because we had already used it at it worked.

rule download_files:
    output:
        expand("data/{id}_{read}.fastq", id=["ERR246968",
        "ERR246969",
        "ERR246970",
        "ERR246971",
        "ERR246972",
        "ERR246973",
        "ERR246974",
        "ERR246975"], read=["_1", "_2"])
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
# Before you ran the script, there are important things that you may miss. When accessing variables within shell
# in snakemake, as in this line for fq in "${{accessions[@]}}", we use double curly brackets like {{}}. I wll explain
# why this is the case later. Meanwhile, just trust that it works at the moment
# So far, we have looked at what a rule is, its structure, and how we can use expand to work with many output files

# second rule: fastqc
# In some cases, you want the rule to start from some input files. In our case, we want the fastqc rule to start from the raw reads,
# process these files using certain commands, then give us some input. So we understand our inputs, know what we want to output, and how 
# to generate these outputs. In this case, we add an additional functionality/component of the rule, input. So our updated rule would be 
# something like this:
rule fastqc:
    input:
    output:
    shell:
        """
        """
# With this structure, we can then populate the rule using the command we previously used in the advanced.sh
# We will also use expand to handle both our inputs and outputs. Ideally, fastqc outputs a html file and zipped folder for each file it handles
# When using expand, we can actually define some of extension outside the rule, especially if these are things we will constantly reuse.

ids=["ERR246968",
        "ERR246969",
        "ERR246970",
        "ERR246971",
        "ERR246972",
        "ERR246973",
        "ERR246974",
        "ERR246975"]
read=["1", "2"]
rule fastqc:
    input:
        expand("data/{id}_{read}.fastq", id=ids, read=read)
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
# Before you run this script, let me introduce you to rule all. This is a special rule that you will always include at the top of your snakefile
# if you have more than one rule that you expect to depend on each other. If you have more than two rules but fail to include rule all at the top,
# snakemake will only process the rule that is on top. In defining rule all, theoritically you call include all components of rule as we have done 
# previously. However, in this case, let us only focus on the input for now. We will build on this as we progress.
# The difference with rule all is that the input here is the output you expect from the rule you expect to run last. In our case, the input
# for rule all would be the output we expect from running fastqc. Snakemake will, therefore, look at what we specify in rule all, then 
# figure our which rule should be creating these files. It will check if these files are missing and if missing, execute the exact rule that 
# should be creating these files. The fun part is that if this particular rule, say fastqc depends on another rule say, rule download_files,
# snakemake will check to see if the files rule fastqc needs are there and if not there, execute the rule that should be generating those files. This creates a
# kind of cyclic DAG that is stepwise. Isn't that amazing?
# We will continue to uncover the intricacies of snakemake as we move a long. For now, just define rule all at the top of your snakefile and include
# the output we expect from rule fastqc
rule all:
    input:
        expand("results/fastqc/{id}_{read}_fastqc.{ext}", id=ids, ext=["html", "zip"], read=read)
# Rule 3: Fastp

# We will use the same concepts we learned from earlier practices, to create rule fastp. We will just adopt the command we created earlier to do this.
# Remember, that fastp gives out a json file, a html file, and trimmed reads (1, 2). So we will specify all these in our rule:
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
# Take notes about the double curly brackets {{}} that we are adding when accessing bash variables within the shell. We use {{}} because single curly bracked 
# {} are reserved for snakemake functionalities. We will see this later.

# rule 4: Downloading and indexing reference genome

# Since we know how to download and index a reference genome using bash, we will just get a snakemake version of it
# We must appreciate that indexing with bwa results in files with the folllowing extensions: "amb", "ann", "bwt", "pac", "sa". So 
# we need to specify them in our output
ref_name="GCF_000001405.26_GRCh38_genomic.fna.gz"

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
        bwa index "$ref_path"
        """
# A few notes before we run this script. In this case, we have defined an external variable ref_name="GCF_000001405.26_GRCh38_genomic.fna.gz" that we want 
# to access in the snakemake rule. This time, we will invoke the f string to help us achieve this. 
# Also, take note that are including {{}} in the expand. This is bceuse both f strig and expand() use {}, so the second {} is to escpae f string's
# utility.
# We have also introduced another component of snakemake rules called params. This is not the best place to use it as we will see in other examples.
# However, we are using params to define the reference variable we had set earlier in the global scope. This will allow us to access this variable
# seamlessly in the shell. Important, when accessing snakemake rule functionality and variables in the shell, we use the dot notation.
# This means, that we can access the variable ref_name under params using params.ref_name. We will exessively make use of this idea in subsequent rules.

# rule 4: mapping
# As we have done previously, we will simply adopt the commands for mapping we had in advanced.sh. This time, we will do some slight modifications. We
# expect to work with the ref file that we downloaded in the rule download_index_ref and use the files from the rule fastp. 
# So, we will just copy the outputs of these rules and put them as input in our rule map. This is one way of connecting the rules. So, here, we are making 
# rule map depend on rule download_index_ref. In the output we expected to have sorted bams and files with bam.bai.
# In the shell script, I have just combined all the possible commands we need to get from raw fastq files to sorted bam files. We had done this in advanced.sh
# Also take note of how we are defining and accessing variables for input and params functionalities. We are using input.var and params.var that we learned earlier.
# This also tells us that we can define as many variables as possible then access them in our shell script as much as we want.
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
# rule variant call
 # We will finish our workflow by defining the rule variant call. For this rule, we will continue to see how we can access variables set 
 # within inputs or params in the shell script. We will also continue to learn how we can make rules depend on each other by taking the output
 # from one rule and using it as input for the next rule. This basically give snakemake some hints on the order of how things should run in our pipeline
 # In the rule variant_call, we will use the bam outputs from the rule map as inputs. We will also be inputing our index reference genome
 # as an input file here. The output we expect are filtered variant call files.
 # For this rule, I would wish us to split it into two: one that does the variant call and the other that does the variant filtering.
 # This is just to tell the program to call the variants then filter just in case there were some bam files in the directory we will be working prior to our 
 # processing.


rule variant_call:
    input:
        bams=expand("results/map/{id}.sorted.{ext}", id=ids, ext=["bam", "bam.bai"]),
        ref=f"ref/{ref_name}",
        index_files=expand(f"ref/{ref_name}.{{ext}}", ext=["amb", "ann", "bwt", "pac", "sa"])
    output:
        expand("results/variants_1/{id}.vcf", id=ids)
    params:
        vcf_call: "-Ov -f",
        bcf_tool_call: "-mv  -Ov"
    shell:
        """
        map_dir="results/map"
        variants="results/variants_1"
        for so_bam in "${{map_dir}}"/*.sorted.bam; do
        vcf="$variants_1/$(basename "$so_bam" .sorted.bam).vcf"
        bcftools mpileup {params.vcf_call} {input.ref} "${so_bam}" | bcftools call {params.bcf_tool_call} -o "${vcf}"
        done
        """
# Just a slight caution, when I was using the reference genome file, I got an error [E::fai_build3_core] Cannot index files compressed with gzip, please use bgzip
# To address this error, I just unzipped the file then bgzipped it. I have since added a bashcript that checks the zipping format of the file then if
# not bgzipped, it just unzips it the bgzips it back.

# Let us then finish the pipeline by writing the variant filtering rule:

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
# We maintan everything in this rule. However, I want you to take note of what is happening with the params this time.
# The "{params.vcf_filter}" is now enclosed in double quotes, so Bash correctly passes it as a single argument to bcftools.
# This prevents Bash from misinterpreting && as a shell logical operator.

# Now we have finally finished our mini workshop. By no way have I given you all that is there to be know about snakemake. But I hope I have
# piqued your interest to continue learning more using other materials. The snakemake documentation provides a good support so I would 
# recommend that you check that out. Continue practicing on your own. I have included a test for you in the test.smk file to help you develop
# a new workflow and publish in your github. I will be doing more of this tutorials in the future and if you are interested, hang with me in my 
# socials. Also, if you have recommendations on what you want to learn or questions about what we have covered, you can write to me directlt through
# my socials and email. Good luck.