import os
import sys
from snakemake.io import glob_wildcards
import yaml
import csv
from parser import parser


configfile : parser(config["file"])

# Define the output directory
output_dir = "output"

# Define the final rule that specifies the targets to generate
rule all:
    input:
        expand("output/fastq/{sample}/{sample}.clean_1.fastq.gz", sample = config["Samples"].keys()),
        expand("output/fastq/{sample}/{sample}.clean_2.fastq.gz", sample = config["Samples"].keys()),
        expand("output/kraken/{sample}/{sample}.kraken_taxonomy.txt", sample = config["Samples"].keys()),
        expand("output/kraken/{sample}/{sample}.kraken_output.txt", sample = config["Samples"].keys()),
        expand("output/bracken/{sample}.bracken_output.tsv", sample = config["Samples"].keys()),
        expand("output/sample_validation/{sample}.output.txt", sample = config["Samples"].keys()),
        expand("output/assembly/{sample}/{sample}.fasta", sample = config["Samples"].keys())
        
rule process_file_pair:
    conda:
        "env/conda-qc.yaml"
    threads:
        3
    input:
        fwd = lambda wildcards: config["Samples"][wildcards.id[:7]]["R1Fastq"],
        rev = lambda wildcards: config["Samples"][wildcards.id[:7]]["R2Fastq"]
    output:
        r1_clean = os.path.join(output_dir, "fastq/{id}/{id}.clean_1.fastq.gz"),
        r2_clean = os.path.join(output_dir, "fastq/{id}/{id}.clean_2.fastq.gz"),
        html = os.path.join(output_dir, "fastq/{id}/{id}.fastp.html"),
        json = os.path.join(output_dir, "fastq/{id}/{id}.fastp.json")
    log:
        "output/fastq/{id}/{id}.fastp.log.txt"
    shell:
        """
        fastp -i {input.fwd} -I {input.rev} \
              -o {output.r1_clean} -O {output.r2_clean} \
              -w 3 -h {output.html} -j {output.json}
        """

rule run_kraken2:
    conda:
        "env/conda-kraken2.yaml"
    input:
        clean_fwd=lambda wildcards: os.path.join(output_dir, "fastq", wildcards.id, f"{wildcards.id}.clean_1.fastq.gz"),
        clean_rev=lambda wildcards: os.path.join(output_dir, "fastq", wildcards.id, f"{wildcards.id}.clean_2.fastq.gz")
    output:
        kraken_report = os.path.join(output_dir, "kraken", "{id}", "{id}.kraken_taxonomy.txt"),
        kraken_output = os.path.join(output_dir, "kraken", "{id}", "{id}.kraken_output.txt"),
    threads: 4
    shell:
        """
        kraken2 --db /workspace/Gene-pipeline/databases/k2/minikraken2_v2_8GB_201904_UPDATE --threads 4 --report {output.kraken_report} --output {output.kraken_output} \
            --paired {input.clean_fwd} {input.clean_rev}
        """

rule run_bracken:
    conda:
        "env/conda-bracken.yaml"
    input:
       kraken_report_file= lambda wildcards: os.path.join(output_dir, "kraken", wildcards.id, f"{wildcards.id}.kraken_taxonomy.txt")
    output:
        bracken_output= os.path.join(output_dir, "bracken", "{id}.bracken_output.tsv")
    shell:
        """
            bracken -i {input.kraken_report_file} -d /workspace/Gene-pipeline/databases/k2/minikraken2_v2_8GB_201904_UPDATE -o {output.bracken_output}
            
        """

        
rule sample_validation:
    input: 
        bracken_output_file = lambda wildcards: os.path.join(output_dir, "bracken", f"{wildcards.id}.bracken_output.tsv"),
    output:
        sample_validation_output = os.path.join(output_dir, "sample_validation", "{id}.output.txt")
    params:
        specie = lambda wildcards: config["Samples"][wildcards.id]["specie"],
        sample_id = lambda wildcards: wildcards.id
    shell:    
         """
         python samp_val.py --input {input.bracken_output_file} --output {output.sample_validation_output} {params.specie} {params.sample_id}
         
         """

rule assembly:
    conda:
         "env/conda-assembly.yaml"
    input:
        clean_fwd = os.path.join(output_dir, "fastq", "{id}", "{id}.clean_1.fastq.gz"),
        clean_rev = os.path.join(output_dir, "fastq", "{id}", "{id}.clean_1.fastq.gz"),
    output:
        assembly_output = directory("output/assembly/{id}.fasta")
    shell:
        """
        shovill --trim --R1 {input.clean_fwd} --R2 {input.clean_rev} --outdir {output.assembly_output} 
        
        """


        
