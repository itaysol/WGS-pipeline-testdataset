import os
import sys
from snakemake.io import glob_wildcards
import yaml
import csv
import itertools
import pandas as pd
from scripts.parser import *
from datetime import datetime




configfile : parser(config["file"])
comparisonGroupTuples, cgSpecieDict, cgCounter,specieGenomeSizeDict, cgSampleDict = setComparisonGroups(config)
output_dir = "output"
krakenDB = "/workspace/WGS-pipeline/databases/k2"




rule process_file_pair:
    wildcard_constraints:
        id = r'[^/]+'
    conda:
        "env/conda-qc.yaml"
    threads:
        3
    priority:
        50
    input:
        fwd = lambda wildcards: config["Samples"][wildcards.id]["R1Fastq"],
        rev = lambda wildcards: config["Samples"][wildcards.id]["R2Fastq"]
    output:
        r1_clean = os.path.join(output_dir, "fastq","{id}","{id}.clean_fwd.fastq.gz"),
        r2_clean = os.path.join(output_dir, "fastq","{id}","{id}.clean_rev.fastq.gz"),
        html = os.path.join(output_dir, "fastq","{id}","{id}.fastp.html"),
        json = os.path.join(output_dir, "fastq","{id}","{id}.fastp.json")
    log:
        os.path.join(output_dir, "fastq","{id}","{id}.fastp.log.txt")
    shell:
        """
        fastp -i {input.fwd} -I {input.rev} \
              --out1 {output.r1_clean} --out2 {output.r2_clean} \
              -w 3 -h {output.html} -j {output.json}
        """

rule run_kraken2:
    conda:
        "env/conda-kraken_and_bracken.yaml"
    params:
        db = krakenDB
    input:
        clean_fwd=lambda wildcards: os.path.join(output_dir, "fastq", wildcards.id, f"{wildcards.id}.clean_fwd.fastq.gz"),
        clean_rev=lambda wildcards: os.path.join(output_dir, "fastq", wildcards.id, f"{wildcards.id}.clean_rev.fastq.gz")
    priority:
        25
    output:
        kraken_report = os.path.join(output_dir,"readSpeciesID","kraken", "{id}", "{id}.kraken_taxonomy.txt"),
        kraken_output = os.path.join(output_dir,"readSpeciesID","kraken", "{id}", "{id}.kraken_output.txt"),
    threads: 
         4
    shell:
        """
        kraken2 --db {params.db} --threads 4 --report {output.kraken_report} --output {output.kraken_output}\
        --paired {input.clean_fwd} {input.clean_rev}
        """

rule run_bracken:
    conda:
        "env/conda-kraken_and_bracken.yaml"
    params:
        db = krakenDB
    input:
       kraken_report_file= lambda wildcards: os.path.join(output_dir,"readSpeciesID", "kraken", wildcards.id, f"{wildcards.id}.kraken_taxonomy.txt")
    output:
        bracken_output= os.path.join(output_dir,"readSpeciesID", "bracken", "{id}.bracken_output.tsv")
    shell:
        """
            bracken -i {input.kraken_report_file} -d {params.db} -o {output.bracken_output}
            
        """

        
rule sample_validation:
    input: 
        bracken_output_file = lambda wildcards: os.path.join(output_dir,"readSpeciesID", "bracken", f"{wildcards.id}.bracken_output.tsv"),
    output:
        sample_validation_output = os.path.join(output_dir, "readSpeciesID", "sample_validation", "{id}.output.txt")
    params:
        specie = lambda wildcards: config["Samples"][wildcards.id[:7]]["specie"],
        sample_id = lambda wildcards: wildcards.id
    shell:    
         """
         python scripts/samp_val.py --input {input.bracken_output_file} --output {output.sample_validation_output} {params.specie} {params.sample_id}
         
         """
