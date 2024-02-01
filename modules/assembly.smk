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

rule assembly:
    wildcard_constraints:
        id = r'[^/]+'
    conda:
         "env/conda-assembly.yaml"
    priority:
        50
    threads: 
        8
    input:
        clean_fwd = os.path.join(output_dir, "fastq", "{id}", "{id}.clean_fwd.fastq.gz"),
        clean_rev = os.path.join(output_dir, "fastq", "{id}", "{id}.clean_rev.fastq.gz"),
    params:
        genome_size = lambda wildcards: specieGenomeSizeDict[config["Samples"][wildcards.id[:7]]["specie"]]
    output:
        assembly_dir = directory(os.path.join(output_dir, "assembly","{id}_assembly")),
        assembly_output = os.path.join(output_dir, "assembly","{id}_assembly","contigs.fa")
    shell:
        """
        shovill --cpus {threads} --R1 {input.clean_fwd} --R2 {input.clean_rev} --force --gsize {params.genome_size}M --outdir {output.assembly_dir} --assembler skesa
        
        """

rule whole_genome_krak_brack:
    wildcard_constraints:
        id = r'[^/]+'
    conda:
        "env/conda-kraken_and_bracken.yaml"
    params:
        specie = lambda wildcards: config["Samples"][wildcards.id[:7]]["specie"],
        sample_id = lambda wildcards: wildcards.id,
        db = krakenDB
    input:
        contigs_file = os.path.join(output_dir, "assembly","{id}_assembly","contigs.fa")
    output:
        wgv_kraken_report = os.path.join(output_dir, "assemblySpeciesID","wgkb","{id}", "{id}.kraken_taxonomy.txt"),
        wgv_kraken_output = os.path.join(output_dir, "assemblySpeciesID","wgkb", "{id}", "{id}.kraken_output.txt"),
        wgv_bracken_output = os.path.join(output_dir, "assemblySpeciesID","wgkb","{id}", "{id}.bracken_output.txt"),
    shell:
        """
        kraken2 --db {params.db} \
         --threads 4 --report {output.wgv_kraken_report} --output {output.wgv_kraken_output} {input.contigs_file}

        bracken -i {output.wgv_kraken_report} -d {params.db}\
         -o {output.wgv_bracken_output}

        """
rule whole_genome_validation:
    wildcard_constraints:
        id = r'[^/]+'
    input:
        bracken_output =  os.path.join(output_dir, "assemblySpeciesID","wgkb","{id}", "{id}.bracken_output.txt"),
    params:
        specie = lambda wildcards: config["Samples"][wildcards.id[:7]]["specie"],
        sample_id = lambda wildcards: wildcards.id
    output:
        wgv_sample_validation_output = os.path.join(output_dir ,"assemblySpeciesID" ,"wgv" ,"{id}.validation_output.txt")
    shell:
        """
        python scripts/samp_val.py --input {input.bracken_output} --output {output.wgv_sample_validation_output} {params.specie} {params.sample_id}  
        """
rule run_mlst:
    conda:
        "env/conda-mlst.yaml"
    params:
        comparisonGroup = lambda wildcards: config["Samples"][wildcards.id[:7]]["comparisonGroup"]
    input:
        contigs_file = os.path.join(output_dir, "assembly" , "{id}_assembly","contigs.fa")
    output:
        mlst_output = os.path.join(output_dir,"assemblySpeciesID" ,"mlst","{id}.contigs.mlst.tsv"),
        rmlst_output = os.path.join(output_dir,"assemblySpeciesID", "rmlst","{id}.tsv")
    shell:
        """
        mlst {input.contigs_file} > {output.mlst_output}
        python scripts/rmlst.py --file {input.contigs_file} > {output.rmlst_output}
        
        """
rule resistome_and_virulome:
    conda:
        "env/conda-resistome_and_virulome.yaml"
    params:
        comparisonGroup = lambda wildcards: config["Samples"][wildcards.id[:7]]["comparisonGroup"]
    input:
        contigs_file = os.path.join(output_dir, "assembly","{id}_assembly","contigs.fa"),
    output:
        abricate_card = os.path.join(output_dir,"abricate","card","{id}.abricate.card.tsv"),
        abricate_vfdb = os.path.join(output_dir,"abricate","vfdb","{id}.abricate.vfdb.tsv"),
        AMR = os.path.join(output_dir,"AMR","{id}.amrfinderplus.tsv")
    shell:
        """
        abricate --db card --minid 60 --mincov 60 {input.contigs_file} > {output.abricate_card}
        abricate --db vfdb --minid 60 --mincov 60 {input.contigs_file} > {output.abricate_vfdb}
        amrfinder -u
        amrfinder --nucleotide {input.contigs_file} --plus --threads 4 -o {output.AMR}
        """
