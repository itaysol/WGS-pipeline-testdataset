import os
import sys
from snakemake.io import glob_wildcards
import yaml
import csv
import itertools
from parser import parser,setComparisonGroups


configfile : parser(config["file"])
comparisonGroupTuples = setComparisonGroups(config)


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
        expand("output/assembly/comparisonGroup{comparisonGroup}/{sample}_assembly", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand("output/wgkb/comparisonGroup{comparisonGroup}/{sample}/{sample}.kraken_taxonomy.txt", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand("output/wgkb/comparisonGroup{comparisonGroup}/{sample}/{sample}.kraken_output.txt", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand("output/wgkb/comparisonGroup{comparisonGroup}/{sample}/{sample}.bracken_output.txt", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand("output/wgv/{sample}/{sample}.validation_output.txt", sample = config["Samples"].keys()),
        expand("output/mlst/comparisonGroup{comparisonGroup}/{sample}/{sample}.contigs.mlst.tsv", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand("output/rmlst/comparisonGroup{comparisonGroup}/{sample}.tsv", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand("output/abricate/comparisonGroup{comparisonGroup}/{sample}/{sample}.abricate.card.tsv", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand("output/AMR/comparisonGroup{comparisonGroup}/{sample}.amrfinderplus.tsv", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        


rule process_file_pair:
    conda:
        "env/conda-qc.yaml"
    threads:
        3
    input:
        fwd = lambda wildcards: config["Samples"][wildcards.id[:7]]["R1Fastq"],
        rev = lambda wildcards: config["Samples"][wildcards.id[:7]]["R2Fastq"]
    output:
        r1_clean = temporary(os.path.join(output_dir, "fastq/{id}/{id}.clean_1.fastq.gz")),
        r2_clean = temporary(os.path.join(output_dir, "fastq/{id}/{id}.clean_2.fastq.gz")),
        html = temporary(os.path.join(output_dir, "fastq/{id}/{id}.fastp.html")),
        json = temporary(os.path.join(output_dir, "fastq/{id}/{id}.fastp.json"))
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
        kraken2 --db /workspace/WGS-pipeline/databases/k2 --threads 4 --report {output.kraken_report} --output {output.kraken_output} \
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
            bracken -i {input.kraken_report_file} -d /workspace/WGS-pipeline/databases/k2 -o {output.bracken_output}
            
        """

        
rule sample_validation:
    input: 
        bracken_output_file = lambda wildcards: os.path.join(output_dir, "bracken", f"{wildcards.id}.bracken_output.tsv"),
    output:
        sample_validation_output = os.path.join(output_dir, "sample_validation", "{id}.output.txt")
    params:
        specie = lambda wildcards: config["Samples"][wildcards.id[:7]]["specie"],
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
    params:
        comparisonGroup = lambda wildcards: config["Samples"][wildcards.id[:7]]["comparisonGroup"]
    output:
        assembly_output = directory("output/assembly/{comparisonGroup}/{id}_assembly")
    shell:
        """
        echo Comparison Group: {params.comparisonGroup}
        shovill --R1 {input.clean_fwd} --R2 {input.clean_rev} --outdir {output.assembly_output} --assembler skesa
        
        """

rule whole_genome_krak_brack:
    conda:
        "env/conda-whole_genome_validation.yaml"
    params:
       # comparisonGroup = lambda wildcards: config["Samples"][wildcards.id[:7]]["comparisonGroup"],
        specie = lambda wildcards: config["Samples"][wildcards.id[:7]]["specie"],
        sample_id = lambda wildcards: wildcards.id
    input:
        contigs_file = os.path.join(output_dir, "assembly","{comparisonGroup}","{id}_assembly")
    output:
        wgv_kraken_report = os.path.join(output_dir, "wgkb","{comparisonGroup}" ,"{id}", "{id}.kraken_taxonomy.txt"),
        wgv_kraken_output = os.path.join(output_dir, "wgkb","{comparisonGroup}", "{id}", "{id}.kraken_output.txt"),
        wgv_bracken_output = os.path.join(output_dir, "wgkb", "{comparisonGroup}","{id}", "{id}.bracken_output.txt"),
    shell:
        """
        kraken2 --db /workspace/WGS-pipeline/databases/k2 \
         --threads 4 --report {output.wgv_kraken_report} --output {output.wgv_kraken_output} {input.contigs_file}/contigs.fa

        bracken -i {output.wgv_kraken_report} -d /workspace/WGS-pipeline/databases/k2\
         -o {output.wgv_bracken_output}

        """
rule whole_genome_validation:
    input:
        bracken_output =  os.path.join(output_dir, "wgkb", "{id}", "{id}.bracken_output.txt"),
    params:
        specie = lambda wildcards: config["Samples"][wildcards.id[:7]]["specie"],
        sample_id = lambda wildcards: wildcards.id
    output:
        wgv_sample_validation_output = os.path.join(output_dir, "wgv", "{id}.validation_output.txt")
    shell:
        """
            python samp_val.py --input {input.bracken_output} --output {output.wgv_sample_validation_output} {params.specie} {params.sample_id}  
    
        """
rule run_mlst:
    conda:
        "env/conda-mlst.yaml"
    params:
        comparisonGroup = lambda wildcards: config["Samples"][wildcards.id[:7]]["comparisonGroup"]
    input:
        contigs_file = os.path.join(output_dir, "assembly","{comparisonGroup}", "{id}_assembly")
    output:
        mlst_output = os.path.join(output_dir, "mlst","{comparisonGroup}", "{id}", "{id}.contigs.mlst.tsv"),
        rmlst_output = os.path.join(output_dir, "rmlst", "{comparisonGroup}","{id}.tsv")
    shell:
        """
        mlst {input.contigs_file}/contigs.fa > {output.mlst_output}
        bash ./rMLST/rmlst.sh {input.contigs_file}/contigs.fa > {output.rmlst_output}
        
        """
rule resistome_and_virulome:
    conda:
        "env/conda-resistome_and_virulome.yaml"
    params:
        comparisonGroup = lambda wildcards: config["Samples"][wildcards.id[:7]]["comparisonGroup"]
    input:
        contigs_file = os.path.join(output_dir, "assembly","{comparisonGroup}","{id}_assembly"),
    output:
        abricate_card = os.path.join(output_dir,"abricate","{comparisonGroup}","{id}","{id}.abricate.card.tsv"),
        abricate_vfdb = os.path.join(output_dir,"abricate","{comparisonGroup}","{id}.abricate.vfdb.tsv"),
        AMR = os.path.join(output_dir,"AMR","{comparisonGroup}","{id}.amrfinderplus.tsv")
    shell:
        """
        abricate --db card --minid 60 --mincov 60 {input.contigs_file}/contigs.fa > {output.abricate_card}
        abricate --db vfdb --minid 60 --mincov 60 {input.contigs_file}/contigs.fa > {output.abricate_vfdb}
        amrfinder -u
        amrfinder --nucleotide {input.contigs_file}/contigs.fa --plus --threads 4 -o {output.AMR}
        
        """



    

            
    
    


        
