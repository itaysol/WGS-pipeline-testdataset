import os
import sys
from snakemake.io import glob_wildcards
import yaml
import csv
import itertools
from scripts.parser import *
from datetime import datetime



configfile : parser(config["file"])
comparisonGroupTuples = setComparisonGroups(config)

# Define the output directory
#timeStamp = datetime.now().strftime("%Y.%m.%d.%H.%M.%S")
#in order to resume a failed run, we need to use the same output directory name - so I added it manually
#timeStamp = "20231117000901"
output_dir = "output"

# Define the final rule that specifies the targets to generate
rule all:
    input:
        expand(output_dir+"/fastq/{sample}/{sample}.clean_fwd.fastq.gz", sample = config["Samples"].keys()),
        expand(output_dir+"/fastq/{sample}/{sample}.clean_rev.fastq.gz", sample = config["Samples"].keys()),
        expand(output_dir+"/kraken/{sample}/{sample}.kraken_taxonomy.txt", sample = config["Samples"].keys()),
        expand(output_dir+"/kraken/{sample}/{sample}.kraken_output.txt", sample = config["Samples"].keys()),
        expand(output_dir+"/bracken/{sample}.bracken_output.tsv", sample = config["Samples"].keys()),
        expand(output_dir+"/sample_validation/{sample}.output.txt", sample = config["Samples"].keys()),
        expand(output_dir+"/assembly/comparisonGroup{comparisonGroup}/{sample}_assembly", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        #expand(output_dir+"/assemblyContigsOnly/comparisonGroup{comparisonGroup}", comparisonGroup = [item[0] for item in comparisonGroupTuples]),
        expand(output_dir+"/wgkb/comparisonGroup{comparisonGroup}/{sample}/{sample}.kraken_taxonomy.txt", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/wgkb/comparisonGroup{comparisonGroup}/{sample}/{sample}.kraken_output.txt", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/wgkb/comparisonGroup{comparisonGroup}/{sample}/{sample}.bracken_output.txt", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/wgv/{sample}/{sample}.validation_output.txt", sample = config["Samples"].keys()),
        expand(output_dir+"/mlst/comparisonGroup{comparisonGroup}/{sample}/{sample}.contigs.mlst.tsv", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/rmlst/comparisonGroup{comparisonGroup}/{sample}.tsv", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/abricate/comparisonGroup{comparisonGroup}/{sample}/{sample}.abricate.card.tsv", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/AMR/comparisonGroup{comparisonGroup}/{sample}.amrfinderplus.tsv", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/cgMLST/comparisonGroup{comparisonGroup}/{sample}/schema", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/cgMLST/comparisonGroup{comparisonGroup}/{sample}/Allelecall", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/cgMLST/comparisonGroup{comparisonGroup}/{sample}/MST", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/cgMLST/comparisonGroup{comparisonGroup}/{sample}/grapetree/cgMLST100.tre", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),       
        
        
rule process_file_pair:
    conda:
        "env/conda-qc.yaml"
    threads:
        3
    input:
        fwd = lambda wildcards: config["Samples"][wildcards.id[:7]]["R1Fastq"],
        rev = lambda wildcards: config["Samples"][wildcards.id[:7]]["R2Fastq"]
    output:
        r1_clean = os.path.join(output_dir, "fastq/{id}/{id}.clean_fwd.fastq.gz"),
        r2_clean = os.path.join(output_dir, "fastq/{id}/{id}.clean_rev.fastq.gz"),
        html = temporary(os.path.join(output_dir, "fastq/{id}/{id}.fastp.html")),
        json = temporary(os.path.join(output_dir, "fastq/{id}/{id}.fastp.json"))
    log:
        output_dir+"/fastq/{id}/{id}.fastp.log.txt"
    shell:
        """
        fastp -i {input.fwd} -I {input.rev} \
              --out1 {output.r1_clean} --out2 {output.r2_clean} \
              -w 3 -h {output.html} -j {output.json}
        """

rule run_kraken2:
    conda:
        "env/conda-kraken_and_bracken.yaml"
    input:
        clean_fwd=lambda wildcards: os.path.join(output_dir, "fastq", wildcards.id, f"{wildcards.id}.clean_fwd.fastq.gz"),
        clean_rev=lambda wildcards: os.path.join(output_dir, "fastq", wildcards.id, f"{wildcards.id}.clean_rev.fastq.gz")
    output:
        kraken_report = os.path.join(output_dir, "kraken", "{id}", "{id}.kraken_taxonomy.txt"),
        kraken_output = os.path.join(output_dir, "kraken", "{id}", "{id}.kraken_output.txt"),
    threads: 4
    shell:
        """
        kraken2 --db /workspace/WGS-pipeline/databases/k2 --threads 4 --report {output.kraken_report} --output {output.kraken_output}\
        --paired {input.clean_fwd} {input.clean_rev}
        """

rule run_bracken:
    conda:
        "env/conda-kraken_and_bracken.yaml"
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
         python scripts/samp_val.py --input {input.bracken_output_file} --output {output.sample_validation_output} {params.specie} {params.sample_id}
         
         """

rule assembly:
    conda:
         "env/conda-assembly.yaml"
    priority:
        50
    input:
        clean_fwd = os.path.join(output_dir, "fastq", "{id}", "{id}.clean_fwd.fastq.gz"),
        clean_rev = os.path.join(output_dir, "fastq", "{id}", "{id}.clean_rev.fastq.gz"),
    params:
        comparisonGroup = lambda wildcards: config["Samples"][wildcards.id[:7]]["comparisonGroup"]
    output:
        assembly_output = directory(output_dir+"/assembly/{comparisonGroup}/{id}_assembly")
    shell:
        """
        shovill --R1 {input.clean_fwd} --R2 {input.clean_rev} --outdir {output.assembly_output} --assembler skesa
        
        """

rule whole_genome_krak_brack:
    conda:
        "env/conda-kraken_and_bracken.yaml"
    params:
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
        python scripts/samp_val.py --input {input.bracken_output} --output {output.wgv_sample_validation_output} {params.specie} {params.sample_id}  
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
rule create_training_file:
    conda:
        "env/conda-prodigal.yaml"
    input:
        contigs_file = os.path.join(output_dir, "assembly","{comparisonGroup}","{id}_assembly")
    output:
       training_file = os.path.join(output_dir,"trainingFiles","{comparisonGroup}","{id}.trn")
    shell:
        """
        prodigal -i {input.contigs_file}/contigs.fa -t {output.training_file} -p single

        """    


rule adhoc_cgMLST:
    conda:
        "env/conda-chewBBACA.yaml"
    params:
        comparisonGroup = lambda wildcards: config["Samples"][wildcards.id[:7]]["comparisonGroup"]
    input:
        training_file = os.path.join(output_dir,"trainingFiles","{comparisonGroup}","{id}.trn")
    output:
        createSchema = directory(os.path.join(output_dir,"cgMLST","{comparisonGroup}","{id}","schema")),
        allelecall = directory(os.path.join(output_dir,"cgMLST","{comparisonGroup}","{id}","Allelecall")),
        MST = directory(os.path.join(output_dir,"cgMLST","{comparisonGroup}","{id}","MST"))

    shell:
       """
        for comp_group_dir in output/assembly/comparisonGroup*/; do
            comp_group=$(basename "$comp_group_dir")
            output_dir="output/assemblyContigsOnly/$comp_group"

            # Create the output directory if it doesn't exist
            mkdir -p "$output_dir"

            # Copy contigs.fa from each sample assembly to the output directory
            for assembly_dir in "$comp_group_dir"/*/; do
                sample_id=$(basename "$assembly_dir")
                cp "$assembly_dir/contigs.fa" "$output_dir/${{sample_id}}_contigs.fa"
            done
        done

        chewBBACA.py CreateSchema -i output/assemblyContigsOnly/comparisonGroup{params.comparisonGroup}/ -o {output.createSchema} --ptf {input.training_file} --cpu 20

        # Loop until the condition is met
        until chewBBACA.py AlleleCall -i output/assemblyContigsOnly/comparisonGroup{params.comparisonGroup} -g {output.createSchema}/schema_seed/ -o {output.allelecall} --cpu 14; do
            if [ ! -f {output.allelecall}/results_statistics.tsv ]; then
                continue
            fi

            if python -c "import sys, os; sys.exit(0 if os.path.exists('{output.allelecall}/results_statistics.tsv') and all(line.split('\\t')[1].strip() == '0' for line in open('{output.allelecall}/results_statistics.tsv')) else 1)"; then
                break
            fi
        done

        cut -f 2 {output.allelecall}/paralogous_loci.tsv | sed 's/|/\\n/g' {output.allelecall}/paralogous_loci.tsv
        chewBBACA.py RemoveGenes -i {output.allelecall}/results_alleles.tsv -g {output.allelecall}/paralogous_loci.tsv -o {output.allelecall}/results_alleles.no_paralogs.tsv
        chewBBACA.py ExtractCgMLST -i {output.allelecall}/results_alleles.no_paralogs.tsv -o {output.MST}
        """

rule grapetree:
    conda:
        "env/conda-grapetree.yaml"
    input:
        grapetree_input = os.path.join(output_dir,"cgMLST","{comparisonGroup}","{id}","MST")
    output:
        grapetree_output = os.path.join(output_dir,"cgMLST","{comparisonGroup}","{id}","grapetree","cgMLST100.tre")
    shell:
        """
        grapetree --profile {input.grapetree_input}/cgMLST100.tsv > {output.grapetree_output}

        """

    