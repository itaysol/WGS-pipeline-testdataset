import os
import sys
from snakemake.io import glob_wildcards

IDS, = glob_wildcards("fastq/{id}_1.fastq.gz")

# Define the output directory
output_dir = "output"

#Define the specie we want to check
specie="ppp"
if len(sys.argv) > 5:
    specie = sys.argv[5][7:]

print(specie)
# Define the final rule that specifies the targets to generate
rule all:
    input:
          expand("output/fastq/{id}/{id}.clean_1.fastq.gz", id=IDS), expand("output/fastq/{id}/{id}.clean_2.fastq.gz", id=IDS),
          expand("output/kraken/{id}/{id}.kraken_taxonomy.txt",id=IDS), expand("output/kraken/{id}/{id}.kraken_output.txt",id=IDS),
          expand("output/bracken/{id}.bracken_output.tsv",id=IDS)

rule process_file_pair:
    conda:
        "env/conda-qc.yaml"
    threads:
        3
    input:
        fwd = "fastq/{id}_1.fastq.gz",
        rev = "fastq/{id}_2.fastq.gz"
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
        clean_fwd="output/fastq/{id}/{id}.clean_1.fastq.gz",
        clean_rev="output/fastq/{id}/{id}.clean_2.fastq.gz"
    output:
        kraken_report=os.path.join(output_dir, "kraken/{id}/{id}.kraken_taxonomy.txt"),
        kraken_output=os.path.join(output_dir, "kraken/{id}/{id}.kraken_output.txt"),
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
       kraken_report_file= "output/kraken/{id}/{id}.kraken_taxonomy.txt"
    output:
        bracken_output= os.path.join(output_dir, "bracken/{id}.bracken_output.tsv")
    shell:
        """
            bracken -i {input.kraken_report_file} -d /workspace/Gene-pipeline/databases/k2/minikraken2_v2_8GB_201904_UPDATE -o {output.bracken_output}
        """
rule sample_validation:
    input: 
        bracken_output_file= "output/bracken/{id}.bracken_output.tsv"
        
