import os
from snakemake.io import glob_wildcards

IDS, = glob_wildcards("fastq/{id}_R1.fastq.gz")


# # Define the input directory containing the FASTQ files
# input_file1 = "fastq/SRR7204427_R1.fastq.gz"
# input_file2 = "fastq/SRR7204427_R2.fastq.gz"

# Define the output directory
output_dir = "output"

# Define the final rule that specifies the targets to generate
rule all:
    input:
          expand("fastq/{id}_R1.fastq.gz", id=IDS), expand("fastq/{id}_R2.fastq.gz", id=IDS)

rule process_file_pair:
    input:
        "fastq/{id}_R1.fastq.gz", "fastq/{id}_R2.fastq.gz"
    output:
        r1_clean = os.path.join(output_dir, "fastq/{id}_R1.fastq.gz"),
        r2_clean = os.path.join(output_dir, "fastq/{id}_R2.fastq.gz"),
        html = os.path.join(output_dir, "fastq/{id}.fastp.html")
    shell:
        """
        fastp -i {"fastq/{id}_R1.fastq.gz"} -I {"fastq/{id}_R2.fastq.gz"} \
              -o {output.r1_clean} -O {output.r2_clean} \
              -w 3 -h {output.html}
        """





