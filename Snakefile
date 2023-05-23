import os
from snakemake.io import glob_wildcards

IDS, = glob_wildcards("fastq/{id}_1.fastq.gz")

print(IDS)
# # Define the input directory containing the FASTQ files
# input_file1 = "fastq/SRR7204427_R1.fastq.gz"
# input_file2 = "fastq/SRR7204427_R2.fastq.gz"

# Define the output directory
output_dir = "output"

# Define the final rule that specifies the targets to generate
rule all:
    input:
          expand("output/fastq/{id}/{id}.clean_1.fastq.gz", id=IDS), expand("output/fastq/{id}/{id}.clean_2.fastq.gz", id=IDS)

rule process_file_pair:
    input:
        fwd = "fastq/{id}_1.fastq.gz",
        rev = "fastq/{id}_2.fastq.gz"
    output:
        r1_clean = os.path.join(output_dir, "fastq/{id}/{id}.clean_1.fastq.gz"),
        r2_clean = os.path.join(output_dir, "fastq/{id}/{id}.clean_2.fastq.gz"),
        html = os.path.join(output_dir, "fastq/{id}/{id}.fastp.html"),
        json = os.path.join(output_dir, "fastq/{id}/{id}.fastp.json")
    shell:
        """
        fastp -i {input.fwd} -I {input.rev} \
              -o {output.r1_clean} -O {output.r2_clean} \
              -w 3 -h {output.html} -j {output.json}
        """





