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
sampleToGroupDict = {}
for group, samples in cgSampleDict.items():
    for sample in samples:
        sampleToGroupDict[sample] = group
output_dir = "output"
krakenDB = "/workspace/WGS-pipeline/databases/k2"


# Define the final rule that specifies the targets to generate
rule all:
    input:  
        expand(output_dir+"/fastq/{sample}/{sample}.clean_fwd.fastq.gz", sample = config["Samples"].keys()),
        expand(output_dir+"/fastq/{sample}/{sample}.clean_rev.fastq.gz", sample = config["Samples"].keys()),
        expand(output_dir+"/fastq/{sample}/{sample}.fastp.html", sample = config["Samples"].keys()),
        expand(output_dir+"/fastq/{sample}/{sample}.fastp.json", sample = config["Samples"].keys()),
        expand(output_dir+"/readSpeciesID/kraken/{sample}/{sample}.kraken_taxonomy.txt", sample = config["Samples"].keys()),
        expand(output_dir+"/readSpeciesID/kraken/{sample}/{sample}.kraken_output.txt", sample = config["Samples"].keys()),
        expand(output_dir+"/readSpeciesID/bracken/{sample}.bracken_output.tsv", sample = config["Samples"].keys()),
        expand(output_dir+"/readSpeciesID/sample_validation/{sample}.output.txt", sample = config["Samples"].keys()),
        expand(output_dir+"/assembly/{sample}_assembly/contigs.fa", sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/assemblySpeciesID/wgkb/{sample}/{sample}.kraken_taxonomy.txt", sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/assemblySpeciesID/wgkb/{sample}/{sample}.kraken_output.txt",  sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/assemblySpeciesID/wgkb/{sample}/{sample}.bracken_output.txt",  sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/assemblySpeciesID/wgv/{sample}.validation_output.txt", sample = config["Samples"].keys()),
        expand(output_dir+"/assemblySpeciesID/mlst/{sample}.contigs.mlst.tsv", sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/assemblySpeciesID/rmlst/{sample}.tsv", sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/abricate/card/{sample}.abricate.card.tsv", sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/abricate/vfdb/{sample}.abricate.vfdb.tsv", sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/AMR/{sample}.amrfinderplus.tsv", sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/cgMLST/comparisonGroup{comparisonGroup}/contigs_dir", comparisonGroup = [item[0] for item in comparisonGroupTuples]),
        expand(output_dir+"/cgMLST/comparisonGroup{comparisonGroup}/trainingFiles/trn_file.trn",comparisonGroup = [item[0] for item in comparisonGroupTuples]),
        expand(output_dir+"/cgMLST/comparisonGroup{comparisonGroup}/schema", comparisonGroup = [item[0] for item in comparisonGroupTuples]),
        expand(output_dir+"/cgMLST/comparisonGroup{comparisonGroup}/Allelecall",  comparisonGroup = [item[0] for item in comparisonGroupTuples]),
        expand(output_dir+"/cgMLST/comparisonGroup{comparisonGroup}/MST", comparisonGroup = [item[0] for item in comparisonGroupTuples]),
        expand(output_dir+"/cgMLST/comparisonGroup{comparisonGroup}/grapetree/cgMLST100.tre", comparisonGroup = [item[0] for item in comparisonGroupTuples])
      
include: "modules/preprocess.smk"
include: "modules/assembly.smk"
include: "modules/cgMLST.smk"