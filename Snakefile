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


# Define the final rule that specifies the targets to generate
rule all:
    input:  
        expand(output_dir+"/fastq/{sample}/{sample}.clean_fwd.fastq.gz", sample = config["Samples"].keys()),
        expand(output_dir+"/fastq/{sample}/{sample}.clean_rev.fastq.gz", sample = config["Samples"].keys()),
        expand(output_dir+"/readSpeciesID/kraken/{sample}/{sample}.kraken_taxonomy.txt", sample = config["Samples"].keys()),
        expand(output_dir+"/readSpeciesID/kraken/{sample}/{sample}.kraken_output.txt", sample = config["Samples"].keys()),
        expand(output_dir+"/readSpeciesID/bracken/{sample}.bracken_output.tsv", sample = config["Samples"].keys()),
        expand(output_dir+"/readSpeciesID/sample_validation/{sample}.output.txt", sample = config["Samples"].keys()),
        expand(output_dir+"/assembly/comparisonGroup{comparisonGroup}/{sample}_assembly/contigs.fa", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/assemblySpeciesID/wgkb/comparisonGroup{comparisonGroup}/{sample}/{sample}.kraken_taxonomy.txt", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/assemblySpeciesID/wgkb/comparisonGroup{comparisonGroup}/{sample}/{sample}.kraken_output.txt", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/assemblySpeciesID/wgkb/comparisonGroup{comparisonGroup}/{sample}/{sample}.bracken_output.txt", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/assemblySpeciesID/wgv/{sample}/{sample}.validation_output.txt", sample = config["Samples"].keys()),
        expand(output_dir+"/assemblySpeciesID/mlst/comparisonGroup{comparisonGroup}/{sample}.contigs.mlst.tsv", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/assemblySpeciesID/rmlst/comparisonGroup{comparisonGroup}/{sample}.tsv", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/abricate/comparisonGroup{comparisonGroup}/card/{sample}.abricate.card.tsv", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/abricate/comparisonGroup{comparisonGroup}/vfdb/{sample}.abricate.vfdb.tsv", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/AMR/comparisonGroup{comparisonGroup}/{sample}.amrfinderplus.tsv", zip, comparisonGroup = [item[0] for item in comparisonGroupTuples], sample = [item[1] for item in comparisonGroupTuples]),
        expand(output_dir+"/cgMLST/comparisonGroup{comparisonGroup}/contigs_dir", comparisonGroup = [item[0] for item in comparisonGroupTuples]),
        expand(output_dir+"/cgMLST/comparisonGroup{comparisonGroup}/trainingFiles,trn_dir",comparisonGroup = [item[0] for item in comparisonGroupTuples]),
        expand(output_dir+"/cgMLST/comparisonGroup{comparisonGroup}/schema", comparisonGroup = [item[0] for item in comparisonGroupTuples]),
        expand(output_dir+"/cgMLST/comparisonGroup{comparisonGroup}/Allelecall",  comparisonGroup = [item[0] for item in comparisonGroupTuples]),
        expand(output_dir+"/cgMLST/comparisonGroup{comparisonGroup}/MST", comparisonGroup = [item[0] for item in comparisonGroupTuples]),
        expand(output_dir+"/cgMLST/comparisonGroup{comparisonGroup}/grapetree/cgMLST100.tre", comparisonGroup = [item[0] for item in comparisonGroupTuples])
      
include: "modules/preprocess.smk"
include: "modules/assembly.smk"
include: "modules/cgMLST.smk"