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

rule create_contig_dirs:
    input:
        assembly = lambda wildcards: expand(os.path.join(output_dir, "assembly", f"{wildcards.comparisonGroup}", "{sample}_assembly", "contigs.fa"), sample=cgSampleDict[f"{wildcards.comparisonGroup}"[-1]]),
    output:
        contigs_dir = directory(os.path.join(output_dir,"cgMLST","{comparisonGroup}","contigs_dir"))
    shell:
        """
        python scripts/extract_contigs.py {input.assembly[0]} {output.contigs_dir} 
        
        """

rule create_training_file:
    conda:
        "env/conda-prodigal.yaml"
    params:
        specie = lambda wildcards: cgSpecieDict[wildcards.comparisonGroup[-1]],
        specie_no_spaces = lambda wildcards: cgSpecieDict[wildcards.comparisonGroup[-1]].replace(" ","") 
    output:
        specie_zip = directory(os.path.join(output_dir,"cgMLST","{comparisonGroup}","trainingFiles")),
        training_file = os.path.join(output_dir,"cgMLST","{comparisonGroup}","trainingFiles","trn_file.trn")
    shell:
        """
        datasets download genome taxon "{params.specie}"  --reference --include genome --filename {output.specie_zip}/{params.specie_no_spaces}.zip
        unzip {output.specie_zip}/{params.specie_no_spaces}.zip -d {output.specie_zip}
        
        # Extract the path of the FNA file from dataset_catalog.json
        fna_file_path=$(cat {output.specie_zip}/ncbi_dataset/data/dataset_catalog.json | jq -r '.assemblies[].files[] | select(.fileType == "GENOMIC_NUCLEOTIDE_FASTA") | .filePath')
        
        # Run prodigal using the extracted FNA file path
        prodigal -i {output.specie_zip}/ncbi_dataset/data/$fna_file_path -t {output.training_file} -p single
        
        """    


rule adhoc_cgMLST:
    conda:
        "env/conda-chewBBACA.yaml"
    input:
        contigs_dir = os.path.join(output_dir,"cgMLST","{comparisonGroup}","contigs_dir"),
        training_file = os.path.join(output_dir,"cgMLST","{comparisonGroup}","trainingFiles","trn_file.trn")
    output:
        createSchema = directory(os.path.join(output_dir,"cgMLST","{comparisonGroup}","schema")),
        allelecall = directory(os.path.join(output_dir,"cgMLST","{comparisonGroup}","Allelecall")),
        MST = directory(os.path.join(output_dir,"cgMLST","{comparisonGroup}","MST"))
    shell:
      """

        chewBBACA.py CreateSchema -i {input.contigs_dir} -o {output.createSchema} --ptf {input.training_file} --cpu 20

        # Loop until the condition is met
        until chewBBACA.py AlleleCall -i {input.contigs_dir} -g {output.createSchema}/schema_seed/ -o {output.allelecall} --cpu 14; do
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
        grapetree_input = os.path.join(output_dir,"cgMLST","{comparisonGroup}","MST")
    output:
        grapetree_output = os.path.join(output_dir,"cgMLST","{comparisonGroup}","grapetree","cgMLST100.tre")
    shell:
        """
        grapetree --profile {input.grapetree_input}/cgMLST100.tsv > {output.grapetree_output}

        """

    