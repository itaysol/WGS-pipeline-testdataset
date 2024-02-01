import os
import shutil
import sys
import ast
import json
import re

def extract_contigs(mapping,sampleToGroupDict):
    for input_file in mapping:
        # Extract the sample name from the input file path
        print(f"Input File: {input_file}")
        sample = input_file.split('_')[0].split('/')[2]
        print(f"Extracted Sample: {sample}")

        # Determine the comparison group directory based on the sample name
        comparison_group_dir = os.path.join("output","cgMLST",f'comparisonGroup{sampleToGroupDict[sample]}',"contigs_dir")

        # Create the output directory if it doesn't exist
        os.makedirs(comparison_group_dir, exist_ok=True)

        # Copy the file to the output directory and rename it
        output_contigs_file = os.path.join(comparison_group_dir, f'{sample}.contigs.fa')
        shutil.copy(input_file, output_contigs_file)
        print(f"Copied {input_file} to {output_contigs_file}")

def create_input_list():
    input_list = []
    for i in range(1,len(sys.argv)):
        input_list.append(sys.argv[i])
    return input_list

if __name__ == '__main__':
    mapping = create_input_list()
    print(os.environ['SAMPLE_TO_GROUP_DICT'])

    # Load the sample_to_group_dict from the environment variable
    sample_to_group_dict_str = os.environ.get('SAMPLE_TO_GROUP_DICT', '{}')

    # Add quotes around keys and values
    sample_to_group_dict_str_fixed = re.sub(r'([a-zA-Z0-9_]+):', r'"\1":', sample_to_group_dict_str)

    # Convert the string to a dictionary using JSON parsing
    sample_to_group_dict = json.loads(sample_to_group_dict_str_fixed)
    extract_contigs(mapping, sample_to_group_dict)
