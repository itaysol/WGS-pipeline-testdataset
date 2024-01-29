import os
import shutil
import sys

def extract_contigs(input_directory, output_directory):
    # Iterate over files in the input directory
    for sample_file in os.listdir(input_directory):
        # Create the full path to the input file
        input_file_path = os.path.join(input_directory, sample_file, "contigs.fa")
        
        # Check if it's a file and ends with contigs.fa
        if os.path.isfile(input_file_path):
            # Extract the sample name from the file name
            sample = sample_file.split('_')[0]

            # Create the output directory if it doesn't exist
            output_group_dir = os.path.join(output_directory)
            os.makedirs(output_group_dir, exist_ok=True)

            # Copy the file to the output directory and rename it
            output_contigs_file = os.path.join(output_group_dir, f'{sample}.contigs.fa')
            shutil.copy(input_file_path, output_contigs_file)
            print(f"Copied {input_file_path} to {output_contigs_file}")

# Example usage
if __name__ == '__main__':
    # Pass the input and output directories as command line arguments
    input_directory = sys.argv[1][:32]
    output_directory = sys.argv[2]
    extract_contigs(input_directory, output_directory)
