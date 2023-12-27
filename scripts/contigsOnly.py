import os
import shutil
import sys

def copy_contigs(input_dir, output_dir):
   def copy_contigs(input_dir, output_dir):
    # List all subdirectories in the input directory (assumed to be sample assemblies)
    samples = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]

    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Copy contigs.fa from each sample assembly to the output directory
    for sample in samples:
        assembly_dir = os.path.join(input_dir, sample)
        contigs_file = os.path.join(assembly_dir, 'contigs.fa')

        # Check if the contigs.fa file exists for the current sample
        if os.path.exists(contigs_file):
            # Copy contigs.fa to the output directory, naming it after the sample
            output_contigs_file = os.path.join(output_dir, f'{sample}_contigs.fa')
            shutil.copy(contigs_file, output_contigs_file)
            print(f"Copied {contigs_file} to {output_contigs_file}")
        else:
            print(f"contigs.fa not found for {sample}")

# Example usage
if __name__ == '__main__':
    input_directory = sys.argv[1]
    output_directory = sys.argv[2]
    copy_contigs(input_directory, output_directory)
