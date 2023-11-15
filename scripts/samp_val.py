import os
import pandas as pd
import sys
# BRACKEN_TXT = r"bracken.txt"


def main(input_file, organism, sampleid, output_path) -> None:
    current_directory = os.getcwd()
    file_name = output_path
    output_path = os.path.join(current_directory, file_name)
    data = pd.read_csv(input_file, delimiter='\t')
    data = data.iloc[:, [0, -1]]
    other_species_criteria = True
    max_index = data['fraction_total_reads'].idxmax()  # row of the organism with max abundance
    max_name = data.loc[max_index, 'name']  # name of the organism with max abundance
    max_abundance = data.loc[max_index, 'fraction_total_reads']
    my_specie_criteria = (max_name == organism) and (max_abundance >= 0.7)
    selected_rows = []
    for index, row in data.iterrows():
        if row['fraction_total_reads'] > 0.06 and row['name'] != organism:
            selected_rows.append(row)
            other_species_criteria = False  # found a specie with > 6% abundance
    my_row = data[data['name'] == organism]
    my_abundance = my_row['fraction_total_reads']
    my_row = pd.DataFrame(data[data['name'] == organism])
    other_rows = pd.DataFrame(selected_rows)
    if my_specie_criteria and other_species_criteria:
        with open(output_path, 'w') as output_file:
            output_file.write(sampleid + '. PASSED.\n\n')
    elif not my_abundance.empty and my_abundance < 0.7:
        with open(output_path, 'w') as output_file:
            output_file.write(sampleid + '. FAILED. ' + organism + ' abundance is less than 0.7\n\n')
    elif not other_species_criteria:
        with open(output_path, 'w') as output_file:
            output_file.write(sampleid + '. FAILED. Other species have abundance > 0.06\n\n')
    else:
        with open(output_path, 'w') as output_file:
            output_file.write(sampleid + '. FAILED.\n\n')
    my_row.to_csv(output_path, mode='a', index=False, sep='\t')
    # append the second DataFrame to the output file (append mode)
    other_rows.to_csv(output_path, mode='a', header=False, index=False, sep='\t')
    return output_path

if __name__ == '__main__':
    if len(sys.argv) < 4:
        raise Exception("Invalid number of arguments. Expected Bracken file, User specie and Sample ID.")
    BRACKEN_TXT = sys.argv[2]
    specie = sys.argv[5]+" "+ sys.argv[6]
    sample_id = sys.argv[7]
    output_path = sys.argv[4]
    main(BRACKEN_TXT, specie, sample_id, output_path)



