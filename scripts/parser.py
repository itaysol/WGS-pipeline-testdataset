import csv
import yaml
import itertools
import pandas as pd
# determine the type of the config file
def determineConfigFileType(file_path):
    type = file_path.split(".")[1]
    if type == "yaml" or type == "csv":
        return type
    else:
        print("Invalid type of config file. Please enter a csv or yaml file")
        return None

# turn csv files to yaml files
def turnCsvToYaml(file_path):
    csv_file = file_path
    yaml_file = 'user_config_file.yaml'
    data = {}
    with open(csv_file, 'r') as csvfile:
        csvreader = csv.DictReader(csvfile)
        for row in csvreader:
            sample = row['sample']
            if sample not in data:
                data[sample] = {}
            for key, value in row.items():
                if key != 'sample':
                    data[sample][key] = value

    output = {'__use_yte__': True, 'Samples': data}
    with open(yaml_file, 'w') as yamlfile:
        yaml.dump(output, yamlfile)
    return yaml_file

def setConfigFile(file_path):
    type = determineConfigFileType(file_path)
    if type == None:
        exit()
    elif type == "csv" :
         return turnCsvToYaml(file_path)
    else:
        return file_path

def parser(file_path):
    return setConfigFile(file_path)

def setComparisonGroups(file_path):
    comparisonGroupList = []
    comparisonGroupDict = {}
    comparisonGroupSpecies = {}
    comparisonGroupSampCounter = {}
    genomeSizeDict = {}
    cgSampleDict = {}

    for sample, sample_val in zip(file_path['Samples'], file_path['Samples'].values()):
        comparisonGroup = sample_val.get('comparisonGroup')
        specie = sample_val.get('specie')
        if comparisonGroup in cgSampleDict:
            cgSampleDict[comparisonGroup].append(sample)
        else:
            cgSampleDict[comparisonGroup]=[sample]
        if comparisonGroup in comparisonGroupSampCounter:
            comparisonGroupSampCounter[comparisonGroup]+=1
        else:
            comparisonGroupSampCounter[comparisonGroup]=1

        # Check if the comparison group already has a species assigned
        if comparisonGroup in comparisonGroupSpecies:
            existingSpecie = comparisonGroupSpecies[comparisonGroup]
            if existingSpecie != specie:
                raise ValueError(f"Error: Inconsistent species in comparison group {comparisonGroup}")
        else:
            comparisonGroupSpecies[comparisonGroup] = specie

        comparisonGroupList.append((comparisonGroup, sample))
        comparisonGroupDict[comparisonGroup] = specie
    
    df2 = pd.read_csv("https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/overview.txt",sep="\t",low_memory=False)
    df2_bacteria = df2[df2['Kingdom'].isin(['Bacteria'])]
    df2_bacteria_of_interest = df2_bacteria[df2_bacteria['#Organism/Name'].isin(comparisonGroupSpecies.values())]
    for (spec,size) in zip(df2_bacteria_of_interest['#Organism/Name'], df2_bacteria_of_interest['Size (Mb)']):
       genomeSizeDict[spec] = size
    return comparisonGroupList, comparisonGroupDict, comparisonGroupSampCounter, genomeSizeDict, cgSampleDict







