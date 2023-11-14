import csv
import yaml
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
#    comparisonGroupsList = [sample.get('comparisonGroup') for sample in file_path['Samples']]
    comparisonGroupDict = []
    for sample,sample_val in zip(file_path['Samples'],file_path['Samples'].values()):
        comparisonGroupDict.append((sample_val.get('comparisonGroup'),sample))
    return comparisonGroupDict




