from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
import pandas as pd
from CSV_Load import CSV_Loader


data = CSV_Loader("tested_molecules-1.csv")

descriptor_names = [desc[0] for desc in Descriptors.descList]
calc = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)

non_inhibitor_value = []
inhibitor_value =[]

for index in range(len(data)):
    mol = Chem.MolFromSmiles(data["SMILES"][index])
    if data["ALDH1_inhibition"][index] == 0:
        non_inhib_desc_values= calc.CalcDescriptors(mol)
        non_inhibitor_value.append(non_inhib_desc_values)
    else:
        inhib_desc_values = calc.CalcDescriptors(mol)
        inhibitor_value.append(inhib_desc_values)

mean_non_inhibitors = [sum(col) / len(col) for col in zip(*non_inhibitor_value)]
mean_inhibitors = [sum(col) / len(col) for col in zip(*inhibitor_value)]

with open('means.txt', 'w') as file: 
    file.write('mean_non_inhibitors=' + str(mean_non_inhibitors) + '\n')
    file.write('mean_inhibitors=' + str(mean_inhibitors) + '\n')