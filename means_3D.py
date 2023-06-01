# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 14:37:50 2023

@author: 20203129
"""
# for 2D and 3D enzo
from rdkit import Chem
from rdkit.Chem import Descriptors, Descriptors3D, AllChem
from rdkit.ML.Descriptors import MoleculeDescriptors
from CSV_Load import CSV_Loader
from scipy import stats
import pandas as pd
import numpy as np

data = CSV_Loader("tested_molecules-1.csv")
num_conformers = 5
# Create a list of descriptor names
# descriptor_names_3d = ['TPSA', 'Asphericity', 'Eccentricity', 
#                        'InertialShapeFactor', 'NPR2', 'PMI1', 'PMI2', 
#                        'PMI3', 'RadiusOfGyration', 'SpherocityIndex']
descriptor_names_3d = ['Asphericity']
descriptor_names = [desc[0] for desc in Descriptors.descList]

# Create a descriptor calculator 
calc = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)

non_inhibitor_value = []
inhibitor_value =[]

for index in range(len(data)):
    mol = Chem.AddHs(Chem.MolFromSmiles(data["SMILES"][index]))
    if data["ALDH1_inhibition"][index] == 0:
        # Make tuple with all 2d descriptors:
        non_inhib_desc_values= calc.CalcDescriptors(mol)
        
        # Calculate 3d descriptors:
        AllChem.EmbedMultipleConfs(mol, num_conformers)
        driedee_list = []
        asphericities = []
        for conf_id in range(num_conformers):
            conf = mol.GetConformer(conf_id)
            asphericity = Descriptors3D.Asphericity(mol, confId=conf_id)
            asphericities.append(asphericity)
            # Calculate the mean asphericity
            mean_asphericity = sum(asphericities) / num_conformers
            driedee_list.append(mean_asphericity)
            
            driedee_tuple = tuple(driedee_list)
        non_inhibitor_value.append(non_inhib_desc_values+driedee_tuple)
    else:
        inhib_desc_values = calc.CalcDescriptors(mol)
        # Calculate 3d descriptors:
        AllChem.EmbedMultipleConfs(mol, num_conformers)
        driedee_list = []
        asphericities = []
        for conf_id in range(num_conformers):
            conf = mol.GetConformer(conf_id)
            asphericity = Descriptors3D.Asphericity(mol, confId=conf_id)
            asphericities.append(asphericity)
            # Calculate the mean asphericity
            mean_asphericity = sum(asphericities) / num_conformers
            driedee_list.append(mean_asphericity)
            
            driedee_tuple = tuple(driedee_list)
            inhibitor_value.append(inhib_desc_values+driedee_tuple)
        
        
# for index in range(len(data)):
#     mol = Chem.MolFromSmiles(data["SMILES"][index])
#     AllChem.EmbedMultipleConfs(mol, num_conformers)
#     # Calculate the asphericity for each conformer
#     asphericities = []
#     for conf_id in range(num_conformers):
#         conf = mol.GetConformer(conf_id)
#         asphericity = Descriptors3D.Asphericity(mol, confId=conf_id)
#         asphericities.append(asphericity)
#         # Calculate the mean asphericity
#         mean_asphericity = sum(asphericities) / num_conformers
#     if data["ALDH1_inhibition"][index] == 0:
#         non_inhibitor_value.append(mean_asphericity)
#     else:
#         inhibitor_value.append(mean_asphericity)
        
mean_non_inhibitors = [sum(col) / len(col) for col in zip(*non_inhibitor_value)]
mean_inhibitors = [sum(col) / len(col) for col in zip(*inhibitor_value)]

# Perform t-test per descriptor
mean_ttest_results = []
significant_count = 0 
super_significant_count = 0
for i, descriptor in enumerate(descriptor_names+descriptor_names_3d):
    non_inhib_values = [desc[i] for desc in non_inhibitor_value]
    inhib_values = [desc[i] for desc in inhibitor_value]
    mean_ttest_result = stats.ttest_ind(non_inhib_values, inhib_values)
    mean_ttest_results.append(mean_ttest_result)
    if mean_ttest_result.pvalue < 0.05:
        significant_count += 1
    if mean_ttest_result.pvalue < 0.01:
        super_significant_count += 1

# Perform t-test for medians
median_ttest_results = []
for i, descriptor in enumerate(descriptor_names):
    non_inhib_values = [desc[i] for desc in non_inhibitor_value]
    inhib_values = [desc[i] for desc in inhibitor_value]
    median_ttest_result = stats.ttest_ind(non_inhib_values, inhib_values)
    median_ttest_results.append(median_ttest_result)
# Calculate median values
median_non_inhibitors = [np.median(col) for col in zip(*non_inhibitor_value)]
median_inhibitors = [np.median(col) for col in zip(*inhibitor_value)]

# Create a DataFrame with the mean, median, and descriptor names
df = pd.DataFrame({'Descriptor': descriptor_names,
                   'mean_non_inhibitors': mean_non_inhibitors,
                   'mean_inhibitors': mean_inhibitors,
                   'median_non_inhibitors': median_non_inhibitors,
                   'median_inhibitors': median_inhibitors,
                   'T-Statistic Mean': [result.statistic for result in mean_ttest_results],
                   'p-value Mean': [result.pvalue for result in mean_ttest_results],
                   'T-Statistic Median': [result.statistic for result in median_ttest_results],
                   'p-value Median': [result.pvalue for result in median_ttest_results]})

# Add significance columns
df['Significance Mean'] = df['p-value Mean'].apply(lambda p: 'jaaaaaaaaaaaaa, goed verschilletje hiero' if p < 0.05 else 'nope, deze niet')
df['p-value < 0.05 Count'] = significant_count
df['Super_Significance Mean'] = df['p-value Mean'].apply(lambda p: 'woooooooooooooooooooooooooooooooooooooooooooooooow' if p < 0.01 else 'nope')
df['p-value < 0.01 Count'] = super_significant_count
df['Significance Median'] = df['p-value Median'].apply(lambda p: 'jaaaaaaaaaaaaa, goed verschilletje hiero' if p < 0.05 else 'nope, deze niet')

# Lijstjes voor Jeroen:
super_significant_descriptors_mean_list = df[df['p-value Mean'] < 0.01]['Descriptor'].tolist()
super_significant_descriptors_median_list = df[df['p-value Median'] < 0.01]['Descriptor'].tolist()

df.to_csv('means_table_3D.csv', index=False)

