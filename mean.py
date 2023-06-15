from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from CSV_Load import CSV_Loader
from scipy import stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def Mean_Median_Desc(filepath):
    data = CSV_Loader(filepath)

    descriptor_names = [desc[0] for desc in Descriptors.descList]
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)

    non_inhibitor_value = []
    inhibitor_value = []

    for index in range(len(data)):
        mol = Chem.MolFromSmiles(data["SMILES"][index])
        if data["ALDH1_inhibition"][index] == 0:
            non_inhib_desc_values = calc.CalcDescriptors(mol)
            non_inhibitor_value.append(non_inhib_desc_values)
        else:
            inhib_desc_values = calc.CalcDescriptors(mol)
            inhibitor_value.append(inhib_desc_values)

    mean_non_inhibitors = [sum(col) / len(col) for col in zip(*non_inhibitor_value)]
    mean_inhibitors = [sum(col) / len(col) for col in zip(*inhibitor_value)]

    # Perform t-test per descriptor
    ttest_results = []
    significant_count = 0
    super_significant_count = 0

    for i, descriptor in enumerate(descriptor_names):
        non_inhibitor = [desc[i] for desc in non_inhibitor_value]
        inhib_values = [desc[i] for desc in inhibitor_value]
        ttest_result = stats.ttest_ind(non_inhibitor, inhib_values)
        ttest_results.append(ttest_result)
        if ttest_result.pvalue < 0.05:
            significant_count += 1
        if ttest_result.pvalue < 0.01:
            super_significant_count += 1

    # Calculate median values
    median_non_inhibitors = [np.median(col) for col in zip(*non_inhibitor_value)]
    median_inhibitors = [np.median(col) for col in zip(*inhibitor_value)]

    
    # Create a DataFrame with the mean, median, and descriptor names
    df = pd.DataFrame({'Descriptor': descriptor_names,
                       'mean_non_inhibitors': mean_non_inhibitors,
                       'mean_inhibitors': mean_inhibitors,
                       'median_non_inhibitors': median_non_inhibitors,
                       'median_inhibitors': median_inhibitors,
                       'T-Statistic': [result.statistic for result in ttest_results],
                       'p-value': [result.pvalue for result in ttest_results],
                        })
    # Add significance columns
    df['Significance'] = df['p-value'].apply(
        lambda p: 'goed verschilletje hiero' if p < 0.05 else 'nope, deze niet')
    df['p-value < 0.05 Count'] = significant_count
    df['Super_Significance'] = df['p-value'].apply(
        lambda p: 'woow' if p < 0.01 else 'nope')
    df['p-value < 0.01 Count'] = super_significant_count


    # Lijstjes voor Jeroen:
    significant_descriptors = df[df['p-value'] < 0.05]['Descriptor'].tolist()
    super_significant_descriptors = df[df['p-value'] < 0.01]['Descriptor'].tolist()

    df.to_csv('means_median_stat_table.csv', index=False)

    # make second table with only the descriptors with significant difference 
    df2 = df [df['p-value']<0.01]
    df2 = df2[['Descriptor', 'mean_non_inhibitors','mean_inhibitors', 'median_non_inhibitors', 'median_inhibitors', 'p-value',]]
    df2.to_csv('descriptors_low_Pvalue.csv', index=False)

    return significant_descriptors, super_significant_descriptors
    
significant_descriptors, super_significant_descriptors = Mean_Median_Desc('tested_molecules_v2.csv')

print(len(super_significant_descriptors))

