from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
import pandas as pd
from CSV_Load import CSV_Loader
import numpy as np
from scipy import stats

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

# Perform t-test per descriptor
ttest_results = []
significant_count = 0 
super_significant_count = 0
for i, descriptor in enumerate(descriptor_names):
    non_inhib_values = [desc[i] for desc in non_inhibitor_value]
    inhib_values = [desc[i] for desc in inhibitor_value]
    ttest_result = stats.ttest_ind(non_inhib_values, inhib_values)
    ttest_results.append(ttest_result)
    if ttest_result.pvalue < 0.05:
        significant_count += 1
    if ttest_result.pvalue < 0.01:
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
                   'T-Statistic Mean': [result.statistic for result in ttest_results],
                   'p-value Mean': [result.pvalue for result in ttest_results],
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

df.to_csv('means_table.csv', index=False)


