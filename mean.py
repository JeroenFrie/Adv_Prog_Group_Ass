from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
import pandas as pd
from CSV_Load import CSV_Loader
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

# Create a DataFrame with the mean values and descriptor names
df = pd.DataFrame({'Descriptor': descriptor_names, 'mean_non_inhibitors': mean_non_inhibitors, 'mean_inhibitors': mean_inhibitors,'T-Statistic': [result.statistic for result in ttest_results], 'p-value': [result.pvalue for result in ttest_results]})
df['Significance'] = df['p-value'].apply(lambda p: 'jaaaaaaaaaaaaa, goed verschilletje hiero' if p < 0.05 else 'nope, deze niet')
df['p-value < 0.05 Count'] = significant_count
df['Super_Significance'] = df['p-value'].apply(lambda p: 'woooooooooooooooooooooooooooooooooooooooooooooooow' if p < 0.01 else 'nope')
df['p-value < 0.05 Count'] = super_significant_count


# Sort the DataFrame based on p-values in ascending order
df_sorted = df.sort_values(by='p-value')

# Save the sorted DataFrame as a CSV file
df_sorted.to_csv('ttest_results_sorted.csv', index=False)
