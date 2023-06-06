from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from CSV_Load import CSV_Loader
from scipy import stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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

    # Perform t-test per descriptor for means
    ttest_results_mean = []
    significant_count_mean = 0
    super_significant_count_mean = 0

    for i, descriptor in enumerate(descriptor_names):
        mean_non_inhib = [desc[i] for desc in non_inhibitor_value]
        mean_inhib_values = [desc[i] for desc in inhibitor_value]
        ttest_result_mean = stats.ttest_ind(mean_non_inhib, mean_inhib_values)
        ttest_results_mean.append(ttest_result_mean)
        if ttest_result_mean.pvalue < 0.05:
            significant_count_mean += 1
        if ttest_result_mean.pvalue < 0.01:
            super_significant_count_mean += 1

    # Calculate median values
    median_non_inhibitors = [np.median(col) for col in zip(*non_inhibitor_value)]
    median_inhibitors = [np.median(col) for col in zip(*inhibitor_value)]

    # Perform t-test per descriptor for medians
    ttest_results_median = []
    significant_count_median = 0
    super_significant_count_median = 0
    for i, descriptor in enumerate(descriptor_names):
        median_non_inhibitors = [desc[i] for desc in non_inhibitor_value]
        median_inhibitors = [desc[i] for desc in inhibitor_value]
        median_ttest_result = stats.ttest_ind(median_non_inhibitors, median_inhibitors)
        ttest_results_median.append(median_ttest_result)
        if ttest_result_median.pvalue < 0.05:
            significant_count_median += 1
        if ttest_result_median.pvalue < 0.01:
            super_significant_count_median += 1


    # Create a DataFrame with the mean, median, and descriptor names
    df = pd.DataFrame({'Descriptor': descriptor_names,
                       'mean_non_inhibitors': mean_non_inhibitors,
                       'mean_inhibitors': mean_inhibitors,
                       'median_non_inhibitors': median_non_inhibitors,
                       'median_inhibitors': median_inhibitors,
                       'T-Statistic Mean': [result.statistic for result in ttest_results_mean],
                       'p-value Mean': [result.pvalue for result in ttest_results_mean],
                       'T-Statistic Median': [result.statistic for result in median_ttest_results],
                       'p-value Median': [result.pvalue for result in median_ttest_results]})

    # Add significance columns
    df['Significance Mean'] = df['p-value Mean'].apply(
        lambda p: 'jaaaaaaaaaaaaa, goed verschilletje hiero' if p < 0.05 else 'nope, deze niet')
    df['p-value < 0.05 Count mean'] = significant_count_mean
    df['Super_Significance Mean'] = df['p-value Mean'].apply(
        lambda p: 'woooooooooooooooooooooooooooooooooooooooooooooooow' if p < 0.01 else 'nope')
    df['p-value < 0.01 Count'] = super_significant_count_mean

    df['Significance Median'] = df['p-value Median'].apply(
        lambda p: 'jaaaaaaaaaaaaa, goed verschilletje hiero' if p < 0.05 else 'nope, deze niet')
    df['p-value < 0.05 Count median'] = significant_count_median
    df['Super_Significance Median'] = df['p-value Median'].apply(
        lambda p: 'woooooooooooooooooooooooooooooooooooooooooooooooow' if p < 0.01 else 'nope')
    df['p-value < 0.01 Count'] = super_significant_count_median

    # Lijstjes voor Jeroen:
    super_significant_descriptors_mean_list = df[df['p-value Mean'] < 0.01]['Descriptor'].tolist()
    super_significant_descriptors_median_list = df[df['p-value Median'] < 0.01]['Descriptor'].tolist()

    df.to_csv('means_table.csv', index=False)

    # Boxplots for mean
    # Create a DataFrame with the descriptor values for significant descriptors
    significant_descriptor_values = pd.DataFrame({desc: [desc_val[i] for desc_val in non_inhibitor_value + inhibitor_value] for i, desc in enumerate(descriptor_names) if desc in super_significant_descriptors_mean_list})

    # Create boxplots for the significant descriptors
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=significant_descriptor_values)
    plt.xticks(rotation=90)
    plt.xlabel('Descriptor')
    plt.ylabel('Descriptor Value')
    plt.title('Boxplots for P<0.01 Mean')
    plt.tight_layout()

    # Save the boxplots to a file
    plt.savefig('boxplots_mean.png')

    # Write the descriptor values to a CSV file
    significant_descriptor_values.to_csv('descriptor_values_boxplots_mean.csv', index=False)

    # Now for the median values
    significant_descriptor_values_median = pd.DataFrame({desc: [desc_val[i] for desc_val in non_inhibitor_value + inhibitor_value] for i, desc in enumerate(descriptor_names) if desc in super_significant_descriptors_median_list})

    plt.figure(figsize=(12, 6))
    sns.boxplot(data=significant_descriptor_values_median)
    plt.xticks(rotation=90)
    plt.xlabel('Descriptor')
    plt.ylabel('Descriptor Value')
    plt.title('Boxplots for P<0.01 Median')
    plt.tight_layout()

    # Save the boxplots to a file
    plt.savefig('boxplots_median.png')

    # Write the descriptor values to a CSV file
    significant_descriptor_values.to_csv('descriptor_values_boxplots_median.csv', index=False)


    return super_significant_descriptors_mean_list, super_significant_descriptors_median_list
    
