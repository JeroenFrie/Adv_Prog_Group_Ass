from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from CSV_Load import CSV_Loader
from scipy import stats
import pandas as pd
import numpy as np
<<<<<<< HEAD
import matplotlib.pyplot as plt

=======

#import matplotlib.pyplot as plt
#import seaborn as sns
>>>>>>> a13ea30601e82abf22fccc873aff4006451760e8


def Mean_Median_Desc(filepath):
    data = CSV_Loader(filepath)
    
    descriptor_names = [desc[0] for desc in Descriptors.descList]
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)

    non_inhibitor_values = []
    inhibitor_values = []

    for index in range(len(data)):
        mol = Chem.MolFromSmiles(data["SMILES"][index])
        if data["ALDH1_inhibition"][index] == 0:
            non_inhib_desc_values = calc.CalcDescriptors(mol)
            non_inhibitor_values.append(non_inhib_desc_values)
        else:
            inhib_desc_values = calc.CalcDescriptors(mol)
            inhibitor_values.append(inhib_desc_values)
    
    # Convert the lists of descriptor values to NumPy arrays
    non_inhibitor_values = np.array(non_inhibitor_values)
    inhibitor_values = np.array(inhibitor_values)
    
    # Calculate the z-scores for the molecules
    non_inhibitor_z_scores = stats.zscore(non_inhibitor_values, axis=0)
    inhibitor_z_scores = stats.zscore(inhibitor_values, axis=0)
    
    # Identify outlier molecules based on a threshold
    threshold = 3
    non_inhibitor_outliers = np.abs(non_inhibitor_z_scores) > threshold
    inhibitor_outliers = np.abs(inhibitor_z_scores) > threshold
      
    # Replace outliers with 'nan'
    non_inhibitor_values[non_inhibitor_outliers] = np.nan
    inhibitor_values[inhibitor_outliers] = np.nan

    # Calculate the mean
    mean_non_inhibitors = np.nanmean(non_inhibitor_values, axis = 0)
    mean_inhibitors = np.nanmean(inhibitor_values, axis = 0)
    
    # Calculate median 
    median_non_inhibitors = np.nanmedian(non_inhibitor_values, axis =0)
    median_inhibitors = np.nanmedian(inhibitor_values, axis = 0)

    # Perform t-test per descriptor
    ttest_results = []
    significant_count = 0
    super_significant_count = 0

    for i, descriptor in enumerate(descriptor_names):
        non_inhibitor = non_inhibitor_values [:,i]
        inhibitor = inhibitor_values [:,i]
        
        #remove 'nan' values 
        non_inhibitor_clean = non_inhibitor [~ np.isnan(non_inhibitor)]
        inhibitor_clean = inhibitor [~ np.isnan(inhibitor)]
        
        ttest_result = stats.ttest_ind(non_inhibitor_clean, inhibitor_clean)
        ttest_results.append(ttest_result)
        if ttest_result.pvalue < 0.05:
            significant_count += 1
        if ttest_result.pvalue < 0.01:
            super_significant_count += 1

    
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
 
    # figure for the report
    fig,ax = plt.subplots()
    ax.plot(df2['Descriptor'], df2['p-value'], marker = 'o', linestyle= '')

    ax.set_xlabel('descriptors')
    ax.set_ylabel('p-value')
    ax.set_title('p-values plot')
    plt.xticks(rotation=90)
    plt.savefig('pvalue_plot.png')

    # table for the report
    df2.to_csv('p-value_table.csv', index=False)

    return significant_descriptors, super_significant_descriptors
    
<<<<<<< HEAD
significant_descriptors, super_significant_descriptors = Mean_Median_Desc('tested_molecules_v2.csv')

print(len(super_significant_descriptors))

=======
Mean_Median_Desc('tested_molecules_v3.csv')
>>>>>>> a13ea30601e82abf22fccc873aff4006451760e8
