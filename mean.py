from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from CSV_Load import CSV_Loader
from scipy import stats
import pandas as pd
import numpy as np

#import matplotlib.pyplot as plt
#import seaborn as sns

def remove_outliers_zscore (data, threshold=3):
     z_scores = stats.zscore(data)
     outliers = np.abs(z_scores) > threshold
     clean_data = data [outliers == False]
     return clean_data

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
    super_significant_descriptors_mean_list = df[df['p-value'] < 0.01]['Descriptor'].tolist()
    super_significant_descriptors_median_list = df[df['p-value'] < 0.01]['Descriptor'].tolist()

    df.to_csv('means_median_stat_table.csv', index=False)

    # make second table with only the descriptors with significant difference 
    df2 = df [df['p-value']<0.01]
    df2 = df2[['Descriptor', 'mean_non_inhibitors','mean_inhibitors', 'median_non_inhibitors', 'median_inhibitors', 'p-value',]]
    df2.to_csv('descriptors_low_Pvalue.csv', index=False)

    # Boxplots for mean
    # Create a DataFrame with the descriptor values for significant descriptors
   # significant_descriptor_values = pd.DataFrame({desc: [desc_val[i] for desc_val in non_inhibitor_value + inhibitor_value] for i, desc in enumerate(descriptor_names) if desc in super_significant_descriptors_mean_list})

    # Create boxplots for the significant descriptors
 #   plt.figure(figsize=(12, 6))
 #   sns.boxplot(data=significant_descriptor_values)
 #   plt.xticks(rotation=90)
 #   plt.xlabel('Descriptor')
 #   plt.ylabel('Descriptor Value')
 #   plt.title('Boxplots for P<0.01 Mean')
 #   plt.tight_layout()

    # Save the boxplots to a file
  #  plt.savefig('boxplots_mean.png')

    # Write the descriptor values to a CSV file
  #  significant_descriptor_values.to_csv('descriptor_values_boxplots_mean.csv', index=False)

    # Now for the median values
  #  significant_descriptor_values_median = pd.DataFrame({desc: [desc_val[i] for desc_val in non_inhibitor_value + inhibitor_value] for i, desc in enumerate(descriptor_names) if desc in super_significant_descriptors_median_list})

  #  plt.figure(figsize=(12, 6))
  #  sns.boxplot(data=significant_descriptor_values_median)
  #  plt.xticks(rotation=90)
  #  plt.xlabel('Descriptor')
  #  plt.ylabel('Descriptor Value')
  #  plt.title('Boxplots for P<0.01 Median')
  #  plt.tight_layout()

    # Save the boxplots to a file
  #  plt.savefig('boxplots_median.png')

    # Write the descriptor values to a CSV file
   # significant_descriptor_values.to_csv('descriptor_values_boxplots_median.csv', index=False)


    return super_significant_descriptors_mean_list, super_significant_descriptors_median_list
    
Mean_Median_Desc('tested_molecules-1.csv')