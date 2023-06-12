# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 14:37:50 2023

@author: 20203129
"""
# for 2D and 3D
from rdkit import Chem
from rdkit.Chem import Descriptors, Descriptors3D, AllChem
from rdkit.ML.Descriptors import MoleculeDescriptors
from CSV_Load import CSV_Loader
from scipy import stats
import pandas as pd
import numpy as np

# Descriptor calculator for 3D
def drieD_descriptors (num_conformers, mol):
    drieD_list = []
    # List for every descriptor
    asphericities = []
    eccentricities = []
    inertialShapeFactors = []
    NPR1s = []
    NPR2s = []
    PMI1s =[]
    PMI2s =[]
    PMI3s =[]
    RadiusOfGyrations = []
    SpherocityIndexs =[]
    
    # Iteratate for every conformation
    for conf_id in range(num_conformers):
        conf = mol.GetConformer(conf_id)
        
        # Calculate value of every descriptor
        asphericity     = Descriptors3D.Asphericity(mol, confId=conf_id)
        eccentricity    = Descriptors3D.Eccentricity(mol, confId=conf_id)
        inertialShapeFactor = Descriptors3D.InertialShapeFactor(mol, confId=conf_id)
        NPR1 = Descriptors3D.NPR1(mol, confId=conf_id)
        NPR2 = Descriptors3D.NPR2(mol, confId=conf_id)
        PMI1 = Descriptors3D.PMI1(mol, confId=conf_id)
        PMI2 = Descriptors3D.PMI2(mol, confId=conf_id)
        PMI3 = Descriptors3D.PMI3(mol, confId=conf_id)
        RadiusOfGyration = Descriptors3D.RadiusOfGyration(mol, confId=conf_id)
        SpherocityIndex  = Descriptors3D.SpherocityIndex(mol, confId=conf_id)
        
        # Add descriptor value to list
        asphericities.append(asphericity)        
        eccentricities.append(eccentricity)   
        inertialShapeFactors.append(inertialShapeFactor)
        NPR1s.append(NPR1)
        NPR2s.append(NPR2)
        PMI1s.append(PMI1)
        PMI2s.append(PMI2)
        PMI3s.append(PMI3)
        RadiusOfGyrations.append(RadiusOfGyration)
        SpherocityIndexs.append(SpherocityIndex)
        
    # Calculate the mean per descriptor for the total number of conformations
    drieD_descriptors_list = [asphericities, eccentricities, 
                              inertialShapeFactors, NPR1s, NPR2s, 
                              PMI1s, PMI2s, PMI3s, RadiusOfGyrations,
                              SpherocityIndexs ]
    for descriptor in drieD_descriptors_list:
        mean = sum(descriptor) / num_conformers
        drieD_list.append(mean)  
    return drieD_list

<<<<<<< HEAD
non_inhibitor_value = []
inhibitor_value =[]

info = ({
    'SMILES':[],
    'ALDH1_inhibition' :[],
    'Asphericity':[],
    'Eccentricity':[],
    'InertialShapeFactor':[],
    'NPR1':[],
    'NPR2':[],
    'PMI1':[],
    'PMI2':[],
    'PMI3':[],
    'RadiusOfGyration':[],
    'SpherocityIndex':[]
               })
df_3d_descriptors  = pd.DataFrame(info)

for index in range(len(data)):
    mol = Chem.AddHs(Chem.MolFromSmiles(data["SMILES"][index])) #AddHs adds the hydrogen atoms
   
    #non inhibitor;
    if data["ALDH1_inhibition"][index] == 0:     
        # Make tuple with all 2d descriptors:
        non_inhib_desc_values= calc.CalcDescriptors(mol)
        
        # Calculate 3d descriptors:
        AllChem.EmbedMultipleConfs(mol, num_conformers)
        driedee_list = drieD_descriptors(num_conformers, mol)        
        
        # Make tuple with all 3d descriptors
        driedee_tuple = tuple(driedee_list)
        # Add both tuples to non inhibotor list
        non_inhibitor_value.append(non_inhib_desc_values + driedee_tuple)
        
        #Add 3D descriptor values to the 3D dataframe 
        info_list = [data["SMILES"][index],0]+driedee_list
        df_3d_descriptors.loc[len(df_3d_descriptors)] = info_list
=======
def Mean_median_desc(filepath):
    data = CSV_Loader(filepath)
    num_conformers = 20
    
    # Create a list of descriptor names
    descriptor_names_3d = ['asphericities', 'eccentricities', 'inertialShapeFactors',
                           'NPR1s', 'NPR2s', 'PMI1s', 'PMI2s', 'PMI3s', 
                           'RadiusOfGyrations', 'SpherocityIndexs']
    descriptor_names_2d = [desc[0] for desc in Descriptors.descList]
    descriptor_names_all = descriptor_names_2d + descriptor_names_3d
    
    # Create a descriptor calculator for 2D
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names_2d)
    
    # Create list of the non-inhibitors and inhibitors
    non_inhibitor_values = []
    inhibitor_values =[]
    
    
    for index in range(len(data)):
        # Convert smiles to molecules
        mol = Chem.AddHs(Chem.MolFromSmiles(data["SMILES"][index])) #AddHs adds the hydrogen atoms
       
        # Non inhibitor
        if data["ALDH1_inhibition"][index] == 0:     
            # Make tuple with all 2d descriptors:
            non_inhib_desc_values= calc.CalcDescriptors(mol)
            
<<<<<<< HEAD
        # Make tuple with all 3d descriptors
        driedee_tuple = tuple(driedee_list)
        # Add both tuples to non inhibotor list
        inhibitor_value.append(inhib_desc_values + driedee_tuple)
        
        #Add 3D descriptor values to the 3D dataframe 
        info_list = [data["SMILES"][index],1]+driedee_list
        df_3d_descriptors.loc[len(df_3d_descriptors)] = info_list

# Calculate mean values        
mean_non_inhibitors = [sum(col) / len(col) for col in zip(*non_inhibitor_value)]
mean_inhibitors = [sum(col) / len(col) for col in zip(*inhibitor_value)]

# Perform t-test per descriptor for means
mean_ttest_results = []
significant_count = 0 
super_significant_count = 0
for i, descriptor in enumerate(descriptor_names_all):
    non_inhib_values = [desc[i] for desc in non_inhibitor_value]
    inhib_values = [desc[i] for desc in inhibitor_value]
    mean_ttest_result = stats.ttest_ind(non_inhib_values, inhib_values)
    mean_ttest_results.append(mean_ttest_result)
    if mean_ttest_result.pvalue < 0.05:
        significant_count += 1
    if mean_ttest_result.pvalue < 0.01:
        super_significant_count += 1

# Calculate median values
median_non_inhibitors = [np.median(col) for col in zip(*non_inhibitor_value)]
median_inhibitors = [np.median(col) for col in zip(*inhibitor_value)]

# Perform t-test per descriptor for medians
median_ttest_results = []
for i, descriptor in enumerate(descriptor_names_all):
    non_inhib_values = [desc[i] for desc in non_inhibitor_value]
    inhib_values = [desc[i] for desc in inhibitor_value]
    median_ttest_result = stats.ttest_ind(non_inhib_values, inhib_values)
    median_ttest_results.append(median_ttest_result)

# Create a DataFrame with the mean, median, and descriptor names
df = pd.DataFrame({'Descriptor': descriptor_names_all,
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
df_3d_descriptors.to_csv('3D_descriptor_values.csv',index=False)
=======
            # Calculate 3d descriptors:
            AllChem.EmbedMultipleConfs(mol, num_conformers)
            driedee_list = drieD_descriptors (num_conformers, mol)        
            
            # Make tuple with all 3d descriptors
            driedee_tuple = tuple(driedee_list)
            
            # Add both tuples to non inhibotor list
            non_inhibitor_values.append(non_inhib_desc_values + driedee_tuple)
        
        # Inhibitor
        else:
            # Make tuple with all 2d descriptors:
            inhib_desc_values = calc.CalcDescriptors(mol)
            
            # Calculate 3d descriptors:
            AllChem.EmbedMultipleConfs(mol, num_conformers)
            driedee_list = drieD_descriptors (num_conformers, mol)  
                
            # Make tuple with all 3d descriptors
            driedee_tuple = tuple(driedee_list)
            
            # Add both tuples to non inhibotor list
            inhibitor_values.append(inhib_desc_values + driedee_tuple)
    
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
    ttest_results= []
    significant_count = 0
    super_significant_count = 0
   
    for i, descriptor in enumerate(descriptor_names_all):
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
    df = pd.DataFrame({'Descriptor': descriptor_names_all,
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
    
    df.to_csv('means_median_stat_table_3D_2.csv', index=False)
    
    # make second table with only the descriptors with significant difference 
    df2 = df [df['p-value']<0.01]
    df2 = df2[['Descriptor', 'mean_non_inhibitors','mean_inhibitors', 'median_non_inhibitors', 'median_inhibitors', 'p-value',]]
    df2.to_csv('descriptors_low_Pvalue_3D.csv', index=False)

Mean_median_desc('tested_molecules-1.csv')
