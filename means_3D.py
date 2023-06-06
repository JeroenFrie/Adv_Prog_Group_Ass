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
    non_inhibitor_value = []
    inhibitor_value =[]
    
    
    for index in range(len(data)):
        # Convert smiles to molecules
        mol = Chem.AddHs(Chem.MolFromSmiles(data["SMILES"][index])) #AddHs adds the hydrogen atoms
       
        # Non inhibitor
        if data["ALDH1_inhibition"][index] == 0:     
            # Make tuple with all 2d descriptors:
            non_inhib_desc_values= calc.CalcDescriptors(mol)
            
            # Calculate 3d descriptors:
            AllChem.EmbedMultipleConfs(mol, num_conformers)
            driedee_list = drieD_descriptors (num_conformers, mol)        
            
            # Make tuple with all 3d descriptors
            driedee_tuple = tuple(driedee_list)
            
            # Add both tuples to non inhibotor list
            non_inhibitor_value.append(non_inhib_desc_values + driedee_tuple)
        
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
            inhibitor_value.append(inhib_desc_values + driedee_tuple)
    
    # Calculate mean values        
    mean_non_inhibitors_list = [sum(col) / len(col) for col in zip(*non_inhibitor_value)]
    mean_inhibitors_list = [sum(col) / len(col) for col in zip(*inhibitor_value)]
    
    # Calculate median values
    median_non_inhibitors_list = [np.median(col) for col in zip(*non_inhibitor_value)]
    median_inhibitors_list = [np.median(col) for col in zip(*inhibitor_value)]
    
    # Perform t-test per descriptor for means
    ttest_results_mean = []
    significant_count_mean = 0 
    super_significant_count_mean = 0
    
    # Perform t-test per descriptor for medians
    ttest_results_median = []
    significant_count_median = 0
    super_significant_count_median = 0
   
    for i, descriptor in enumerate(descriptor_names_all):
        mean_non_inhibitors = [desc[i] for desc in non_inhibitor_value]
        mean_inhibitors = [desc[i] for desc in inhibitor_value]
        median_non_inhibitors = [desc[i] for desc in non_inhibitor_value]
        median_inhibitors = [desc[i] for desc in inhibitor_value]
        
        ttest_result_mean = stats.ttest_ind(mean_non_inhibitors, mean_inhibitors)
        ttest_result_median = stats.ttest_ind(median_non_inhibitors, median_inhibitors)
        
        ttest_results_mean.append(ttest_result_mean)
        ttest_results_median.append(ttest_result_median)
       
        if ttest_result_mean.pvalue < 0.05:
            significant_count_mean += 1
        if ttest_result_mean.pvalue < 0.01:
            super_significant_count_mean += 1
        if ttest_result_median.pvalue < 0.05:
            significant_count_median += 1
        if ttest_result_median.pvalue < 0.01:
            super_significant_count_median += 1         
    
    # Create a DataFrame with the mean, median, and descriptor names
    df = pd.DataFrame({'Descriptor': descriptor_names_all,
                       'mean_non_inhibitors': mean_non_inhibitors_list,
                       'mean_inhibitors': mean_inhibitors_list,
                       'median_non_inhibitors': median_non_inhibitors_list,
                       'median_inhibitors': median_inhibitors_list,
                       'T-Statistic Mean': [result.statistic for result in ttest_results_mean],
                       'p-value Mean': [result.pvalue for result in ttest_results_mean],
                       'T-Statistic Median': [result.statistic for result in ttest_results_median],
                       'p-value Median': [result.pvalue for result in ttest_results_median]})
    
    # Add significance columns
    df['Super_Significance Mean'] = df['p-value Mean'].apply(lambda p: 'woooooooooow' if p < 0.01 else 'nope, deze niet')
    df['p-value < 0.01 Count'] = super_significant_count_mean
    df['Significance Median'] = df['p-value Median'].apply(lambda p: 'wooooooow' if p < 0.01 else 'nope, deze niet')
    df['p-value < 0.01 Count'] = super_significant_count_median
    
    # Lijstjes voor Jeroen:
    super_significant_descriptors_mean_list = df[df['p-value Mean'] < 0.01]['Descriptor'].tolist()
    super_significant_descriptors_median_list = df[df['p-value Median'] < 0.01]['Descriptor'].tolist()
    
    df.to_csv('means_table_3D.csv', index=False)

Mean_median_desc('tested_molecules-1.csv')