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
from statistics import mean, median

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
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names_2d) #ChatGPT
    
    # Create list of the non-inhibitors and inhibitors
    non_inhibitor_values = []
    inhibitor_values =[]
    
    
    for index in range(len(data)):
        # Convert smiles to molecules
        mol = Chem.AddHs(Chem.MolFromSmiles(data["SMILES"][index])) #AddHs adds the hydrogen atoms #ChatGPT
       
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
            non_inhibitor_values.append(non_inhib_desc_values + driedee_tuple)
        
        # Inhibitor
        else:
            # Make tuple with all 2D descriptors:
            inhib_desc_values = calc.CalcDescriptors(mol)
            
            # Calculate 3D descriptors:
            AllChem.EmbedMultipleConfs(mol, num_conformers)
            driedee_list = drieD_descriptors (num_conformers, mol)  
                
            # Make tuple with all 3d descriptors
            driedee_tuple = tuple(driedee_list)
            
            # Add both tuples to non inhibotor list
            inhibitor_values.append(inhib_desc_values + driedee_tuple)   
   
    # Calculate the mean
    mean_non_inhibitors = mean(non_inhibitor_values)
    mean_inhibitors = mean(inhibitor_values)
   
    # Calculate median 
    median_non_inhibitors = median(non_inhibitor_values)
    median_inhibitors = median(inhibitor_values)
    
    # Perform t-test per descriptor
    ttest_results= []
    significant_count = 0
    super_significant_count = 0
   
    for i, descriptor in enumerate(descriptor_names_all):
        non_inhibitor = non_inhibitor_values [:,i]
        inhibitor = inhibitor_values [:,i]
        
        ttest_result = stats.ttest_ind(non_inhibitor, inhibitor)
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
    df['Super_Significance'] = df['p-value'].apply(
        lambda p: 'Significant difference' if p < 0.01 else 'nope')
    df['p-value < 0.01 Count'] = super_significant_count
    
    df.to_csv('statistics_table_2D_3D.csv', index=False)
    
    # Make second table with only the descriptors with significant difference 
    df2 = df [df['p-value']<0.01]
    df2 = df2[['Descriptor', 'mean_non_inhibitors','mean_inhibitors', 'median_non_inhibitors', 'median_inhibitors', 'p-value',]]
    df2.to_csv('descriptors_low_pvalue_2D_3D.csv', index=False)

Mean_median_desc('tested_molecules-1.csv')
