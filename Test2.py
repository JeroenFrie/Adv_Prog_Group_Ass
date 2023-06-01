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
num_conformers = 2

# Create a list of descriptor names
descriptor_names_3d = ['asphericities', 'eccentricities', 'inertialShapeFactors',
                       'NPR1s', 'NPR2s', 'PMI1s', 'PMI2s', 'PMI3s', 
                       'RadiusOfGyrations', 'SpherocityIndexs']
descriptor_names = [desc[0] for desc in Descriptors.descList]
descriptor_names_all = descriptor_names + descriptor_names_3d

# Create a descriptor calculator 
calc = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)

non_inhibitor_value = []
inhibitor_value =[]


def drieD_descriptors (num_conformers, mol):
    driedee_list = []
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
        
    # Calculate the mean per discriptor for the total number of conformations
    drieD_descriptors_list = [asphericities, eccentricities, 
                              inertialShapeFactors, NPR1s, NPR2s, 
                              PMI1s, PMI2s, PMI3s, RadiusOfGyrations,
                              SpherocityIndexs ]
    for descriptor in drieD_descriptors_list:
        mean = sum(descriptor) / num_conformers
        driedee_list.append(mean)  
    return driedee_list



for index in range(len(data)):
    mol = Chem.AddHs(Chem.MolFromSmiles(data["SMILES"][index]))
    if data["ALDH1_inhibition"][index] == 0:
        # Make tuple with all 2d descriptors:
        non_inhib_desc_values= calc.CalcDescriptors(mol)
        
        # Calculate 3d descriptors:
        AllChem.EmbedMultipleConfs(mol, num_conformers)
        driedee_list = drieD_descriptors (num_conformers, mol)        
        
        driedee_tuple = tuple(driedee_list)
        non_inhibitor_value.append(non_inhib_desc_values + driedee_tuple)
    else:
        inhib_desc_values = calc.CalcDescriptors(mol)
        # Calculate 3d descriptors:
        AllChem.EmbedMultipleConfs(mol, num_conformers)
        driedee_list = drieD_descriptors (num_conformers, mol)  
            
        driedee_tuple = tuple(driedee_list)
        inhibitor_value.append(inhib_desc_values + driedee_tuple)
        
mean_non_inhibitors = [sum(col) / len(col) for col in zip(*non_inhibitor_value)]
mean_inhibitors = [sum(col) / len(col) for col in zip(*inhibitor_value)]



