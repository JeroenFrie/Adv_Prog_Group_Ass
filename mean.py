from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
import pandas as pd
from CSV_Load import CSV_Loader
from statistics import mean

data = CSV_Loader("tested_molecules-1.csv")

non_inhibitor_mass = []
inhibitor_mass =[]

for index in range(len(data)):
    mol = Chem.MolFromSmiles(data["SMILES"[index]])
    if data["ALDH1_inhibition"][index] == 0:
        non_inhibitor_mass.append(Discriptors.MolWt(mol))
    else:
        inhibitor_mass.append(Discriptors.MolWt(mol))

mean_mass_non_inhibitors = mean(non_inhibitor_mass)
mean_mass_inhibitors = mean(inhibitor_mass)