#import rdkit.Chem.rdMolDescriptors
from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
import pandas as pd
from CSV_Load import CSV_Loader

Molecule_DF = CSV_Loader("tested_molecules-1.csv")

index_list = []
for index in range(len(Molecule_DF)):
    if Molecule_DF["ALDH1_inhibition"][index] == 0:
        index_list.append(index)
Molecule_DF_In = Molecule_DF.drop(index=index_list)

Mol_list = []
for row in range(len(Molecule_DF)):
    Mol_list.append(Chem.MolFromSmiles(Molecule_DF["SMILES"][row]))

#Rot_bonds = []
#for MOL in range(len(Mol_list)):
 #   Rot_bonds.append(rdkit.Chem.rdMolDescriptors.CalcNumRotatableBonds((Mol_list[MOL])))
#Molecule_DF["Rot_Bonds"] = Rot_bonds

desc_list = [n[0] for n in Descriptors._descList]
short_desc = [i for i in desc_list if not i.startswith("fr_")]

for desc in short_desc:
    temp_list = []
    calc = MoleculeDescriptors.MolecularDescriptorCalculator([desc])
    for MOL in Mol_list:
        val_desc = calc.CalcDescriptors(MOL)[0]
        temp_list.append(val_desc)
    Molecule_DF[desc] = temp_list

Corr_Mol_Des = Molecule_DF.corr("spearman", numeric_only=True)

High_Corr = []
for row in range(len(Corr_Mol_Des)):
    temp_list = []
    Corr_val = Corr_Mol_Des["ALDH1_inhibition"][row]
    if Corr_val >= 0.1 and Corr_val != 1 or Corr_val <= -0.1 and Corr_val != 1:
        temp_list.append(row)
        temp_list.append(Corr_val)
        High_Corr.append(temp_list)


index_corr_list=[]
for corr_mol in range(len(High_Corr)):
    row = High_Corr[corr_mol][0]
    index_corr_list.append(Corr_Mol_Des.iloc[[row]].index[0])

for name in short_desc:
    if name not in index_corr_list:
        Molecule_DF = Molecule_DF.drop(name, axis=1)

desc_df = Molecule_DF.iloc[:, 2:len(Molecule_DF.columns)]
desc_corr_df = desc_df.corr("spearman", numeric_only=True)

desc_high_corr = []
column_name_list = []
for column_val in range(len(desc_corr_df.columns)):
    column_name_list.append(desc_corr_df.columns[column_val])
    for row in range(len(desc_corr_df)):
        Corr_val = desc_corr_df.iloc[row, column_val]
        temp_list = []
        print(Corr_val)
        if Corr_val >= 0.9 and Corr_val != 1 and desc_corr_df.iloc[[row]].index[0] not in\
                column_name_list or Corr_val <= -0.9 and Corr_val != 1 and desc_corr_df.iloc[[row]].index[0] not in\
                column_name_list:
            temp_list.append(desc_corr_df.columns[column_val])
            temp_list.append(desc_corr_df.iloc[[row]].index[0])
            desc_high_corr.append(temp_list)

refine_desc_corr = []
for desc_in in range(len(desc_high_corr)):
    if desc_high_corr[desc_in][0] not in refine_desc_corr:
        refine_desc_corr.append(desc_high_corr[desc_in][0])

print(desc_high_corr)
print(len(desc_high_corr))

print(refine_desc_corr)
print(len(refine_desc_corr))

for name in index_corr_list:
    if name in refine_desc_corr:
        Molecule_DF = Molecule_DF.drop(name, axis=1)

print(Molecule_DF)