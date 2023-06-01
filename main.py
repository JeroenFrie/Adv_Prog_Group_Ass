# import rdkit.Chem.rdMolDescriptors
from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
from mean import Mean_Median_Desc
import pandas as pd
from CSV_Load import CSV_Loader

Molecule_DF = CSV_Loader("tested_molecules-1.csv")

def Inter_Corr(Dataframe):
    desc_high_corr = []
    column_name_list = []
    for column_val in range(len(Dataframe.columns)):
        column_name_list.append(Dataframe.columns[column_val])
        for row in range(len(Dataframe)):
            corr_val = Dataframe.iloc[row, column_val]
            temp_list = []
            if corr_val >= 0.9 and corr_val != 1 and Dataframe.iloc[[row]].index[0] not in \
                    column_name_list or corr_val <= -0.9 and corr_val != 1 and Dataframe.iloc[[row]].index[0] not in \
                    column_name_list:
                temp_list.append(Dataframe.columns[column_val])
                temp_list.append(Dataframe.iloc[[row]].index[0])
                desc_high_corr.append(temp_list)

    refine_mean_corr = []
    for desc_in in range(len(desc_high_corr)):
        if desc_high_corr[desc_in][0] not in refine_mean_corr:
            refine_mean_corr.append(desc_high_corr[desc_in][0])
    return refine_mean_corr


index_list = []
for index in range(len(Molecule_DF)):
    if Molecule_DF["ALDH1_inhibition"][index] == 0:
        index_list.append(index)
Molecule_DF_In = Molecule_DF.drop(index=index_list)

Mol_list = []
for row in range(len(Molecule_DF)):
    Mol_list.append(Chem.MolFromSmiles(Molecule_DF["SMILES"][row]))

desc_list = [n[0] for n in Descriptors._descList]
short_desc = [i for i in desc_list if not i.startswith("fr_")]

for desc in short_desc:
    temp_list = []
    calc = MoleculeDescriptors.MolecularDescriptorCalculator([desc])
    for MOL in Mol_list:
        val_desc = calc.CalcDescriptors(MOL)[0]
        temp_list.append(val_desc)
    Molecule_DF[desc] = temp_list

Re_Molecule_DF = Molecule_DF
Corr_Mol_Des = Molecule_DF.corr("spearman", numeric_only=True)

High_Corr = []
for row in range(len(Corr_Mol_Des)):
    temp_list = []
    Corr_val = Corr_Mol_Des["ALDH1_inhibition"][row]
    if Corr_val >= 0.1 and Corr_val != 1 or Corr_val <= -0.1 and Corr_val != 1:
        temp_list.append(row)
        temp_list.append(Corr_val)
        High_Corr.append(temp_list)

index_corr_list = []
for corr_mol in range(len(High_Corr)):
    row = High_Corr[corr_mol][0]
    index_corr_list.append(Corr_Mol_Des.iloc[[row]].index[0])

for name in short_desc:
    if name not in index_corr_list:
        Re_Molecule_DF = Re_Molecule_DF.drop(name, axis=1)

desc_df = Re_Molecule_DF.iloc[:, 2:len(Re_Molecule_DF.columns)]
desc_corr_df = desc_df.corr("spearman", numeric_only=True)

refine_desc_corr = Inter_Corr(desc_corr_df)

for name in refine_desc_corr:
    Re_Molecule_DF = Re_Molecule_DF.drop(name, axis=1)

Desc_lists = Mean_Median_Desc("tested_molecules-1.csv")
Desc_Mean_list = Desc_lists[0]
Desc_Median_list = Desc_lists[1]

Mean_Median_DF = Molecule_DF
desc_short_list = []
for name in short_desc:
    if name not in Desc_Mean_list:
        Mean_Median_DF = Mean_Median_DF.drop(name, axis=1)
    else:
        desc_short_list.append(name)


desc_mean_df = Mean_Median_DF.iloc[:, 2:len(Mean_Median_DF.columns)]
Mean_Median_Corr = desc_mean_df.corr("spearman", numeric_only=True)

refine_mean_corr = Inter_Corr(Mean_Median_Corr)

for name in refine_mean_corr:
    desc_short_list.remove(name)

for name in desc_short_list:
    if name not in Re_Molecule_DF.columns:
        Re_Molecule_DF[name] = Mean_Median_DF[name]


print(Re_Molecule_DF)