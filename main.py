from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
import pandas as pd
from CSV_Load import CSV_Loader
import sklearn.preprocessing as sp
from sklearn import decomposition, linear_model

Molecule_DF = CSV_Loader("tested_molecules_v3.csv")
Molecule_DF = Molecule_DF.drop(index=909)
New_Index = list(range(len(Molecule_DF)))

Molecule_DF["New_Index"] = New_Index

Molecule_DF = Molecule_DF.set_index("New_Index")

DrieD_Mol_DF = CSV_Loader("3D_descriptor_values.csv")


def Corr_Calc(Dataframe):
    High_Corr = []
    for row in range(len(Dataframe)):
        temp_list = []
        Corr_val = Dataframe["ALDH1_inhibition"][row]
        if Corr_val >= 0.3 and Corr_val != 1 or Corr_val <= -0.3 and Corr_val != 1:
            temp_list.append(row)
            temp_list.append(Corr_val)
            High_Corr.append(temp_list)

    index_corr_list = []
    for corr_mol in range(len(High_Corr)):
        row = High_Corr[corr_mol][0]
        index_corr_list.append(Dataframe.iloc[[row]].index[0])
    return index_corr_list

def Inter_Corr(Dataframe):
    desc_high_corr = []
    column_name_list = []
    if len(Dataframe.columns) > 1:
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

Mol_list = []
for row in range(len(Molecule_DF)):
    Mol_list.append(Chem.MolFromSmiles(Molecule_DF["SMILES"][row]))

desc_list = [n[0] for n in Descriptors._descList]
#short_desc = [i for i in desc_list if not i.startswith("fr_")]

for desc in desc_list:
    temp_list = []
    calc = MoleculeDescriptors.MolecularDescriptorCalculator([desc])
    for MOL in Mol_list:
        val_desc = calc.CalcDescriptors(MOL)[0]
        temp_list.append(val_desc)
    Molecule_DF[desc] = temp_list

Re_Molecule_DF = Molecule_DF

Corr_Mol_Des = Molecule_DF.corr("spearman", numeric_only=True)

index_corr_list = Corr_Calc(Corr_Mol_Des)

for name in desc_list:
    if name not in index_corr_list:
        Re_Molecule_DF = Re_Molecule_DF.drop(name, axis=1)

desc_df = Re_Molecule_DF.iloc[:, 2:len(Re_Molecule_DF.columns)]
desc_corr_df = desc_df.corr("spearman", numeric_only=True)

refine_desc_corr = Inter_Corr(desc_corr_df)

for name in refine_desc_corr:
    Re_Molecule_DF = Re_Molecule_DF.drop(name, axis=1)

Corr_Drie_D = DrieD_Mol_DF.corr("spearman",numeric_only=True)

DrieD_corr_index = Corr_Calc(Corr_Drie_D)

for name in DrieD_Mol_DF.columns:
    if name not in DrieD_corr_index and name != "SMILES" and name != "ALDH1_inhibition":
        DrieD_Mol_DF = DrieD_Mol_DF.drop(name, axis=1)

Column_Values = DrieD_Mol_DF.columns
Pre_Inter_Corr = DrieD_Mol_DF.iloc[:, 2:len(DrieD_Mol_DF.columns)]

Inter_Corr_Val = Inter_Corr(Pre_Inter_Corr)

for name in Column_Values:
    if name in Inter_Corr_Val and name != "SMILES" and name != "ALDH1_inhibition":
        DrieD_Mol_DF = DrieD_Mol_DF.drop(name, axis=1)

for name in DrieD_Mol_DF.columns:
    if name not in Re_Molecule_DF.columns:
        Re_Molecule_DF[name] = DrieD_Mol_DF[name]

similarity_dataframe = CSV_Loader("similarity_data.csv")
similarity_dataframe = similarity_dataframe.drop(index=909)
New_Index = list(range(len(similarity_dataframe)))

similarity_dataframe["New_Index"] = New_Index

similarity_dataframe = similarity_dataframe.set_index("New_Index")
for i in range(3, 7):
    Re_Molecule_DF = pd.concat([Re_Molecule_DF, similarity_dataframe.iloc[:, i]], axis=1)

scaler_type = sp.StandardScaler()
scale_data_df = Re_Molecule_DF
scale_data_df = scale_data_df.drop("SMILES", axis=1)
scale_data_df = scale_data_df.drop("ALDH1_inhibition", axis=1)
scaler_type.fit(scale_data_df)
scaled_data = scaler_type.transform(scale_data_df)

standard_scaled = pd.DataFrame(scaled_data, columns=scale_data_df.columns)

standard_scaled.insert(0, "ALDH1_inhibition", Re_Molecule_DF["ALDH1_inhibition"])
standard_scaled.insert(0, "SMILES", Re_Molecule_DF["SMILES"])

standard_scaled.set_index("SMILES")

print(standard_scaled)
standard_scaled.to_csv("Descriptors_Vals_2D_3D.csv", index=False)