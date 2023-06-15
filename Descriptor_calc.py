from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
import pandas as pd
from CSV_Load import CSV_Loader
import sklearn.preprocessing as sp
from sklearn import decomposition, linear_model

Molecules_DF = CSV_Loader("untested_molecules.csv")
Descriptor_DF = CSV_Loader("Descriptors_Vals_2D_3D.csv")
Similarity_Vals_DF = CSV_Loader("Similarity_val_un.csv")
DrieD_DF = CSV_Loader("3D_descriptor_values_un.csv")

Column_list = list(Descriptor_DF.columns)

for i in range(2):
    Column_list.pop(0)

ToD_Desc = Column_list[0:Column_list.index("PMI1")]
Sim_Desc = Column_list[Column_list.index("PMI1")+1:len(Column_list)]

Mol_list = []
for row in range(len(Molecules_DF)):
    Mol_list.append(Chem.MolFromSmiles(Molecules_DF["SMILES"][row]))

for desc in ToD_Desc:
    temp_list = []
    calc = MoleculeDescriptors.MolecularDescriptorCalculator([desc])
    for MOL in Mol_list:
        val_desc = calc.CalcDescriptors(MOL)[0]
        temp_list.append(val_desc)
    Molecules_DF[desc] = temp_list

Molecules_DF = pd.concat([Molecules_DF,DrieD_DF["PMI1"]], axis=1)

for Sim in Similarity_Vals_DF.columns:
    if Sim not in Sim_Desc:
        Similarity_Vals_DF = Similarity_Vals_DF.drop(Sim, axis=1)

all_col = Similarity_Vals_DF.iloc[:, :]
Molecules_DF = pd.concat([Molecules_DF, all_col], axis=1)


scaler_type = sp.StandardScaler()
Smiles_df = Molecules_DF["SMILES"]
Molecules_DF = Molecules_DF.drop("SMILES", axis=1)

scaler_type.fit(Molecules_DF)
scaled_data = scaler_type.transform(Molecules_DF)

standard_scaled = pd.DataFrame(scaled_data, columns=Molecules_DF.columns)

standard_scaled.insert(0, "SMILES", Smiles_df["SMILES"])

standard_scaled.set_index("SMILES")

standard_scaled.to_csv("Unknown_Mol_Desc.csv", index=False)

