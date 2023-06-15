from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
import pandas as pd
from CSV_Load import CSV_Loader
import sklearn.preprocessing as sp
from sklearn import decomposition, linear_model

Molecules_DF = CSV_Loader("untested_molecules.csv")
Descriptor_DF = CSV_Loader("Descriptors_Vals_2D_3D.csv")

Column_list = list(Descriptor_DF.columns)

for i in range(2):
    Column_list.pop(0)

ToD_Desc = Column_list[0:Column_list.index("PMI1")]

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

