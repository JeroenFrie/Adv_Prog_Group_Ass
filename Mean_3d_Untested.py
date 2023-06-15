from rdkit import Chem
from rdkit.Chem import Descriptors3D, AllChem
from CSV_Load import CSV_Loader
import pandas as pd

data = CSV_Loader("untested_molecules.csv")

num_conformers = 20

# Create a list of descriptor names
descriptor_names_3d = ['PMI1s']

# Descriptor calculator for 3D
def drieD_descriptors(num_conformers, mol, smiles):
    drieD_list = []
    # List for every descriptor
    PMI1s = []

    # Iteratate for every conformation
    for conf_id in range(num_conformers):
        conf = mol.GetConformer(conf_id)

        # Calculate value of every descriptor
        PMI1 = Descriptors3D.PMI1(mol, confId=conf_id)

        # Add descriptor value to list
        PMI1s.append(PMI1)

    # Calculate the mean per descriptor for the total number of conformations
    drieD_descriptors_list = [PMI1s]
    drieD_list.append(smiles)
    for index in range(len(drieD_descriptors_list)):
        mean = sum(drieD_descriptors_list[index]) / num_conformers
        drieD_list.append(mean)
    return drieD_list

info = ({
    'SMILES': [],
    'PMI1': []
})
df_3d_descriptors = pd.DataFrame(info)

for index in range(len(data)):
    mol = Chem.AddHs(Chem.MolFromSmiles(data["SMILES"][index]))  # AddHs adds the hydrogen atoms
    smile = data["SMILES"][index]
    # Calculate 3d descriptors:
    AllChem.EmbedMultipleConfs(mol, num_conformers)
    driedee_list = drieD_descriptors(num_conformers, mol, smile)
    print(driedee_list)
    df_3d_descriptors.loc[len(df_3d_descriptors)] = driedee_list
    print(df_3d_descriptors)


print(df_3d_descriptors)
#df_3d_descriptors.to_csv('3D_descriptor_values.csv', index=False)
