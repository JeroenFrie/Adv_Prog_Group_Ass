from rdkit import Chem
from rdkit.Chem import Descriptors3D, AllChem
from CSV_Load import CSV_Loader
import pandas as pd

# Loading in the csv containing all untested molecules
data = CSV_Loader("untested_molecules.csv")

# Stating amount of conformers per molecule
num_conformers = 20

# Create a list of descriptor names
descriptor_names_3d = ['PMI1s']

# Descriptor calculator for 3D
def drieD_descriptors(num_conformers, mol, smiles):
    """Function to calculate PMI1 for a given molecule
               ----------
               num_conformers : int
                   variable stating the amount of conformers
               smiles : str
                    string containing SMILE of molecule
               mol : str
                    string containing molecule

               Returns
               -------
               drieD_list : list
                   list containing molecule and value of PMI1 of that molecule
           """
    drieD_list = []
    # List for every descriptor
    PMI1s = []

    if mol.GetNumConformers() < num_conformers:
        PMI1s.append(None)
    else:
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
        if drieD_descriptors_list[index][0] != None:
            mean = sum(drieD_descriptors_list[index]) / num_conformers
            drieD_list.append(mean)
        else:
            drieD_list.append(drieD_descriptors_list[index])

    return drieD_list

info = ({
    'SMILES': [],
    'PMI1': []
})
df_3d_descriptors = pd.DataFrame(info)

# For loop calculating the 3D descriptors for each molecule loaded in before
for index in range(len(data)):
    mol = Chem.AddHs(Chem.MolFromSmiles(data["SMILES"][index]))
    smile = data["SMILES"][index]
    AllChem.EmbedMultipleConfs(mol, num_conformers)
    driedee_list = drieD_descriptors(num_conformers, mol, smile)
    df_3d_descriptors.loc[len(df_3d_descriptors)] = driedee_list

# Export dataframe to csv file
df_3d_descriptors.to_csv('3D_descriptor_values_un.csv', index=False)
