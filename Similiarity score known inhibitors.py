import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import IPythonConsole
import pandas as pd
from rdkit.Chem import PandasTools, DataStructs, MACCSkeys
import matplotlib.pyplot as plt

# Sim_inh0: "CN(C)C1=CC=C(C=C1)C=O": p-Dimethylaminobenzaldehyde hydrobromide    
# Source for inhibitor: https://pubmed.ncbi.nlm.nih.gov/25512087/

# Sim_inh2:"O=C1C=C2NCCC2=CC1=O": 5,6-Indolinedione
# Source for inhibitor: https://pubmed.ncbi.nlm.nih.gov/25512087/ 
# (says indolinedione analogs are good ALDHA1 inhibitors)

# Sim_inh3:"CC1=CC2=C(C=C1)NC(=O)C2=O": 5-methyl-1H-indole-2,3-dione             
# Source for inhibitor: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3954746/

# Sim_inh4:"C=CC(=O)C1=CC2=CC=CC=C2C=C1": 1-(2-Naphthalenyl)-2-propen-1-one  
# Source for inhibitor: (SID:26753168)
# https://pubchem.ncbi.nlm.nih.gov/bioassay/1030#section=Data-Table 
    
# Sim_inh5:"CC(C)CCN1C(=NC2=C1C(=O)N(C(=O)N2C)C)CN3CCN(CC3)C(=O)C4CC4": Nct-501
# Source for inhibitor: https://europepmc.org/article/pmc/5185321

# List which contains the SMILES of the reference inhibitors 
Smiles_ref_inhibitor = ["CN(C)C1=CC=C(C=C1)C=O", 
                         "O=C1C=C2NCCC2=CC1=O", "CC1=CC2=C(C=C1)NC(=O)C2=O", 
                         "C=CC(=O)C1=CC2=CC=CC=C2C=C1", 
                         "CC(C)CCN1C(=NC2=C1C(=O)N(C(=O)N2C)C)CN3CCN(CC3)C(=O)"
                         +"C4CC4"]

def smiles_to_fingerprint(smiles: str) -> rdkit: 
    """ Go from the SMILEs of a molecule to the molecules fingerprint, this 
     describes a 3D molecule in a bit structure
      
    Parameters:
        smiles: notation to describe a molecule in a line form
        (Simplified Molecular Input Line Entry System)
    Returns:
        fp: rdkit.datastruct in which the bit information of a molecular
            fingerprint is stored
     """
    try:
        # Get the molecular object from the smiles string
        mol = Chem.MolFromSmiles(smiles)
        # Get the molecular fingerprint from the molecular object
        fp = MACCSkeys.GenMACCSKeys(mol)

        ## Other method which can be used to get the fingerprint (do not see
        # much difference by eye)
        # radius input tells how many bonds from the central atom will be 
        # considered

        # fp = (AllChem.GetMorganFingerprintAsBitVect(mol,radius=3, nBits=1024))
        return fp
    except:
        return None

# Call function smiles_to_fingerprint for the list of reference ALDH1 inhibitors
reference_inhibitor_fps = [smiles_to_fingerprint(smiles= smiles_ref) 
                           for smiles_ref in Smiles_ref_inhibitor]

# Read the test_set of molecules that either do or do not inhibit ALDH1 and save 
# the SMILES of these molecules to a list
path = "tested_molecules_v3.csv"
mol_data = pd.read_csv(path)
SMILES_df = mol_data["SMILES"].tolist()
# Call function smiles_to_fingerprint for the list of test inhibitors
test_inhibitors_fps = [smiles_to_fingerprint(smiles= smiles_test) 
                       for smiles_test in SMILES_df]


# Loop over the reference inhibitors
for nr_ref in range(len(reference_inhibitor_fps)):
    # Create a empty list to store the tanimoto scores
    tanimoto = []
    # Loop over the test inhibitors and determine the tanimoto scores
    # (similarities) for each reference inhibitors to the test inhibitors.
    for inh in range(len(test_inhibitors_fps)):
        tanimoto_scores = DataStructs.TanimotoSimilarity(
            reference_inhibitor_fps[nr_ref],
            test_inhibitors_fps[inh])
        tanimoto.append(tanimoto_scores)
    # Add the tanimoto list to a new column in the mol_data data frame
    mol_data["Sim_inh"+str(nr_ref)] = tanimoto

# For exploratory data analysis purposes, plot boxplots for the 
# similarity scores  the test inhibitors and each
# reference ALDH1 inhibitor
mol_data.boxplot(column=["Sim_inh0"], by='ALDH1_inhibition')
mol_data.boxplot(column=["Sim_inh1"], by='ALDH1_inhibition')
mol_data.boxplot(column=["Sim_inh2"], by='ALDH1_inhibition')
mol_data.boxplot(column=["Sim_inh3"], by='ALDH1_inhibition')
mol_data.boxplot(column=["Sim_inh4"], by='ALDH1_inhibition')
plt.show()

