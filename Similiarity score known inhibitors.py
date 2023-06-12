import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import IPythonConsole
import pandas as pd
from rdkit.Chem import PandasTools, DataStructs, MACCSkeys
import matplotlib.pyplot as plt
import scipy.stats as stats

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
Smiles_ref_inhibitor = ["CN(C)C1=CC=C(C=C1)C=O", "O=C1C=C2NCCC2=CC1=O", "CC1=CC2=C(C=C1)NC(=O)C2=O", 
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


# # For exploratory data analysis purposes, plot boxplots for the 
# # similarity scores  the test inhibitors and each
# # reference ALDH1 inhibitor
# mol_data.boxplot(column=["Sim_inh0"], by='ALDH1_inhibition')
# mol_data.boxplot(column=["Sim_inh1"], by='ALDH1_inhibition')
# mol_data.boxplot(column=["Sim_inh2"], by='ALDH1_inhibition')
# mol_data.boxplot(column=["Sim_inh3"], by='ALDH1_inhibition')
# mol_data.boxplot(column=["Sim_inh4"], by='ALDH1_inhibition')
# plt.show()

def t_test(df: pd.DataFrame, numeric_var: str, categoric_var: str) -> list :
    """ Perform a t-test for two dataframes that have one numeric and one 
    categorical variable. The t-test will determine if the difference of the 
    numerical data of the categorical groups is significant. 
    
    parameters:
        df: dataframe with the two columns, one for the numeric variable and the
            other for the categorical varable
        numeric_var: columnname of the numeric variable
        categoric_var: columnname of the categoric variable
        
    returns:
        ttest_simularity: dataframe in which contains the name of the reference 
            inhibitor (molecular descriptor), the t-test value and the p-value
    """
    # Create dataframes for which only include either the inhibitors or the 
    # non-inhibitors
    df_inh = df[df[categoric_var]== 1]
    df_non_inh = df[df[categoric_var]== 0]
    # Test whether the groups of the inhibitors and non-inhibitors have equal
    # variances, since this is a parameter that is important for the t-test
    variance_equal_test = stats.levene(df_inh[numeric_var], 
                                       df_non_inh[numeric_var], center='median')
    # In case the p value is smaller than 0.05, the variance can assumed to be
    # equal
    if variance_equal_test[1] < 0.05:
        ttest_simularity = stats.ttest_ind(df_inh[numeric_var], 
                                           df_non_inh[numeric_var], 
                                           equal_var= True)
        return([ttest_simularity[0], ttest_simularity[1]])
    # In case the p value is greater or equal than 0.05, the variance is not
    # assumed to be equal
    elif variance_equal_test[1] >= 0.05:
        ttest_simularity = stats.ttest_ind(df_inh[numeric_var], 
                                           df_non_inh[numeric_var], 
                                           equal_var= False)
        return([ttest_simularity[0], ttest_simularity[1]])

def moods_median_test(df: pd.DataFrame, numeric_var: str,
                       categoric_var: str) -> list :
    """ Perform a t-test for two dataframes that have one numeric and one 
    categorical variable. The t-test will determine if the difference of the 
    numerical data of the categorical groups is significant. 
    
    parameters:
        df: dataframe with the two columns, one for the numeric variable and the
            other for the categorical varable
        numeric_var: columnname of the numeric variable
        categoric_var: columnname of the categoric variable
        
    returns:
        statistic: the mood statistic for the median
        p_value: p value for the mood statistic for the median
    """
    # Create dataframes for which only include either the inhibitors or the 
    # non-inhibitors
    df_inh = df[df[categoric_var]== 1]
    df_non_inh = df[df[categoric_var]== 0]
    # Conduct the moods test for calculating if the difference in significancy
    statistic, p_value, median, table= stats.median_test(df_inh[numeric_var],
                                                          df_non_inh[numeric_var])
    return([statistic, p_value])


# Create an empty dataframe in which the statistics of the t-test and mood test
# will be stored later
ttest_moods_statistics= pd.DataFrame(columns= ["molecular descriptor", "SMILES",
                                         "t-test score", "p-value t-test",
                                           "significant t-test", 
                                           "moods score", 
                                           "p-value mood",
                                           "significant mood"])

for i in range(len(Smiles_ref_inhibitor)):
    # Select only the columns needed for calculating the statistics for 
    # the specific reference inhibitor
    mol_data_select = mol_data[["ALDH1_inhibition", "Sim_inh"+str(i)]].copy()
    # Change the type of the values in the column "ALDH1_inhibition" to strings
    mol_data_select["ALDH1_inhibition"].apply(str)
    # Perform the t-test and mood test for all reference inhibitors
    results_t_test= t_test(mol_data_select, "Sim_inh"+str(i),
                            "ALDH1_inhibition")
    results_moods_test = moods_median_test(mol_data_select, "Sim_inh"+str(i),
                            "ALDH1_inhibition")
    # Add the statistics of the tests to the dataframe that was created above 
    # the for-loop
    ttest_moods_statistics.loc[str(i)] = ["Sim_inh"+str(i), 
                                          Smiles_ref_inhibitor[i],
                                     results_t_test[0], results_t_test[1],
                                       results_t_test[1] < 0.01, 
                                       results_moods_test[0], 
                                       results_moods_test[1], 
                                       results_moods_test[1] < 0.01]
print(ttest_moods_statistics.head())
# Save dataframe to csv file
ttest_moods_statistics.to_csv("ttest_moods_statistics_similarity.csv")



# Only Sim_inh3 had insignificant differences in the ttest and moods test for
# the median, for the other reference inhibitors a mean, median dataframe is 
# created
Smiles_ref_inhibitor_insig_removed = ["CN(C)C1=CC=C(C=C1)C=O", 
                                      "O=C1C=C2NCCC2=CC1=O", 
                                      "CC1=CC2=C(C=C1)NC(=O)C2=O",
                                      "CC(C)CCN1C(=NC2=C1C(=O)N(C(=O)N2C)C)CN3"+
                                      "CCN(CC3)C(=O)C4CC4"]

# Create separate dataframes for the inhibitors and non-inhibitors data
mol_inh = mol_data[mol_data["ALDH1_inhibition"] == 1]
mol_non_inh = mol_data[mol_data["ALDH1_inhibition"] == 0]
# Create a dataframe in which the means and medians of the similarity score
# of reference inhibitors are stored, divided for the inhibitors and 
# non-inhibitors
mean_median_similarity= pd.DataFrame(columns= 
                                     ["SMILES", "mean_inh", 
                                      "mean_non_inh", "median_inh", 
                                      "median_non_inh"])
for i in range(len(Smiles_ref_inhibitor_insig_removed)):
        mean_median_similarity.loc[str(i)] = [Smiles_ref_inhibitor_insig_removed[i], 
                                              mol_inh["Sim_inh"+str(i)].mean(), 
                                               mol_non_inh["Sim_inh"+str(i)].mean(),
                                               mol_inh["Sim_inh"+str(i)].median(),
                                               mol_non_inh["Sim_inh"+str(i)].median()]
print(mean_median_similarity)
# Save dataframe to csv file
mean_median_similarity.to_csv("mean_median_similarity.csv")