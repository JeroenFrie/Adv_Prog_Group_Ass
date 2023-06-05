from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
from mean import Mean_Median_Desc
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.preprocessing as sp
from sklearn import decomposition, linear_model
import seaborn as sns
from CSV_Load import CSV_Loader

Molecule_DF = CSV_Loader("tested_molecules-1.csv")
DrieD_Mol_DF = CSV_Loader("3D_descriptor_values.csv")

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

df_PCA = pd.merge(Molecule_DF, DrieD_Mol_DF)

df_PCA_B = df_PCA
df_PCA_B = df_PCA_B.drop("SMILES", axis=1)
df_PCA_B = df_PCA_B.drop("ALDH1_inhibition", axis=1)

scaler_type = sp.MinMaxScaler()
scaler_type.fit(df_PCA_B)
scaled_data = scaler_type.transform(df_PCA_B)

minmax_scaled = pd.DataFrame(scaled_data, columns=df_PCA_B.columns)
#sns.boxplot(minmax_scaled.iloc[0:len(minmax_scaled), 0:len(minmax_scaled.columns)])

pca = decomposition.PCA(n_components=25)
principal_components = pca.fit_transform(minmax_scaled)

component_names = [f"PC{i+1}" for i in range(principal_components.shape[1])]
principal_components = pd.DataFrame(principal_components, columns=component_names)

principal_components.insert(0, "SMILES", df_PCA["SMILES"])
principal_components.insert(0, "ALDH1_inhibition", df_PCA["ALDH1_inhibition"])

loadings = pd.DataFrame(pca.components_.T, columns=component_names, index=minmax_scaled.columns)

fig, axs = plt.subplots(nrows=1)
fig.set_size_inches(12,18)

#sns.scatterplot(data=principal_components, x='PC1', y='PC2', hue="SMILES", palette="deep", ax=axs[0])
sns.scatterplot(data=principal_components, x='PC1', y='PC2', hue="ALDH1_inhibition", palette="deep")

fig.suptitle("Plots of PC scores with highlighted groups",fontweight="bold")
for i in range(1):
    fig.get_axes()[i].set_ylabel("Scores on PC2")
    fig.get_axes()[i].set_xlabel("Scores on PC1")

fig, axs = plt.subplots(nrows=2)
fig.set_size_inches(20,35)

sns.scatterplot(data=loadings, x='PC1', y='PC2', ax=axs[1])
sns.scatterplot(data=principal_components, x='PC1', y='PC2', hue="ALDH1_inhibition", palette="deep", ax=axs[0])

loadings['val'] = loadings.index

def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x']+.02, point['y'], point['val'])

#label_point(loadings.PC1, loadings.PC2, loadings.val, plt.gca())


plt.show()

