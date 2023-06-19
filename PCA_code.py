import pandas as pd
import matplotlib.pyplot as plt
import sklearn.preprocessing as sp
from sklearn import decomposition, linear_model
import seaborn as sns
from CSV_Load import CSV_Loader

df_PCA = CSV_Loader("Descriptors_Vals_2D_3D.csv")

df_PCA_B = CSV_Loader("Descriptors_Vals_2D_3D.csv")
df_PCA_B = df_PCA_B.drop("SMILES", axis=1)
df_PCA_B = df_PCA_B.drop("ALDH1_inhibition", axis=1)

#scaler_type = sp.StandardScaler()
#scaler_type.fit(df_PCA_B)
#scaled_data = scaler_type.transform(df_PCA_B)

#standard_scaled = pd.DataFrame(scaled_data, columns=df_PCA_B.columns)
#sns.boxplot(minmax_scaled.iloc[0:len(minmax_scaled), 0:len(minmax_scaled.columns)])

pca = decomposition.PCA(n_components=10)
principal_components = pca.fit_transform(df_PCA_B)

component_names = [f"PC{i+1}" for i in range(principal_components.shape[1])]
principal_components = pd.DataFrame(principal_components, columns=component_names)

principal_components.insert(0, "SMILES", df_PCA["SMILES"])
principal_components.insert(0, "ALDH1_inhibition", df_PCA["ALDH1_inhibition"])

loadings = pd.DataFrame(pca.components_.T, columns=component_names, index=df_PCA_B.columns)

fig, axs = plt.subplots(nrows=1)
#fig.set_size_inches(12,18)

#sns.scatterplot(data=principal_components, x='PC1', y='PC2', hue="SMILES", palette="deep", ax=axs[0])
sns.scatterplot(data=principal_components, x='PC1', y='PC4', hue="ALDH1_inhibition", palette="deep")

fig.suptitle("Plots of PC scores with highlighted groups",fontweight="bold")
#for i in range(1):
 #   fig.get_axes()[i].set_ylabel("Scores on PC2")
  #  fig.get_axes()[i].set_xlabel("Scores on PC1")

fig, axs = plt.subplots(nrows=2)
fig.set_size_inches(20,35)

sns.scatterplot(data=loadings, x='PC1', y='PC4', ax=axs[1])
sns.scatterplot(data=principal_components, x='PC1', y='PC4', hue="ALDH1_inhibition", palette="deep", ax=axs[0])
print(loadings)
loadings['val'] = loadings.index

def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x']+.02, point['y'], point['val'])

label_point(loadings.PC1, loadings.PC4, loadings.val, plt.gca())


#plt.show()

