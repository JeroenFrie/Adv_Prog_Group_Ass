Desc_lists = Mean_Median_Desc("tested_molecules_v3.csv")

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
