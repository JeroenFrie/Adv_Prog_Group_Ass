import pandas as pd
from CSV_Load import CSV_Loader

df = CSV_Loader("tested_molecules-1.csv")
pd.to_numeric(df["ALDH1_inhibition"])
df.to_csv("tested_molecules-1.csv", index=False)