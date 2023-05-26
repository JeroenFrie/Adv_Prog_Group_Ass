import pandas as pd

def CSV_Loader(filepath):
    csv_df = pd.read_csv(filepath)

    return(csv_df)