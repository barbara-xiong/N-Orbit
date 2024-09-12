import pandas as pd
import glob
import numpy as np

intermediate_path = "/path/to/intermediates/SyntheticV1/"

dfs = []
for f in glob.glob(intermediate_path+"norbits/neighborhood_norbits_*.csv"):
    df = pd.read_csv(f)
    dfs.append(df)
df = pd.concat(dfs)
df

df.to_csv(intermediate_path+"neighborhood_norbits.csv")
