import pandas as pd
import glob
import numpy as np

path = "/path/to/intermediates/SyntheticV1/"

dfs = []
for f in glob.glob(path+"norbits/neighborhood_norbits_*.csv"):
    df = pd.read_csv(f)
    dfs.append(df)
df = pd.concat(dfs)
df

df.to_csv(path+"neighborhood_norbits.csv")
