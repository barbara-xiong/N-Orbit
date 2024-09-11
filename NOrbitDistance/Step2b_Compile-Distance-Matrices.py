import pandas as pd
import glob
import numpy as np

# Input paths and column name settings
intermediate_path = "/path/to/intermediates/SyntheticV1/"
inputs = pd.read_csv(intermediate_path+ "neighborhood_norbits.csv")
cells = pd.read_csv("/path/to/data/synthetic_mrf_neighborhoods_v1.csv")
im_label = "Image"
neighborhood_label = "Neighborhood"

# Compile neighborhood distance vectors into matrix
neighborhood_list = sorted(set(inputs["image_neighborhood"]))
dfs = []
for file in glob.glob(intermediate_path + "neighborhood_dists/*.csv"):
    dfs.append(pd.read_csv(file).drop(["Unnamed: 0"], axis=1))
df = pd.concat(dfs, axis = 1)
df = df.loc[:,~df.columns.duplicated()].copy()
df = df.reindex(sorted(df.columns), axis=1)
df = pd.DataFrame(df.values+np.transpose(df.values), index=list(df.columns),columns=list(df.columns))
df.to_csv(intermediate_path+"neighborhood_distance_matrix.csv")
