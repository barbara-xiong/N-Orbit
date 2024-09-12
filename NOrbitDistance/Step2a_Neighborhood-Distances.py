import pandas as pd
import numpy as np
from iteround import saferound
from scipy.spatial import distance_matrix, minkowski_distance
from munkres import Munkres
import plotly.express as px
from matplotlib import pyplot as plt
import sys
import os
from scipy.optimize import linear_sum_assignment

# If UNIT is image, use UNIT_MODE = "Image". If UNIT is neighborhood, use UNIT_MODE = "Neighborhood"
# To run all images in a single job, use UNIT_MODE = "All" (slow)
UNIT = sys.argv[1]
UNIT_MODE = "Image"

# Path to intermediates
intermediate_path = "/path/to/intermediates/SyntheticV1/"

df = pd.read_csv(intermediate_path + "neighborhood_norbits.csv")
neighborhood_list = sorted(set(df["image_neighborhood"]))

if UNIT_MODE == "Image":
	unit_neighborhoods = [item for item in neighborhood_list if item.startswith(UNIT)]
elif UNIT_MODE == "Neighborhood":
	unit_neighborhoods = [item for item in neighborhood_list if item == UNIT]
elif UNIT_MODE == "All":
	unit_neighborhoods = list(neighborhood_list)

# Create directory if not already exists
if not os.path.exists(intermediate_path+"neighborhood_dists/"):
    os.makedirs(intermediate_path+"neighborhood_dists/")


print(unit_neighborhoods)

def neighborhood_distance(df1_vectors, df2_vectors):
    distance_mat = distance_matrix(df1_vectors, df2_vectors, p = 1)
    distance_mat_orig = np.copy(distance_mat)
    row_ind, column_ind = linear_sum_assignment(distance_mat)
    return np.sum(distance_mat[row_ind,column_ind]), row_ind, column_ind

import warnings
warnings.filterwarnings("ignore")

def main():
    unit_dists = {}
    
    n_vectors = {}
    print("Starting vector writing.")
    for i in neighborhood_list:
            df1_vectors = pd.read_csv(intermediate_path + "neighborhood_vectors/"+i+".csv").set_index("Unnamed: 0")
            n_vectors[i] = df1_vectors.values
    print("Finished vector writing.")
    
    for i in unit_neighborhoods:
        unit_dists[i] = np.zeros((len(neighborhood_list)))
        print(i)
        for j in range(len(neighborhood_list)):
            if i > neighborhood_list[j]:
                dist, row_ind, column_ind = neighborhood_distance(n_vectors[i],n_vectors[neighborhood_list[j]]) 
                unit_dists[i][j] = dist
            else:
                break
        
    unit_dists = pd.DataFrame(unit_dists)
    unit_dists.to_csv(intermediate_path+"neighborhood_dists/" + UNIT+"_ndists.csv")

if __name__ == "__main__":
    main()
