import numpy as np
import pandas as pd
from permutation_test import enrichment_analysis
from scipy import stats
import sys

run = sys.argv[1]

# Load your data
tcn_vectors_df = pd.read_csv("/path/to/cluster/vectors/NR_vectors_predicted.csv")
tcn_vectors_df = tcn_vectors_df.drop(columns="neighborhood")

# Extract the numerical values and ensure the correct dtype
tcn_vectors = tcn_vectors_df.values.astype(np.int32)

# Number of permutations
# If doing multiple runs, num_permutations needs to be the same for all runs
num_permutations = 25000

# Number of cell types (first numCellTypes columns are not permuted)
numCellTypes = 18

# Subsample size
sample_size = int(round(len(tcn_vectors_df)*0.2))

# Perform the enrichment analysis
p_values = enrichment_analysis(tcn_vectors, num_permutations, numCellTypes, sample_size)

# Add p-values to the original dataframe
tcn_vectors_df['pvalue'] = p_values
tcn_vectors_df = tcn_vectors_df.drop_duplicates()
tcn_vectors_df = tcn_vectors_df.sort_values("pvalue")
tcn_vectors_df['qvalue'] = stats.false_discovery_control(tcn_vectors_df['pvalue'], method='bh')

# Save the results to a CSV file
tcn_vectors_df.to_csv("/path/to/trials/NOrbit_Enrichment_Run_"+str(run)+".csv")
