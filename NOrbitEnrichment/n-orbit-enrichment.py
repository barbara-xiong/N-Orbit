import numpy as np
import pandas as pd
from permutation_test import enrichment_analysis
from scipy import stats
import sys

vectors_path = "/path/to/cluster/vectors/NR_vectors_predicted.csv"
trials_path = "/path/to/trials/"
run = sys.argv[1]

# Number of permutations
# If doing multiple runs, num_permutations needs to be the same for all runs
num_permutations = 25000

# Number of cell types (first numCellTypes columns are not permuted)
numCellTypes = 18

# Subsample size
sample_size = int(round(len(tcn_vectors_df)*0.2))

# Load data
tcn_vectors_df = pd.read_csv(vectors_path)
tcn_vectors_df = tcn_vectors_df.drop(columns="neighborhood")
tcn_vectors = tcn_vectors_df.values.astype(np.int32)

# Perform the enrichment analysis
p_values = enrichment_analysis(tcn_vectors, num_permutations, numCellTypes, sample_size)

# Add p-values to the original dataframe
tcn_vectors_df['pvalue'] = p_values
tcn_vectors_df = tcn_vectors_df.drop_duplicates()
tcn_vectors_df = tcn_vectors_df.sort_values("pvalue")
tcn_vectors_df['qvalue'] = stats.false_discovery_control(tcn_vectors_df['pvalue'], method='bh')

# Save the results to a CSV file
tcn_vectors_df.to_csv(results_path + "NOrbit_Enrichment_Run_"+str(run)+".csv")
