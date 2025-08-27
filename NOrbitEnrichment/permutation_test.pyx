import numpy as np
cimport numpy as np

def single_permutation(np.ndarray[np.int32_t, ndim=2] data, int numCellTypes, 
                       np.ndarray[np.int32_t, ndim=2] one_indices):
    
    cdef int n_rows = data.shape[0]
    cdef int n_cols = data.shape[1]
    
    cdef np.ndarray[np.int32_t, ndim=2] permuted_data = np.zeros_like(data)
    cdef np.ndarray[np.int32_t, ndim=2] zero_matrix = np.zeros((n_rows, n_cols - numCellTypes), dtype=np.int32)

    permuted_data[:, :numCellTypes] = data[:, :numCellTypes]
    permuted_data[:, numCellTypes:] = np.copy(zero_matrix)  # Reset to zero matrix for each permutation

    permuted_one_indices = np.random.permutation(one_indices[:, 0])
    
    permuted_data[permuted_one_indices, one_indices[:, 1] + numCellTypes] = 1
    
    return permuted_data

def enrichment_analysis(np.ndarray[np.int32_t, ndim=2] data, int n_permutations, int numCellTypes, int sample_size):
    cdef int n_rows = data.shape[0]
    cdef int n_cols = data.shape[1]
    cdef int perm, i
    
    # Dictionary to count unique rows and their occurrences
    cdef dict row_counts = {}
    cdef dict null_greater = {}
    cdef tuple row_as_tuple

    # Initialize null_greater with all unique vectors from the original data
    for i in range(n_rows):
        row_as_tuple = tuple(data[i])
        if row_as_tuple not in null_greater:
            null_greater[row_as_tuple] = n_permutations

    # Perform permutations
    for perm in range(n_permutations):
        if perm % 100 == 0:
            print(perm)
        # Subsample vectors
        subsample_indices = np.random.randint(0, n_rows, sample_size).astype(np.int32)
        subsample_data = data[subsample_indices]

        # Precompute indices for ones in the subsample
        one_indices = np.argwhere(subsample_data[:, numCellTypes:] == 1).astype(np.int32)

        # Count subsampled unique row occurrences
        row_counts.clear()
        for i in range(sample_size):
            row_as_tuple = tuple(subsample_data[i])
            if row_as_tuple in row_counts:
                row_counts[row_as_tuple] += 1
            else:
                row_counts[row_as_tuple] = 1
        
        permuted_counts = {}  # Clear permuted_counts for each permutation
        permuted_data = single_permutation(subsample_data, numCellTypes, one_indices)
        
        for i in range(sample_size):
            row_as_tuple = tuple(permuted_data[i])
            if row_as_tuple in permuted_counts:
                permuted_counts[row_as_tuple] += 1
            else:
                permuted_counts[row_as_tuple] = 1

        # Compare permuted_counts with row_counts
        for row in row_counts:
            if row in permuted_counts and row_counts[row] > permuted_counts[row]:
                null_greater[row] -= 1

    # Calculate p-values for all unique rows in null_greater
    cdef np.ndarray[np.float32_t, ndim=1] p_values = np.zeros(n_rows, dtype=np.float32)
    for i in range(n_rows):
        row_as_tuple = tuple(data[i])
        p_values[i] = null_greater[row_as_tuple] / n_permutations

    return p_values
