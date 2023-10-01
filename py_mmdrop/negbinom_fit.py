from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np

def negbinom_fit(counts):

    def molecule_info(counts):
        # Print another message to indicate the progress
        print("Calculating variances ...")

        counts_matrix = counts.X
        if (counts_matrix < 0).sum().sum() > 0:
            raise ValueError("Expression matrix contains negative values! Please provide raw UMI counts!")

        if (counts_matrix >= 1).sum().sum() != (counts_matrix > 0).sum().sum():
            raise ValueError("Error: Expression matrix is not integers! Please provide raw UMI counts.")

        counts_matrix_transposed = counts_matrix.transpose().tocsc()
        counts = pd.DataFrame.sparse.from_spmatrix(counts_matrix_transposed, columns=counts.obs_names, index=counts.var_names)

        # Calculate total molecules per gene (tjs)
        tjs = counts.sum(axis=1, skipna=True)

        # Check for undetected genes
        no_detect = (tjs <= 0).sum()
        if no_detect > 0:
            raise ValueError(f"Error: Contains {no_detect} undetected genes.")

        # Calculate total molecules per cell (tis)
        tis = counts.sum(axis=0, skipna=True)

        # Check that all cells have at least one detected molecule
        if (tis <= 0).sum() > 0:
            raise ValueError("Error: All cells must have at least one detected molecule.")

        # Calculate observed dropouts per gene (djs)
        djs = len(counts.columns) - (counts > 0).sum(axis=1)

        # Calculate observed dropouts per cell (dis)
        dis = len(counts.index) - (counts > 0).sum(axis=0)

        # Get the number of cells (nc) and genes (ng)
        nc = len(counts.columns)
        ng = len(counts.index)

        # Calculate the total molecules sampled (total)
        total = tis.sum(skipna=True)

        return {
            'tis': tis.values.tolist(),
            'tjs': tjs.values.tolist(),
            'dis': dis.values.tolist(),
            'djs': djs.values.tolist(),
            'total': total,
            'nc': nc,
            'ng': ng
        }

    vals = molecule_info(counts)
    min_size = 1e-10

    counts_matrix = counts.X
    # Transpose the sparse counts_matrix
    counts_matrix_transposed = counts_matrix.transpose().tocsc()
    counts = pd.DataFrame.sparse.from_spmatrix(counts_matrix_transposed, columns=counts.obs_names, index=counts.var_names)
    
    for keys in vals.keys():
        vals[keys] = np.array(vals[keys], dtype=np.float32)

    batch_size = 1000  # Adjust as needed
    num_genes = len(counts)
    my_rowvars = np.empty(num_genes, dtype=np.float32)

    for i in range(0, num_genes, batch_size):
        batch_end = min(i + batch_size, num_genes)
    
        # Calculate mu_is for the current batch
        mu_is = vals['tjs'][i:batch_end][:, np.newaxis] * vals['tis'] / vals['total']
    
        # Compute my_rowvars for the current batch using vectorized operations
        my_rowvars[i:batch_end] = np.var(counts[i:batch_end] - mu_is, axis=1)
    
    size = (vals['tjs'] ** 2) * (np.sum(vals['tis'] ** 2) / vals['total'] ** 2) / ((vals['nc'] - 1) * my_rowvars - vals['tjs'])

    max_size = 10 * np.max(size)
    size[size < 0] = max_size
    size[size < min_size] = min_size

    return {
        'gene_names' : counts.index.tolist(),
        'var_obs': my_rowvars.tolist(),
        'sizes': size.tolist(),
        'vals': vals,
    }