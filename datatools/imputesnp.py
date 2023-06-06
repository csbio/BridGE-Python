import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None

# INPUTESNP impute SNP data (0,1,2 format)
#
# This function imputes missing alleles in SNP dataset by
# replacing them with major-major.
#
# INPUTS:
#   data - SNP data matrix (0,1,2 format, -1 is missing value)
#   with rows are samples and columns are SNPs
#
# OUTPUTS:
#   data - new SNP data matrix

def imputesnp(data):
    # Iterating over SNPs.
    for column in data:
        # Count occurences of 0s, 1s, and 2s.
        sum0 = (data[column].values == 0).sum()
        sum1 = (data[column].values == 1).sum()
        sum2 = (data[column].values == 2).sum()

        # Replace nans (missing values) with most frequent value. (0, 1, or 2)
        gen = np.argmax([sum0, sum1, sum2])
        data[column] = data[column].replace(to_replace=(np.nan), value=gen)

    # Returns imputed data.
    return data
