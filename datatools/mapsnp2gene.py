import pandas as pd
import pickle
import numpy as np


# MAPSNP2GENE creates snp to gene matrix and saves it to a file.
#
# INPUTS:
#   snpAnnotation: path to Plink snp annotation file in .bim format.
#   geneAnnotation: path to gene annotation file.
#   mappingDistance: snp to gene mapping distance
#   option: saving mode for snp-gene map
#   outfile: file name for saving the results.
#
# OUTPUTS:
#   <outfile>.pkl -This pickle file contains a snp to gene mapping in the DataFrame format.


def mapsnp2gene(snpAnnotation, geneAnnotation, mappingDistance, option, outfile):

    # Creating SNP dataframe from snp annotation file.
    snph = ['chromo', 'snprs', 'tmp1', 'snploc', 'tmp2', 'tmp3']
    sdf = pd.read_csv(snpAnnotation, sep=r"\s+", names=snph, engine='python')
    sdf['chromo'] = pd.to_numeric(sdf['chromo'])

    # Creating gene dataframe from gene annotation file.
    geneh = ['chromo', 'geneloc1', 'geneloc2', 'genes']
    gdf = pd.read_csv(geneAnnotation, sep=r"\s+", names=geneh, engine='python')
    gdf = gdf[gdf.chromo.apply(lambda x: x.isnumeric())]
    gdf['chromo'] = pd.to_numeric(gdf['chromo'])
    gdf.sort_values(by='chromo', inplace=True)

    # Expanding gene window by subtracting and adding from start and end loci.
    gdf['geneloc1'] = gdf['geneloc1'] - mappingDistance
    gdf['geneloc2'] = gdf['geneloc2'] + mappingDistance

    # Doing an outer join to get all genes and snp listed by chromosome.
    cdf = gdf.merge(sdf, how='outer', on='chromo')

    # Filtering out entries where snp isn't located between start and end loci.
    cdf = cdf[cdf.apply(
    lambda x: filter_val(x['geneloc1'], x['geneloc2'], x['snploc']), axis=1)]

    # Creating list of unique rsids from filtered results.
    snplist = cdf['snprs'].drop_duplicates()

    # Option chosen to save to snplist.
    if (option == 'snplist'):

        # Saving SNPlist to pickle file.
        final = open(outfile, 'wb')
        pickle.dump(snplist, final)
        final.close()

    # Option chosen to save to matrix.
    elif (option == 'matrix'):

        # Getting list of unique genes from filtered results.
        genelist = cdf['genes'].drop_duplicates()

        # Creating dataframe of appropriate size, and setting labels.
        sgm = pd.DataFrame(np.zeros((len(snplist), len(genelist))),
                                index=snplist, columns=genelist, dtype=int)


        # Setting snp-gene matrix values to true if snp is within gene window.
        for row in cdf.values:
            sgm[row[3]][row[4]] = 1

        # Saving snp-gene matrix to pickle file.
        final = open(outfile, 'wb')
        pickle.dump(sgm, final)
        final.close()

    else:
        # Output option not recognized.
        print("Return option error, valid options are 'snplist', or 'matrix'")


# Filter to keep snps with loci between gene's lower and upper range.
def filter_val(lower, upper, locus):
    return ((lower <= locus) and (locus <= upper))
