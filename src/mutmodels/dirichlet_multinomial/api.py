"""
high level api functions, made available through __init__.py
"""

import pandas as pd

from mutmodels.common.parsers import read_matrix,align
import mutmodels.dirichlet_multinomial.core as dm

# helper functions
def background_correct(data,background):
    """
    subtract background from counts, then remove nonnegative and round to counts
    """
    data = data.sub(background,axis=0)
    data[data<0] = 0
    data = data.round().astype(int)
    return data

def select_substitutions(data,substitutions):
    selection = []
    for sub in substitutions:
        cur = [i for i in data.index.values if sub in  i]

        selection+= (cur)

    # select the relevant substitutions without changing index order of selected 
    data = data[data.index.isin(selection)]
    return data

# main api functions
def dm_two_condition_SBS(matrix_fn,g1,g2,bg_fn=None,substitutions=None,matrix_type='SBS96',
                         n_bootstraps=1000,sig_level=0.05,dispersion_type='split',
                         studentize=True,stat_type='max',transform=None,rng=None):

    matrix = read_matrix(matrix_fn,matrix_type=matrix_type)
    
    data = matrix.loc[:,g1+g2]
    

    if bg_fn:
        bg = align(pd.read_csv(bg_fn,sep='\t',index_col=0).squeeze("columns"))
        data = background_correct(data,bg)

    if substitutions:
        data = select_substitutions(data,substitutions)

    # get the count array for the test. Test expects rows as replicates, columns as categories
    g1_counts = data[g1].values.T
    g2_counts = data[g2].values.T

    result = dm.dm_two_condition(g1_counts,g2_counts,n_bootstraps = n_bootstraps,
                              sig_level = sig_level,dispersion_type=dispersion_type,
                              studentize=studentize,stat_type=stat_type,
                              transform=transform,rng=rng)

    return result


