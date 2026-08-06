"""
imports the high level api functions so they are accessible
"""

# example:
# from mutmodels.dirichlet_multinomial.api import fit_from_file

import pandas as pd

import mutmodels.dirichlet_multinomial.core as dm
from mutmodels.common.parsers import read_matrix,align
