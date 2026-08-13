"""
imports the high level api functions so they are accessible
"""

from mutmodels.dirichlet_multinomial.api import (
    dm_onesample_SBS,
    dm_twosample_SBS,
)

__all__ = [
    "dm_twosample_SBS",
    "dm_onesample_SBS",
]