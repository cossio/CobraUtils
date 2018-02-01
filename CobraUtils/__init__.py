"""
API to interact with the iCHOv1_K1 model.
"""

from CobraUtils.utils import (
    exchange_reaction, flux_constraint, 
    remove_null_reactions, set_fva_bounds, reduce_model_fva)
