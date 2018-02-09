"""
API to interact with the iCHOv1_K1 model.
"""

from CobraUtils.utils import (
    exchange_reaction, remove_null_reactions,
    metabolite_consumers, metabolite_producers, element_exchange_expression,
    metabolite_consumption_expression, metabolite_production_expression,
    exchanges_secretion, exchanges_consumption,
    active_metabolite_producers, active_metabolite_consumers
)

from CobraUtils.fba import (
    set_fva_bounds, reduce_model_fva, optimize_expression, minmax_expression, flux_constraint
)
