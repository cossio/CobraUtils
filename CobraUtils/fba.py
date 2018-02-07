from cobra.flux_analysis import flux_variability_analysis

from .utils import remove_null_reactions


def set_fva_bounds(cobra_model):
    """
    Sets reaction bounds using FVA.
    """

    fva = flux_variability_analysis(cobra_model, fraction_of_optimum=0)
    
    for rxn in cobra_model.reactions:
        assert rxn.lower_bound <= rxn.upper_bound
        rxn.lower_bound = fva['minimum'][rxn.id]
        rxn.upper_bound = max(fva['maximum'][rxn.id], fva['minimum'][rxn.id])


def reduce_model_fva(cobra_model, tol=0):
    """
    Sets FVA bounds and removes null reactions
    """
    set_fva_bounds(cobra_model)
    remove_null_reactions(cobra_model, tol)


def optimize_expression(cobra_model, expr, sense='max'):
    "Optimize an expression"
    assert sense == 'max' or sense == 'min'

    old_objective = cobra_model.objective
    old_direction = cobra_model.objective_direction

    cobra_model.objective = expr
    cobra_model.objective_direction = sense
    f = cobra_model.optimize().f

    cobra_model.objective = old_objective
    cobra_model.objective_direction = old_direction

    return f


def minmax_expression(cobra_model, expr):
    "Min/max of an expression"
    fmax = optimize_expression(cobra_model, expr, sense='max')
    fmin = optimize_expression(cobra_model, expr, sense='min')
    return fmin, fmax
