"""
General utility functions
"""

import cobra


def exchange_reaction(cobra_model, metabolite):
    """
    Returns the exchange reaction for a metabolite
    in a COBRA model, or None.
    """

    for rxn in metabolite.reactions:
        if not rxn.reactants or not rxn.products:
            return rxn


def flux_constraint(cobra_model, reactions, coefficients_forward, coefficients_backward):
    """
    Adds a linear constrain of the form
        0 <= cf[1] vf[1] + cb[1] vb[1] + ... <= 1
    where vf[i], vb[i] are the forward and backward 
    fluxes of reaction i, and 
    cf = coefficients_forward, cb = coefficients_backward.
    """

    flux_expression = 0
    for rxn in reactions:
        assert rxn in cobra_model.reactions
        assert 0 <= coefficients_backward[rxn] < float('inf')
        assert 0 <= coefficients_forward[rxn]  < float('inf')
        flux_expression += coefficients_backward[rxn] * rxn.backward_variable
        flux_expression += coefficients_forward[rxn]  * rxn.forward_variable
    
    cobra_model.problem.Constraint(flux_expression, lb=0, ub=1)


def remove_null_reactions(cobra_model):
    """
    Remove reactions that have lb == ub == 0.
    """

    rxns = [rxn for rxn in cobra_model.reactions
                if rxn.lower_bound == rxn.upper_bound == 0]
    
    cobra_model.remove_reactions(rxns, remove_orphans=True)


def set_fva_bounds(cobra_model):
    """
    Sets reaction bounds using FVA
    """

    fva = cobra.flux_analysis.flux_variability_analysis(cobra_model)
    
    for rxn in cobra_model.reactions:
        rxn.lower_bound = fva['minimum'][rxn.id]
        rxn.upper_bound = fva['maximum'][rxn.id]
        