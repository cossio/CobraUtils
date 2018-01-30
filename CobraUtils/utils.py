"""
General utility functions
"""


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


