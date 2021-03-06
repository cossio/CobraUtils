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


def remove_null_reactions(cobra_model, tol=0):
    """
    Remove reactions that have lb == ub == 0, or
    |ub - lb| <= tol.
    """
    assert tol >= 0

    rxns = [rxn for rxn in cobra_model.reactions
                if abs(rxn.lower_bound - rxn.upper_bound) <= tol]
    
    cobra_model.remove_reactions(rxns, remove_orphans=True)


def metabolite_producers(metabolite):
    "List of reactions that can produce metabolite"
    return [rxn for rxn in metabolite.reactions if 
            metabolite in rxn.products and rxn.upper_bound > 0 or
            metabolite in rxn.reactants and rxn.lower_bound < 0]

def metabolite_consumers(metabolite):
    "List of reactions that can consume metabolite"
    return [rxn for rxn in metabolite.reactions if 
            metabolite in rxn.reactants and rxn.upper_bound > 0 or
            metabolite in rxn.products and rxn.lower_bound < 0]

def active_metabolite_producers(metabolite, fluxes, tol=1e-10):
    "List of reactions that produce metabolite"
    assert 0 <= tol < float('inf')
    return [rxn for rxn in metabolite.reactions if
            metabolite in rxn.products  and fluxes[rxn.id] >  tol or
            metabolite in rxn.reactants and fluxes[rxn.id] < -tol]

def active_metabolite_consumers(metabolite, fluxes, tol=1e-10):
    "List of reactions that consume metabolite"
    assert 0 <= tol < float('inf')
    return [rxn for rxn in metabolite.reactions if
            metabolite in rxn.products  and fluxes[rxn.id] < -tol or
            metabolite in rxn.reactants and fluxes[rxn.id] >  tol]

def metabolite_production_expression(metabolite):
    """
    Expression giving overall production of a metabolite. 
    Equals the turnover if this is an internal metabolite.
    """
    expr = 0
    for rxn in metabolite.reactions:
        if metabolite in rxn.products:
            expr += rxn.forward_variable * rxn.get_coefficient(metabolite.id)
        else:
            expr -= rxn.reverse_variable * rxn.get_coefficient(metabolite.id)
    return expr

def metabolite_consumption_expression(metabolite):
    """
    Expression giving overall consumption of a metabolite. 
    Equals the turnover if this is an internal metabolite.
    """
    expr = 0
    for rxn in metabolite.reactions:
        if metabolite in rxn.products:
            expr += rxn.reverse_variable * rxn.get_coefficient(metabolite.id)
        else:
            expr -= rxn.forward_variable * rxn.get_coefficient(metabolite.id)
    return expr


def exchanges_secretion(cobra_model):
    "Subset of exchanges that can secrete."
    return [rxn for rxn in cobra_model.exchanges if 
            not rxn.products  and rxn.upper_bound > 0 or
            not rxn.reactants and rxn.lower_bound < 0]

def exchanges_consumption(cobra_model):
    "Subset of exchanges that can consume."
    return [rxn for rxn in cobra_model.exchanges if 
            not rxn.products  and rxn.lower_bound < 0 or
            not rxn.reactants and rxn.upper_bound > 0]


def element_exchange_expression(cobra_model, element):
    expr = 0
    for rxn in cobra_model.exchanges:
        for met, s in rxn.metabolites.items():
            if element in met.elements:
                expr += rxn.flux_expression * s * met.elements[element]
    return expr
