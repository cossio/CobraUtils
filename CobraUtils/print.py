from math import isclose
import pandas, cobra

from .fba import limiting_reactions


def show_limiting_reactions(cobra_model, solution, reactions=None):
    "prints a table of reaction fluxes, highliting the ones that are limiting"

    if reactions == None:
        reactions = cobra_model.exchanges
    
    lim_rxns = limiting_reactions(cobra_model, reactions)


    def highlight_limiting(s):
        '''
        highlight limited fluxes
        '''

        return ['background-color: yellow' if s['id'] in lim_rxns else '' 
                for v in range(len(s))]

    df = pandas.DataFrame([[rxn.id, rxn.lower_bound, rxn.upper_bound, solution.fluxes[rxn.id]]
                           for rxn in reactions if solution.fluxes[rxn.id] != 0],
                          columns=['id', 'lb', 'ub', 'flux']).style.apply(highlight_limiting, 
                                                                                                  axis=1)

    return df
