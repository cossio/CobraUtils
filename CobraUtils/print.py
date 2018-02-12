import pandas, cobra


def show_limiting_exchanges(cobra_model, solution):
    "prints a table of exchange fluxes, highliting the ones that are limiting"

    def highlight_limiting(s):
        '''
        highlight limited fluxes
        '''
        is_lb_lim = s['flux'] == s['lb']
        is_ub_lim = s['flux'] == s['ub']
        return ['background-color: yellow' if is_lb_lim or is_ub_lim else '' for v in range(len(s))]


    df = pandas.DataFrame([[rxn.id, rxn.lower_bound, rxn.upper_bound, solution.fluxes[rxn.id], 
                           abs(1 - rxn.lower_bound / solution.fluxes[rxn.id]) <= 1e-2,
                           abs(1 - rxn.upper_bound / solution.fluxes[rxn.id]) <= 1e-2]
                           for rxn in cobra_model.exchanges if solution.fluxes[rxn.id] != 0],
                            columns=['id', 'lb', 'ub', 'flux', 'limit lb', 'limit ub']).style.apply(highlight_limiting, axis=1)

    return df
