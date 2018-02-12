from optlang.symbolics import Zero
import cobra


def flux_constraint(cobra_model, coefficients_forward, coefficients_reverse):
    """
    Adds a linear constrain of the form
        0 <= cf[1] vf[1] + cb[1] vb[1] + ... <= 1
    where vf[i], vb[i] are the forward and reverse fluxes of reaction i, and
    cf = coefficients_forward, cb = coefficients_backward are dictionaries
    """
    coefficients = dict()
    for (bigg_id, cf) in coefficients_forward.items():
        rxn = cobra_model.reactions.get_by_id(bigg_id)
        coefficients[rxn.forward_variable] = cf
    for (bigg_id, cr) in coefficients_reverse.items():
        rxn = cobra_model.reactions.get_by_id(bigg_id)
        coefficients[rxn.reverse_variable] = cr
        
    constraint = cobra_model.problem.Constraint(0, lb=0, ub=1)
    cobra_model.add_cons_vars(constraint)
    cobra_model.solver.update()
    constraint.set_linear_coefficients(coefficients=coefficients)


def set_enzymatic_objective(cobra_model, coefficients_forward, coefficients_reverse):
    """
    Sets objective of model to minimization of enzymatic mass.
    """
    coefficients = dict()
    for (bigg_id, cf) in coefficients_forward.items():
        rxn = cobra_model.reactions.get_by_id(bigg_id)
        coefficients[rxn.forward_variable] = cf
    for (bigg_id, cr) in coefficients_reverse.items():
        rxn = cobra_model.reactions.get_by_id(bigg_id)
        coefficients[rxn.reverse_variable] = cr
    
    cobra_model.objective = cobra_model.problem.Objective(Zero, 
                                                          direction='min', sloppy=True, 
                                                          name="min_enzymatic")

    cobra_model.objective.set_linear_coefficients(coefficients=coefficients)


def fba_and_min_enzyme(cobra_model, coefficients_forward, coefficients_reverse):
    """
    Performs FBA followed by minimization of enzyme content
    """

    with cobra_model as model:
        model.optimize()
        cobra.util.fix_objective_as_constraint(model)
        set_enzymatic_objective(model, coefficients_forward, coefficients_reverse)
        sol = cobra_model.optimize()
        return sol
