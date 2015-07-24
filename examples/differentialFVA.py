import copy
from cameo import load_model
from cameo.methods import DifferentialFVA
from cameo.basics import production_envelope

model = load_model('../tests/data/EcoliCore.xml')

solution = model.optimize()
surrogate_experimental_fluxes = dict([(key, (val * .8, val * 1.2)) for key, val in solution.x_dict.iteritems()])

constraint_model = copy.copy(model)
for reaction_id, (lb, ub) in surrogate_experimental_fluxes.iteritems():
    reaction = constraint_model.reactions.get_by_id(reaction_id)
    reaction.lower_bound = lb
    reaction.upper_bound = ub

diffFVA = DifferentialFVA(design_space_model=model,
                          reference_model=constraint_model,
                          target='EX_succ_LPAREN_e_RPAREN_',
                          variables=[model.reactions.Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2])

print diffFVA.run()
