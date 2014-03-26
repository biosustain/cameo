import time
from copy import copy
from cameo.io import load_model
from cameo.basics import flux_variability_analysis, production_envelope
from cameo.methods import DifferentialFVA
from cameo.config import ViewFacade
from multiprocessing import Pool

model = load_model('../tests/data/iJO1366.pickle')
reference_model = copy(model)
reference_model.reactions.get_by_id('Ec_biomass_iJO1366_core_53p95M').lower_bound = 0.9823718127270133

# t1 = time.time()
# fva_solution = flux_variability_analysis(model)
# t2 = time.time()
# print "Execution time: %s" % (t2-t1)

# t1 = time.time()
# fva_solution = flux_variability_analysis(model, view=ViewFacade())
# t2 = time.time()
# print "Execution time: %s" % (t2-t1)

# t1 = time.time()
# fva_solution = flux_variability_analysis(model, )
# t2 = time.time()
# print "Execution time: %s" % (t2-t1)



dFVA = DifferentialFVA(design_space_model=model,
        reference_model=reference_model,
        target='EX_succ_e',
        variables=['Ec_biomass_iJO1366_core_53p95M']
        )

dFVA.run()
    