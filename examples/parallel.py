import pickle
import copy
import time
import cameo
import numpy as np
from cameo.util import flux_variability_analysis, flux_variability_analysis_parallel, production_envelope
from cobra.io import read_sbml_model
from cobra.oven.phantomas1234.solver_based_model import to_solver_based_model

from IPython import parallel
from IPython.parallel.util import interactive


@parallel.require(time)
def long_lasting_task(arg, arg2):
    time.sleep(.5)
    return (arg**2, arg/arg2)

# print direct_view.map(long_lasting_task, range(1,10), range(1,10))

@parallel.require(copy)
def test_remotely(model):
    original_objective = copy.copy(model.objective)
    model.objective = model.reactions.get_by_id('ENO')
    result = model.optimize().f
    model.objective = original_objective
    return result

# print rc[0].apply(test_remotely, model)


@parallel.require(flux_variability_analysis)
def test_parallel_fva(level, model=None):
    biomass_rxn = model.reactions.get_by_id('Ec_biomass_iJO1366_core_53p95M')
    biomass_rxn.lower_bound = level
    fva_result = flux_variability_analysis(model, reactions=('Ec_biomass_iJO1366_core_53p95M', 'ENO', 'PGK'))
    return [fva_result[key] for key in ('Ec_biomass_iJO1366_core_53p95M', 'ENO', 'PGK')]

# print rc[0].apply(test_parallel_fva, .3, model)


if __name__ == '__main__':


    from functools import partial
    # print pickle.loads(pickle.dumps(partial(test, 2)))(3)

    client = parallel.Client()
    client.block = True
    # rc[:].use_dill()
    direct_view = client.direct_view()

    @direct_view.remote(block=True)
    def getpid():
        import os
        return os.getpid()
    # import pdb; pdb.set_trace();
    print getpid()
   
    model = to_solver_based_model(pickle.load(open('../tests/data/iJO1366.pickle')))

    # direct_view['model'] = model
    # direct_view.scatter('reactions', model.reactions[0:100])
    # direct_view.execute('fva_sol = flux_variability_analysis(model, reactions)')
    # fva_sol = direct_view.gather('fva_sol')
    # print fva_sol

    # @direct_view.parallel(block=True)
    # @parallel.require(to_solver_based_model, pickle, flux_variability_analysis)
    # def test(reactions):
    #     from cameo.util import flux_variability_analysis
    #     model = to_solver_based_model(pickle.load(open('../tests/data/iJO1366.pickle')))
    #     return flux_variability_analysis(model, reactions=reactions)
    #     # return str([r.id for r in reactions])

    # seq = model.reactions[0:100]
    # print seq
    # test(seq)

    
    # print direct_view.map(test_parallel_fva, scan, [model]*len(scan))

    

    # the following does not work
    # print direct_view.map(partial(test_parallel_fva, model=model), list(np.linspace(0,.8, 10)))

    # direct_view.scatter('scan', scan)
    # direct_view['model'] = model
    # direct_view.execute("result = [test_parallel_fva(val, model) for val in scan]")
    # print direct_view.gather('result')

    # from cameo.util import flux_variability_analysis_parallel
    # flux_variability_analysis_parallel(model)

