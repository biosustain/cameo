from cobra.solvers import solver_dict
from cobra.test import create_test_model
from cobra.test.unit_tests import CobraTestCase, TestReactions
import types
import unittest
import nose
from cobra.test.flux_analysis import TestCobraFluxAnalysis
from optlang import Objective
from optlang.glpk_interface import Model as GlpkModel
from cameo import load_model
from cameo.solver_based_model import to_solver_based_model


ecoli_core_sbmodel = load_model('data/EcoliCore.xml')


class SolverBasedModelTestCase(unittest.TestCase):

    """Test features that are unique to SolverBasedModel,
    e.g.,objective specification and handling
    """

    def setUp(self):
        self.model = ecoli_core_sbmodel.copy()

    def test_objective(self):
        obj = self.model.objective
        self.assertEqual(
            obj.__str__(), 'Maximize\n1.0*Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2')

    def test_change_objective(self):
        self.model.objective = Objective(
            self.model.solver.variables['ENO'] + self.model.solver.variables['PFK'])
        self.assertEqual(self.model.objective.__str__(),
                         'Maximize\n1.0*ENO + 1.0*PFK')


# Test SolverBasedModel with cobrapy test suite

# solver_dict.pop('cplex')
# solver_dict.pop('gurobi')


def setUp(self):
    # Make Model pickable and then load a solver based version of test_pickle
    self.model = to_solver_based_model(create_test_model())

for cls in (CobraTestCase, TestReactions, TestCobraFluxAnalysis):
    cls.setUp = types.MethodType(setUp, cls)


if __name__ == '__main__':
    nose.runmodule()

    # Debug ...
    # from cobra.flux_analysis.variability import flux_variability_analysis
    # from cobra.flux_analysis.single_deletion import single_deletion, single_gene_deletion
    # from cobra.manipulation.modify import initialize_growth_medium
    # solver = GlpkModel()
    # cobra_model = to_solver_based_model(create_test_model(test_pickle), solver)
    # cobra_model.media_compositions['lb']
    # print cobra_model.optimize().f
    # print 1
    # initialize_growth_medium(cobra_model, 'LB')
    # print cobra_model.optimize().f
    # print 2
    # the_loci =  ['STM4081', 'STM0247', 'STM3867', 'STM2952']
    # the_genes = tpiA, metN, atpA, eno = map(cobra_model.genes.get_by_id, the_loci)
    # id_to_name = dict([(x.id, x.name) for x in the_genes])
    # growth_dict = {'fba':{tpiA.id:2.41, metN.id:2.44, atpA.id:1.87, eno.id:1.81},
    #                'moma':{ tpiA.id:1.62, metN.id:2.4, atpA.id:1.40, eno.id:0.33}}


    # element_list = the_loci
    # rates, statuses, problems = single_gene_deletion(cobra_model, element_list=element_list, method='fba')

    # cobra_model.reactions.get_by_id('biomass_iRR1083_metals').lower_bound = 2.4620629521939827
    # the_problem = 'return'
    # fva_out = flux_variability_analysis(cobra_model,
    #                                             the_problem=the_problem,
    # the_reactions=cobra_model.reactions[100:140])

    # objective_metabolite = Metabolite('objective')
    #
