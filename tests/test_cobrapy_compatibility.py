import types
import nose
from cobra.test import create_test_model
from cobra.test.unit_tests import CobraTestCase, TestReactions
from cobra.test.flux_analysis import TestCobraFluxAnalysis
from cameo.solver_based_model import to_solver_based_model, SolverBasedModel


def setUp(self):
    # Make Model pickable and then load a solver based version of test_pickle
    self.model = to_solver_based_model(create_test_model())
    self.model_class = SolverBasedModel

for cls in (CobraTestCase, TestReactions, TestCobraFluxAnalysis):
    cls.setUp = types.MethodType(setUp, cls)


if __name__ == '__main__':
    nose.runmodule()