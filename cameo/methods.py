# Copyright 2013 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import logging
import random
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
from cameo.basics import production_envelope
from deap import creator, base, tools, algorithms
import inspyred


class StrainDesignMethod(object):

    def __init__(self, *args, **kwargs):
        super(StrainDesignMethod, self).__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        self.run(*args, **kwargs)

    def run(self, *args, **kwargs):
        raise NotImplementedError


class HeuristicStrainOptimization(StrainDesignMethod):

    """docstring for HeuristicStrainOptimization"""

    def __init__(self, objective=None, termination_condition=False, *args, **kwargs):
        super(HeuristicStrainOptimization, self).__init__(*args, **kwargs)

    def run(self):
        pass


class HeuristicKnockoutSearch(HeuristicStrainOptimization):
    pass


class OptGene(HeuristicStrainOptimization):

    """Implements OptGene.

    
    }
    """


class DifferentialFVA(StrainDesignMethod):

    """Differential flux variability analysis.__init__.py

    Compare flux ranges of a reference model to and a set of models that
    are parameterized to lie on a grid of evenly space points in the
    production envelope of a model.

    ...
    """

    def __init__(self, design_space_model=None, reference_model=None, target=None, variables=[], points=20):
        super(DifferentialFVA, self).__init__()
        self.design_space_model = design_space_model
        self.reference_model = reference_model
        self.target = target
        self.variables = variables
        self.points = points
        self.envelope = None
        self.grid = None
        self._init()

    def _init(self):
        # Mathematica code
        # maxDistanceBetweenLbAndUb=SortBy[Thread[envelope[[{2,3}]]],Abs[Subtract@@#]&][[-1]];
        # stepSize=(#[[2]]-#[[1]])/(OptionValue["Points"]-1)&[maxDistanceBetweenLbAndUb];
        # points2scan=Flatten[Table[NestWhileList[{envelope[[1,i]],#[[2]]+stepSize}&,{envelope[[1,i]],envelope[[2,i]]},#[[2]]<(envelope[[3,i]]-stepSize)&],{i,1,Length[envelope[[2]]]}],1];
        # points2scan=Union[points2scan,Thread[envelope[[{1,2}]]],Thread[envelope[[{1,3}]]]];

        logger.info('Running DifferentialFVA initialization ...')
        self.envelope = production_envelope(
            self.design_space_model, self.target, self.variables, points=self.points)
        zipped_envelope = zip(*self.envelope.itervalues())
        max_interval_in_envelope = sorted(zipped_envelope, key=lambda xy: abs(xy[0]-xy[1]), reverse=True)[0]
        step_size = (max_interval_in_envelope[1] - max_interval_in_envelope[0])/(self.points - 1)
        self.grid = list()
        for lb, ub in zipped_envelope:
            step = 0
            while step < ub:
                self.grid.append((lb, lb+step))
                step += step_size
        self.grid += zipped_envelope
        logger.info('...')

    def run(self):
        return self.envelope

if __name__ == '__main__':
    from optlang import Objective
    from cobra.io import read_sbml_model
    from cobra.flux_analysis import flux_variability_analysis
    from cameo.io import load_model

    model = load_model(
        '/Users/niko/Arbejder/Dev/cameo/tests/data/EcoliCore.xml')
    # print model.solver.variables
    # print model.optimize().f

    bpcy_objective = Objective('bpcy', model.solver.variables[
                               'Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2'] + model.solver.variables['EX_ac_LPAREN_e_RPAREN_'])
    print bpcy_objective


    # from cobra.oven.aliebrahim.designAnalysis import plot_production_envelope
    # import os
    # print os.getcwd()
    # from IPython import parallel
    # rc = parallel.Client()
    # direct_view = rc.direct_view()


    # model.optimize()
    # model.solution.status
    # max_growth = model.solution.f

    # reference_model = model.copy()
    # biomass_rxn = model.reactions.get_by_id('Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2')
    # reference_model.reactions.get_by_id('Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2').lower_bound = max_growth

    # diffFVA = DifferentialFVA(design_space_model=model,
    #     reference_model=reference_model,
    #     target='EX_succ_LPAREN_e_RPAREN_',
    #     variables=['Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2']
    #     )
    # diffFVA.run()



    # method = StrainDesignMethod(model=model)
    # method.run()
