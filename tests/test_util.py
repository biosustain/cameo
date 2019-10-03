# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import absolute_import, print_function

from functools import partial
from itertools import chain

import pytest
from cobra import Metabolite

from cameo.network_analysis.util import distance_based_on_molecular_formula
from cameo.util import (ProblemCache, RandomGenerator, Singleton, TimeMachine,
                        float_ceil, float_floor, frozendict, generate_colors,
                        partition)

SEED = 1234


@pytest.fixture(scope="function")
def problem_cache_trial(core_model):
    reference = core_model.optimize().fluxes
    n_constraints = len(core_model.solver.constraints)
    n_variables = len(core_model.solver.variables)
    return core_model, reference, n_constraints, n_variables


class TestTimeMachine:
    def test_one_change_list(self):
        tm = TimeMachine()
        l = [1, 2, 3, 4]
        tm(do=partial(l.append, 5), undo=l.pop)
        assert l == [1, 2, 3, 4, 5]
        tm.reset()
        assert l == [1, 2, 3, 4]

    def test_str_handles_different_types_of_stored_operations(self):
        tm = TimeMachine()

        def normal_function():
            pass

        partial_function = partial(str, 1)
        tm(do=normal_function, undo=partial_function)
        assert tm.__str__().split('\n')[2:-1] == ["undo: " + str(str) + " (1,) {}", 'redo: normal_function']

    def test_with_statement(self):
        l = [1, 2, 3, 4]
        with TimeMachine() as tm:
            tm(do=partial(l.append, 33), undo=partial(l.pop))
            tm(do=partial(l.append, 66), undo=partial(l.pop))
            tm(do=partial(l.append, 99), undo=partial(l.pop))
        assert l == [1, 2, 3, 4]


def some_method_that_adds_stuff(model, cache):
    assert isinstance(cache, ProblemCache)

    def create_variable(model, var_id, lb, ub):
        return model.solver.interface.Variable(var_id, ub=ub, lb=lb)

    def update_variable(model, var, lb, ub):
        var.lb = lb
        var.ub = ub

    def create_constraint(model, cid, vars, lb, ub):
        return model.solver.interface.Constraint(vars[0] + vars[1], ub=ub, lb=lb)

    def update_constraint(model, constraint, vars, lb, ub):
        constraint.lb = lb
        constraint.ub = ub

    for i in range(10):
        cache.add_variable("var_%i" % (i + 1), create_variable, update_variable, 10, 15)

    for i in range(9):
        v1 = cache.variables["var_%i" % (i + 1)]
        v2 = cache.variables["var_%i" % (i + 2)]
        cache.add_constraint("c_%i" % (i + 1), create_constraint, update_constraint, [v1, v2], -20, 100)


class TestProblemCache:
    def test_add_variable(self, core_model):
        cache = ProblemCache(core_model)

        def add_var(model, var_id):
            return model.solver.interface.Variable(var_id, ub=0)

        def update_var(model, var):
            return setattr(var, "ub", 1000)
        for i in range(10):
            cache.add_variable("%i" % i, add_var, update_var)

        for i in range(10):
            assert cache.variables["%i" % i] in core_model.solver.variables
            assert cache.variables["%i" % i].ub == 0
            assert core_model.solver.variables["%i" % i].ub == 0

        for i in range(10):
            cache.add_variable("%i" % i, add_var, update_var)
            assert cache.variables["%i" % i].ub == 1000
            assert core_model.solver.variables["%i" % i].ub == 1000

        cache.reset()

        for i in range(10):
            with pytest.raises(KeyError):
                core_model.solver.variables.__getitem__("%i" % i)

    def test_add_constraint(self, core_model):
        cache = ProblemCache(core_model)

        def add_var(model, var_id):
            return model.solver.interface.Variable(var_id, ub=0)

        def add_constraint(m, const_id, var):
            return m.solver.interface.Constraint(var, lb=-10, ub=10, name=const_id)

        def update_constraint(model, const, var):
            return setattr(const, "ub", 1000)

        for i in range(10):
            cache.add_variable("%i" % i, add_var, None)
            cache.add_constraint("c%i" % i, add_constraint, update_constraint, cache.variables["%i" % i])

        for i in range(10):
            assert cache.constraints["c%i" % i] in core_model.solver.constraints
            assert cache.constraints["c%i" % i].ub == 10
            assert cache.constraints["c%i" % i].lb == -10
            assert core_model.solver.constraints["c%i" % i].ub == 10
            assert core_model.solver.constraints["c%i" % i].lb == -10

        for i in range(10):
            cache.add_constraint("c%i" % i, add_constraint, update_constraint, cache.variables["%i" % i])
            assert core_model.solver.constraints["c%i" % i].ub == 1000

        cache.reset()

        for i in range(10):
            with pytest.raises(KeyError):
                core_model.solver.variables.__getitem__("%i" % i)
            with pytest.raises(KeyError):
                core_model.solver.constraints.__getitem__("c%i" % i)

    def test_cache_problem(self, problem_cache_trial):
        core_model, reference, n_constraints, n_variables = problem_cache_trial
        # After the number of variables and constraints remains the same if nothing happens
        assert n_constraints == len(core_model.solver.constraints)
        assert n_variables == len(core_model.solver.variables)

        cache = ProblemCache(core_model)
        some_method_that_adds_stuff(core_model, cache)
        # After running some_method_that_adds_stuff with cache, problem has 10 more variables
        assert n_variables + 10 == len(core_model.solver.variables)
        # And has 9 more more constraints
        assert n_constraints + 9 == len(core_model.solver.constraints)

        cache.reset()
        # After reset cache, the problem should return to its original size
        assert n_constraints == len(core_model.solver.constraints)
        assert n_variables == len(core_model.solver.variables)

    def test_with(self, problem_cache_trial):
        core_model, reference, n_constraints, n_variables = problem_cache_trial
        with ProblemCache(core_model) as cache:
            some_method_that_adds_stuff(core_model, cache)
            # After running some_method_that_adds_stuff with cache, problem has 10 more variables
            assert n_variables + 10 == len(core_model.solver.variables)
            # And has 9 more more constraints
            assert n_constraints + 9 == len(core_model.solver.constraints)

            # If the method runs again, it does not add repeated variables
            some_method_that_adds_stuff(core_model, cache)
            # After running some_method_that_adds_stuff with cache, problem has 10 more variables
            assert n_variables + 10 == len(core_model.solver.variables)
            # And has 9 more more constraints
            assert n_constraints + 9 == len(core_model.solver.constraints)

        # After reset cache, the problem should return to its original size
        assert n_constraints == len(core_model.solver.constraints)
        assert n_variables == len(core_model.solver.variables)


class TestRandomGenerator:
    def test_random(self):
        random = RandomGenerator()
        for _ in range(1000):
            assert random.random() >= 0
            assert random.random() <= 1

    def test_randint(self):
        random = RandomGenerator()
        lower = 0
        upper = 10
        for _ in range(10000):
            assert random.randint(lower, upper) >= lower
            assert random.randint(lower, upper) <= upper

        lower = -10
        upper = 100
        for _ in range(10000):
            assert random.randint(lower, upper) >= lower
            assert random.randint(lower, upper) <= upper

        lower = 5
        upper = 21
        for _ in range(10000):
            assert random.randint(lower, upper) >= lower
            assert random.randint(lower, upper) <= upper

        lower = -5
        upper = 5
        for _ in range(10000):
            assert random.randint(lower, upper) >= lower
            assert random.randint(lower, upper) <= upper

    def test_seeded_methods(self):
        random = RandomGenerator()

        random.seed(SEED)
        value = random.random()
        random.seed(SEED)
        assert value == random.random()

        random.seed(SEED)
        value = random.randint(1, 10)
        random.seed(SEED)
        assert value == random.randint(1, 10)

        random.seed(SEED)
        population = [1, 2, 3, 4, 5]
        value = random.sample(population, 2)
        random.seed(SEED)
        assert value == random.sample(population, 2)

        random.seed(SEED)
        value = random.uniform()
        random.seed(SEED)
        assert value == random.uniform()


class TestUtils:
    def test_color_generation(self):
        for i in range(1, 100):
            color_map = generate_colors(i)
            assert len(color_map) == i
            assert len(color_map) == len(set(color_map.values()))

    def test_partition(self):
        chunks = 3
        iterables = [
            [1, 2, 3, 4, 5, 6, 7, 8, 9],
            {5, 3, 8, 3, 8, 5, 8, 0, 10, 11, 15},
            range(29)
        ]
        for fixture in iterables:
            test_output = partition(fixture, chunks)
            assert len(fixture) == sum(map(len, test_output))
            assert len(test_output) == chunks
            assert list(fixture) == list(chain(*test_output))
            for out_chunk in test_output:
                assert set(out_chunk).issubset(set(fixture))

        bad_input = 5
        with pytest.raises(TypeError):
            partition(bad_input, chunks)

    def test_distance_based_on_molecular_formula(self):  # from network_analysis.util
        met1 = Metabolite("H2O", formula="H2O")
        met2 = Metabolite("H2O2", formula="H2O2")
        met3 = Metabolite("C6H12O6", formula="C6H12O6")

        assert distance_based_on_molecular_formula(met1, met2, normalize=False) == 1
        assert distance_based_on_molecular_formula(met1, met2, normalize=True) == 1. / 7

        assert distance_based_on_molecular_formula(met2, met3, normalize=False) == 20
        assert distance_based_on_molecular_formula(met2, met3, normalize=True) == 20. / 28

        assert distance_based_on_molecular_formula(met1, met3, normalize=False) == 21
        assert distance_based_on_molecular_formula(met1, met3, normalize=True) == 21. / 27

    def test_float_conversions(self):
        val = 1.3456
        new_value = float_floor(val, 2)
        assert new_value == 1.34
        new_value = float_ceil(val, 2)
        assert new_value == 1.35
        new_value = float_floor(val, 1)
        assert new_value == 1.3
        new_value = float_ceil(val, 1)
        assert new_value == 1.4
        new_value = float_floor(val)
        assert new_value == 1
        new_value = float_ceil(val)
        assert new_value == 2
        val = 0.00000
        for i in range(1, 10):
            new_value = float_floor(val, i)
            assert new_value == 0
            new_value = float_ceil(val, i)
            assert new_value == 0


class TestFrozendict:
    def test_frozen_attributes(self):
        frozen_dict = frozendict({"A": 1, "B": 2, "C": 3, "D": 4, "E": [2, 3, 4, 5]})
        with pytest.raises(AttributeError):
            frozen_dict.popitem()
        with pytest.raises(AttributeError):
            frozen_dict.pop("A")
        with pytest.raises(AttributeError):
            frozen_dict.__setitem__("C", 1)
        with pytest.raises(AttributeError):
            frozen_dict.setdefault("K")
        with pytest.raises(AttributeError):
            frozen_dict.__delitem__("A")
        with pytest.raises(AttributeError):
            frozen_dict.update()

        assert hasattr(frozen_dict, "__hash__")


class TestSingleton:
    def test_singleton(self):
        s1 = Singleton()
        s2 = Singleton()
        assert s1 is s2
