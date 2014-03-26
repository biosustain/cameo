import pickle
import inspyred
from cobra.oven.phantomas1234.solver_based_model import to_solver_based_model
from optlang.inspyred_interface import Model, Variable, Objective


with open('../tests/data/iJO1366.pickle') as fhandle:
    model = to_solver_based_model(pickle.load(fhandle))

def eval_individual(individual):
    # print individual.values()
    rxns_to_ko = [model.reactions.get_by_id(key.name) for key, val in individual.iteritems() if val == 1]
    # print rxns_to_ko
    original_bounds = [(reaction.lower_bound, reaction.upper_bound) for reaction in rxns_to_ko]
    for reaction in rxns_to_ko:
        reaction.lower_bound = 0
        reaction.upper_bound = 0
    solution = model.optimize()
    for rxn, bounds in zip(rxns_to_ko, original_bounds):
        rxn.lower_bound = bounds[0]
        rxn.upper_bound = bounds[1]
    # print solution.f
    fitness = solution.f * model.solution.x_dict['EX_ac_e']
#     fitness = solution.f
    return fitness


objective = Objective(eval_individual, direction='max')
# heuristic_modelassertIn(assertIn(first, second, msg=message), second, msg=message) = Model(objective=objective, algorithm='PSO')
heuristic_model = Model(objective=objective)

@inspyred.ec.generators.diversify
def custom_generator(random, args):
    individual = [0 for i in xrange(len(model.reactions))]
    choice = random.choice(xrange(len(model.reactions)))
    individual[choice] = 1
    return individual
heuristic_model._generator = custom_generator

def my_observer(population, num_generations, num_evaluations, args):
    best = max(population)
        # print('{0:6} -- {1} : {2}'.format(num_generations, 
        #                               best.fitness, 
        #                               str(best.candidate)))
    print('{0:6} -- {1}'.format(num_generations, best.fitness))

heuristic_model.observer = my_observer
heuristic_model.add([Variable(reaction.id, lb=0, ub=1, type='binary') for reaction in model.reactions])
final_pop = heuristic_model.optimize()
print best(final_pop)
