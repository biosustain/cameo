from random import Random
from random import sample, randint
import inspyred
from inspyred.ec import SA
from cameo import load_model

model = load_model("../tests/data/EcoliCore.xml")
def eval_individual(individual, args):

    rxns_to_ko = [model.reactions[index] for index in individual]
    original_bounds = [(reaction.lower_bound, reaction.upper_bound) for reaction in rxns_to_ko]
    for reaction in rxns_to_ko:
        reaction.lower_bound = 0
        reaction.upper_bound = 0
    solution = model.optimize()
    for rxn, bounds in zip(rxns_to_ko, original_bounds):
        rxn.lower_bound = bounds[0]
        rxn.upper_bound = bounds[1]
    if solution.status == 'optimal':
        try:
            fitness = (solution.f * model.solution.x_dict['EX_ac_LPAREN_e_RPAREN_'])/abs(model.solution.x_dict['EX_glc_LPAREN_e_RPAREN_'])
        except ZeroDivisionError:
            fitness = 0
    else:
        fitness = 0

    #print "B: %s\tP: %s\tU: %s\tF: %s" % (solution.f, model.solution.x_dict['EX_ac_LPAREN_e_RPAREN_'], model.solution.x_dict['EX_glc_LPAREN_e_RPAREN_'], fitness)

    return fitness

def eval_pop(candidates, args):
    return [eval_individual(i, args) for i in candidates]

def custom_generator(random, args):
    max_size = args.get('max_size', 9)
    individual = sample(xrange(len(model.reactions)), random.randint(1, max_size))
    return individual

def my_observer(population, num_generations, num_evaluations, args):
    best = max(population)
        # print('{0:6} -- {1} : {2}'.format(num_generations,
        #                               best.fitness,
        #                               str(best.candidate)))
    print('%s -- %s -- %s' % (num_generations, best.fitness, best))

@inspyred.ec.variators.mutator
def custom_mutation(random, individual, args):
    new_individual = list()
    for index in individual:
        if random.random() < args.get('mutation_rate', .1):
            new_individual.append(randint(0, len(model.reactions)-1))
        else:
            new_individual.append(index)

    new_individual.append(35)
    return list(set(new_individual))

r = Random()
ga = SA(r)
ga.observer = my_observer
ga.variator = [
    inspyred.ec.variators.n_point_crossover,
    custom_mutation
]
ga.terminator = inspyred.ec.terminators.evaluation_termination

final_pop = ga.evolve(generator=custom_generator,
                      evaluator=eval_pop,
                      pop_size=100,
                      maximize=True,
                      max_evaluations=30000,
                      max_generations=300,
                      mutation_rate=0.4)

print final_pop
print ga.termination_cause
print ga.num_evaluations
