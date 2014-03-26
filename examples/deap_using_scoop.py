# not working, submitted issue at https://code.google.com/p/scoop/issues/detail?id=5&thanks=5&ts=1392730066
import array, random
from deap import creator, base, tools, algorithms
from scoop import futures

creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()
toolbox.register("attr_bool", random.randint, 0, 1)
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, n=10)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

# evalOneMax = lambda individual: [sum(individual)]
def evalOneMax(individual):
    return sum(individual),

toolbox.register("evaluate", evalOneMax)
toolbox.register("mate", tools.cxTwoPoints)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.10)
toolbox.register("select", tools.selTournament, tournsize=3)
# toolbox.register("map", map)

# view['evalOneMax'] = evalOneMax
# toolbox.register("map", view.map_sync)

toolbox.register("map", futures.map)


population = toolbox.population(n=300)
hof = tools.HallOfFame(5)

NGEN=40
for gen in xrange(NGEN):
    algorithms.eaSimple(population, toolbox, cxpb=0.5, mutpb=0.2, ngen=10, halloffame=hof, verbose=False)