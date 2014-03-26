# class Algorithm(object):
#     """docstring for Algorithm"""
#     def __init__(self):
#         super(Algorithm, self).__init__()

# class FluxGaps(Algorithm):
#     """This method implements"""
#     def __init__(self, arg):
#         super(FluxGaps, self).__init__()
#         self.arg = arg

#     def run(self):

#         return flux_gaps

# class Redirector(object):
#     """docstring for Redirector"""
#     def __init__(self, arg):
#         super(Redirector, self).__init__()
#         self.arg = arg
#         

# class GenericStrainDesignMethod(object):
    
#     default_evolver = inspyred.ec.GA
#     default_observer = inspyred.ec.observers.default_observer
    
#     def __init__(self, model=None, evolver=None):
#         self.model = model
#         if evolver is None:
#             # Generate the default EvolutionaryComputation object
#             self.evolver = self.default_evolver
#         else:
#             self.evolver = evolver

#     @staticmethod
#     def _binary_to_reactions(model, candidate):
#         return [reaction for reaction, bool in zip(copy.reactions, candidate) if elem == 1]
    
# #     @inspyred.ec.generators.diversify
#     def generate(self, random, args):
# #          print [random.choice([0, 1]) for _ in range(len(self.model.reactions))][0]
#         return [random.choice([0, 1]) for _ in range(len(self.model.reactions))]
    
#     # TODO: The following decorator should be configurable, e.g., inspyred.ec.evaluators.parallel_evaulation_mp
#     # http://stackoverflow.com/questions/642762/is-it-possible-to-replace-a-function-method-decorator-at-runtime-python
#     # might do the trick
#     @inspyred.ec.evaluators.evaluator
#     def evaluate(self, candidate):
#         return sum(candidate)
# #         copy = self.model.copy()
# #         reactions_to_knock_out = self._binary_to_reactions(copy, candidate)
# #         copy.remove_reactions(reactions_to_knock_out) # that's happening inplace
# #         return copy.optimize().solution.f
    
#     def run(self):
#         rand = Random()
# #         rand.seed(int(time()))
#         evo = self.evolver(rand)
#         final_population = evo.evolve(generator=self.generate,
#                             evaluator=self.evaluate,
#                             pop_size=10, # TODO: should be configurable, also could probably be moved out of the function call
#                             maxizmie=True,
#                             bounder=inspyred.ec.Bounder(0, 1),
#                             observer=self.default_observer
#                             )
#         return [self._binary_to_reactions(self.model, candidate) for candidate in final_population]
