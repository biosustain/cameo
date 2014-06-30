import optlang.interface


class ModelSolveError(Exception):
    def __init__(self, message):
        super(ModelSolveError, self).__init__(message)


class ModelInfeasible(ModelSolveError):
    pass


class ModelUnbounded(ModelSolveError):
    pass


class ModelFeasibleButNotOptimal(ModelSolveError):
    pass


class ModelUndefinedSolution(ModelSolveError):
    pass


_OPTLANG_TO_EXCEPTIONS_DICT = dict((
    (optlang.interface.INFEASIBLE, ModelInfeasible), (optlang.interface.UNBOUNDED, ModelUnbounded),
    (optlang.interface.FEASIBLE, ModelFeasibleButNotOptimal), (optlang.interface.UNDEFINED, ModelUndefinedSolution)))
