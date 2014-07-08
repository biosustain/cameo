[![Build Status](https://travis-ci.org/biosustain/cameo.svg?branch=devel)](https://travis-ci.org/biosustain/cameo)

##cameo - computer assisted metabolic engineering & optimization *....UNDER CONSTRUCTION....*

### Vision
Provide a high-level python library to aid the strain design efforts of the CFB iLoop design group. The library provides a modular architecture that enables the efficient construction of custom analysis workflows.

### Design

**No more copying** making copies of models is expensive, stick to the [_command pattern_](http://en.wikipedia.org/wiki/Command_pattern) using TimeMachine (util.py)
**Database**: no more flatfiles ...
**Caching**: optimization results are cached
**Parallelization**: parallelize algorithms
**Long-running jobs**: hide paralellization
**Wrap third party tools**: for example [fast-tFVA](http://bioinformatics.oxfordjournals.org/content/29/7/903)

### Dependencies
This library dependes on

- [cobrapy](https://github.com/opencobra/cobrapy) for constraint-based modeling
- [optlang](https://github.com/biosustain/optlang) for heuristic optimization and mathematical programming

Computationally heavy methods have been parallelized and can be run on a clusters using the IPython parallelization framework (see example and documetnation for more details). The default fallback is python's multiprocessing library.

