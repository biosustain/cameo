## cameo - computer assisted metabolic engineering & optimization

[![Build Status](https://travis-ci.org/biosustain/cameo.svg?branch=devel)](https://travis-ci.org/biosustain/cameo)
[![Coverage Status](https://coveralls.io/repos/biosustain/cameo/badge.png?branch=devel)](https://coveralls.io/r/biosustain/cameo?branch=devel)
[![Gitter chat](https://badges.gitter.im/biosustain/cameo.png)](https://gitter.im/biosustain/cameo)

### Vision
Provide a high-level python library to aid _in silico_ strain design process in metabolic engineering projects. The library provides a modular architecture that enables the efficient construction of custom analysis workflows.

### Dependencies
This library dependes on

- [cobrapy](https://github.com/opencobra/cobrapy) for constraint-based modeling
- [optlang](https://github.com/biosustain/optlang) for heuristic optimization and mathematical programming

Computationally heavy methods have been parallelized and can be run on a clusters using the IPython parallelization framework (see example and documetnation for more details). The default fallback is python's multiprocessing library.
