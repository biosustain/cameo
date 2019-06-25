from os.path import abspath, dirname, join

import cobra.test
import pytest

from cameo import load_model
from cameo.config import solvers


cameo_directory = abspath(join(dirname(abspath(__file__)), ".."))
cameo_location = abspath(join(cameo_directory, ".."))
data_dir = join(cameo_directory, "tests", "data", "")


@pytest.fixture(scope="session")
def data_directory():
    return data_dir


@pytest.fixture(scope='session')
def model(data_directory):
    m = load_model(join(data_directory, 'EcoliCore.xml'))
    m.solver = "glpk"
    return m


@pytest.fixture(scope="session")
def salmonella():
    return cobra.test.create_test_model('salmonella')


@pytest.fixture(scope="session", params=list(solvers))
def ijo1366(request, data_directory):
    ijo = load_model(join(data_directory, 'iJO1366.xml'), sanitize=False)
    ijo.solver = request.param
    return ijo


@pytest.fixture(scope='session')
def iaf1260(data_directory):
    return load_model(join(data_directory, 'iAF1260.xml'))


@pytest.fixture(scope="session")
def universal_model(data_directory):
    universal = load_model(join(data_directory, 'iJO1366.xml'), sanitize=False)
    universal.remove_reactions(universal.boundary)
    return universal


@pytest.fixture(scope="session", params=list(solvers))
def imm904(request, data_directory):
    imm = load_model(join(data_directory, 'iMM904.xml'))
    imm.solver = request.param
    return imm.copy()


# FIXME: should be possible to scope at least to class
@pytest.fixture(scope="function", params=list(solvers))
def toy_model(request, data_directory):
    toy = load_model(join(data_directory, "toy_model_Papin_2003.xml"))
    toy.solver = request.param
    return toy


# FIXME: should be possible to scope at least to class
@pytest.fixture(scope="function", params=list(solvers))
def core_model(request, data_directory):
    ecoli_core = load_model(join(data_directory, 'EcoliCore.xml'), sanitize=False)
    ecoli_core.solver = request.param
    return ecoli_core
