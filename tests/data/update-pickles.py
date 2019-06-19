import pickle

import cobra
from cobra.test import create_test_model
from optlang import glpk_interface

from cameo import load_model


config = cobra.Configuration()
config.solver = "glpk"


ijo = load_model('iJO1366.xml', glpk_interface)
with open('iJO1366.pickle', 'wb') as out:
    pickle.dump(ijo, out, protocol=2)

salmonella = create_test_model('salmonella')
with open('salmonella.pickle', 'wb') as out:
    pickle.dump(salmonella, out, protocol=2)
