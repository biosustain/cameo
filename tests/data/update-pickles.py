import pickle

import cobra.test
import optlang

from cameo import load_model

ijo = load_model('iJO1366.xml', solver_interface=optlang.glpk_interface)
with open('iJO1366.pickle', 'wb') as out:
    pickle.dump(ijo, out, protocol=2)

salmonella = cobra.test.create_test_model('salmonella')
salmonella.solver = 'glpk'
with open('salmonella.pickle', 'wb') as out:
    pickle.dump(salmonella, out, protocol=2)
