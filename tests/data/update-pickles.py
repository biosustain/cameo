import pickle

import cobra.test

from cameo import load_model

ijo = load_model('iJO1366.xml')
with open('iJO1366.pickle', 'wb') as out:
    pickle.dump(ijo, out)

salmonella = cobra.test.create_test_model('salmonella')
with open('salmonella.pickle', 'wb') as out:
    pickle.dump(salmonella, out)

