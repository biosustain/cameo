from cameo.basics import cycle_free_fva

from cameo import load_model
from cameo.basics import cycle_free_fva, _flux_variability_analysis

ecoli_core = load_model('../tests/data/EcoliCore.xml')
fva_sol = _flux_variability_analysis(ecoli_core)
cycle_free_fva_sol = cycle_free_fva(ecoli_core, sloppy=True)

# print fva_sol['FRD7']
# print fva_sol['SUCDi']
print cycle_free_fva_sol['FRD7']
print cycle_free_fva_sol['SUCDi']
