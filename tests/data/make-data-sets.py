from cameo import load_model
from cameo.flux_analysis import phenotypic_phase_plane
import os

TESTDIR = os.path.dirname(__file__)
CORE_MODEL = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'), sanitize=False)

model = CORE_MODEL.copy()
model.solver = 'glpk'
ppp = phenotypic_phase_plane(model, ['EX_o2_LPAREN_e_RPAREN_'])
ppp.data_frame.to_csv(os.path.join(TESTDIR, 'data/REFERENCE_PPP_o2_EcoliCore.csv'))

model = CORE_MODEL.copy()
ppp2d = phenotypic_phase_plane(model, ['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'])
ppp2d.data_frame.to_csv(os.path.join(TESTDIR, 'data/REFERENCE_PPP_o2_glc_EcoliCore.csv'))

model = CORE_MODEL.copy()
objective = model.add_demand(model.metabolites.ac_c)
model.change_objective(objective)
ppp = phenotypic_phase_plane(model, ['EX_o2_LPAREN_e_RPAREN_'])
ppp.data_frame.to_csv(os.path.join(TESTDIR, 'data/REFERENCE_PPP_o2_EcoliCore_ac.csv'))
