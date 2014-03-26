from cameo.config import default_view
from cameo.io import load_model
from cameo.parallel import SequentialView, MultiprocessingView
from cameo.basics import _flux_variability_analysis, flux_variability_analysis

model = load_model('../tests/data/EcoliCore.xml')
# print _flux_variability_analysis(model)

print default_view

try:
    print flux_variability_analysis(model, view=SequentialView())
except:
    pass

try:
    print flux_variability_analysis(model, view=MultiprocessingView())
except Exception, e:
    print e
    print "Multiprocessing didn't work."

try:
    print flux_variability_analysis(model)
except Exception, e:
    print e
    print "Ipython parallel didn't work"

