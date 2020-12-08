# Copyright 2019 Novo Nordisk Foundation Center for Biosustainability, DTU.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from cameo.flux_analysis.analysis import phenotypic_phase_plane
from cameo.visualization.plotting.with_plotly import PlotlyPlotter

plotter = PlotlyPlotter()


def test_ppp_plotly(model):
    """Test if at least it doesn't raise an exception."""
    production_envelope = phenotypic_phase_plane(
        model,
        variables=[model.reactions.Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2],
        objective=model.metabolites.succ_e,
    )
    # pass a grid argument so that it doesn't open a browser tab
    production_envelope.plot(plotter, height=400, grid=[])
