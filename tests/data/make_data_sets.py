from os.path import dirname, join

import cobra
from cobra.io import read_sbml_model

from cameo import load_model
from cameo.io import sanitize_ids
from cameo.flux_analysis import phenotypic_phase_plane
from cameo.strain_design.deterministic.flux_variability_based import DifferentialFVA


config = cobra.Configuration()
TESTDIR = dirname(__file__)


def generate_ppp(model):
    tmp = model.copy()
    ppp = phenotypic_phase_plane(tmp, ["EX_o2_LPAREN_e_RPAREN_"])
    ppp.data_frame.to_csv(
        join(TESTDIR, "REFERENCE_PPP_o2_EcoliCore.csv"), float_format="%.3f"
    )

    tmp = model.copy()
    ppp2d = phenotypic_phase_plane(
        tmp, ["EX_o2_LPAREN_e_RPAREN_", "EX_glc_LPAREN_e_RPAREN_"]
    )
    ppp2d.data_frame.to_csv(
        join(TESTDIR, "REFERENCE_PPP_o2_glc_EcoliCore.csv"), float_format="%.3f"
    )

    tmp = model.copy()
    objective = tmp.add_boundary(tmp.metabolites.ac_c, type="demand")
    tmp.objective = objective
    ppp = phenotypic_phase_plane(tmp, ["EX_o2_LPAREN_e_RPAREN_"])
    ppp.data_frame.to_csv(
        join(TESTDIR, "REFERENCE_PPP_o2_EcoliCore_ac.csv"), float_format="%.3f"
    )


def generate_diff_fva(model):
    diff_fva = DifferentialFVA(model, model.reactions.EX_succ_lp_e_rp_, points=5)
    result = diff_fva.run()
    df = result.nth_panel(0)
    df.index.name = None
    df.to_csv(join(TESTDIR, "REFERENCE_DiffFVA1.csv"), float_format="%.3f")

    reference_model = model.copy()
    biomass_rxn = reference_model.reactions.Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2
    biomass_rxn.lower_bound = 0.3
    target = reference_model.reactions.EX_succ_lp_e_rp_
    target.lower_bound = 2
    result = DifferentialFVA(
        model, target, reference_model=reference_model, points=5
    ).run()
    df = result.nth_panel(0)
    df.index.name = None
    df.to_csv(join(TESTDIR, "REFERENCE_DiffFVA2.csv"), float_format="%.3f")


if __name__ == "__main__":
    core_model = load_model(join(TESTDIR, "EcoliCore.xml"), sanitize=False)
    core_model.solver = "glpk"
    generate_ppp(core_model)
    config.solver = "glpk"
    core_model = read_sbml_model(join(TESTDIR, "EcoliCore.xml"))
    sanitize_ids(core_model)
    generate_diff_fva(core_model)
