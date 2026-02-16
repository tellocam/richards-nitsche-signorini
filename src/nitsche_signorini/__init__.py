"""
HDiv x L2 mixed FEM for Richards equation with Nitsche-Signorini BCs.
"""

from .common import clamp_cf, ez_cf, neg_part, pow_cf
from .constitutive import (
    ConstitutiveSharp,
    K_of_theta_cf,
    h_from_theta_cf,
    h_from_theta_np,
    theta_from_h_np,
)
from .materials import CLAY, SAND, SANDY_LOAM, SoilParams
from .splitting import (
    build_diffusion_form,
    build_diffusion_form_signorini,
    build_godunov_flux,
    build_mesh_and_spaces,
    build_monolithic_form_signorini,
    compute_cfl_dt,
    godunov_subcycle,
)

__all__ = [
    "SoilParams",
    "SANDY_LOAM",
    "CLAY",
    "SAND",
    "ConstitutiveSharp",
    "clamp_cf",
    "pow_cf",
    "ez_cf",
    "neg_part",
    "K_of_theta_cf",
    "h_from_theta_cf",
    "h_from_theta_np",
    "theta_from_h_np",
    "compute_cfl_dt",
    "build_mesh_and_spaces",
    "build_godunov_flux",
    "godunov_subcycle",
    "build_diffusion_form",
    "build_diffusion_form_signorini",
    "build_monolithic_form_signorini",
]
