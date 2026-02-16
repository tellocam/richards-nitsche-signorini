"""
Physical material properties for van Genuchten-Mualem hydraulic model.

Hydraulic-only simplification of groundHeatRecovery materials module.
"""

from dataclasses import dataclass

__all__ = [
    "SoilParams",
    "SANDY_LOAM",
    "CLAY",
    "SAND",
]


@dataclass
class SoilParams:
    """Hydraulic soil parameters for van Genuchten-Mualem model.

    Parameters
    ----------
    theta_s : float
        Saturated water content [-].
    theta_r : float
        Residual water content [-].
    alpha : float
        van Genuchten parameter [1/m].
    n : float
        van Genuchten parameter [-], must be > 1.
    ell : float
        Mualem tortuosity exponent [-], typically 0.5.
    K_s : float
        Saturated hydraulic conductivity [m/s].
    """

    theta_s: float
    theta_r: float
    alpha: float  # [1/m]
    n: float  # [-], n > 1
    ell: float  # [-], Mualem exponent
    K_s: float  # [m/s]


# ── Presets ──────────────────────────────────────────────────────────────

# Values from Gatti et al. 2024 (Table in Section 4, p/K_S tests)
# and groundHeatRecovery (converting kint to K_s = rho*g*kint/mu)

SANDY_LOAM = SoilParams(
    theta_s=0.40,
    theta_r=0.00,
    alpha=2.0,       # [1/m]
    n=3.0,
    ell=0.5,
    K_s=1e-4,        # [m/s]
)

CLAY = SoilParams(
    theta_s=0.40,
    theta_r=0.04,
    alpha=0.2,        # [1/m]
    n=1.5,
    ell=0.5,
    K_s=1e-6,         # [m/s]
)

SAND = SoilParams(
    theta_s=0.40,
    theta_r=0.00,
    alpha=2.0,         # [1/m]
    n=3.0,
    ell=0.5,
    K_s=1e-4,          # [m/s]
)
