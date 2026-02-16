"""
Sharp (non-regularized) van Genuchten-Mualem constitutive law.

Hydraulic-only simplification of groundHeatRecovery's ConstitutiveSharp.
All methods accept pressure head h as CoefficientFunction and SoilParams,
returning CoefficientFunction.

Module-level functions ``K_of_theta_cf``, ``h_from_theta_cf`` operate on
theta (not h) and work with proxy functions and ``.Other()`` -- needed for
Godunov upwind and operator-splitting workflows.

NumPy helpers ``h_from_theta_np``, ``theta_from_h_np`` operate on DOF
arrays for fast element-level conversion between h and theta.
"""

import numpy as np
from ngsolve import CoefficientFunction as CF, IfPos

from .common import clamp_cf, pow_cf
from .materials import SoilParams

__all__ = [
    "ConstitutiveSharp",
    "K_of_theta_cf",
    "h_from_theta_cf",
    "h_from_theta_np",
    "theta_from_h_np",
]


class ConstitutiveSharp:
    """Sharp van Genuchten-Mualem constitutive law (no smoothing).

    All CoefficientFunction outputs use IfPos for sharp transitions at h = 0,
    distinguishing saturated (h >= 0) from unsaturated (h < 0) regimes.

    Parameters
    ----------
    eps_Se : float
        Clamp bounds for effective saturation (avoid 0 and 1 exactly).
    eps_a : float
        Guard for 1 - S_e^{1/m} term in Mualem formula.
    K_min : float
        Floor for hydraulic conductivity [m/s].
    """

    def __init__(self, *, eps_Se=1e-12, eps_a=1e-12, K_min=1e-16):
        self.eps_Se = eps_Se
        self.eps_a = eps_a
        self.K_min = K_min

    def Se(self, h: CF, p: SoilParams) -> CF:
        """Effective saturation S_e(h) [-].

        S_e = [1 + (alpha * max(-h, 0))^n]^(-m),  m = 1 - 1/n.
        Returns 1 for h >= 0 (saturated zone).

        Parameters
        ----------
        h : CF
            Pressure head [m].
        p : SoilParams
            Soil parameters.

        Returns
        -------
        CF
        """
        n = p.n
        m = 1.0 - 1.0 / n
        psi = IfPos(-h, -h, 0.0)  # suction = max(-h, 0)
        Se = (1.0 + (p.alpha * psi) ** n) ** (-m)
        return clamp_cf(Se, self.eps_Se, 1.0 - self.eps_Se)

    def theta(self, h: CF, p: SoilParams) -> CF:
        """Volumetric water content theta(h) [-].

        theta = theta_r + (theta_s - theta_r) * S_e(h).

        Parameters
        ----------
        h : CF
            Pressure head [m].
        p : SoilParams

        Returns
        -------
        CF
        """
        return p.theta_r + (p.theta_s - p.theta_r) * self.Se(h, p)

    def kr(self, h: CF, p: SoilParams) -> CF:
        """Mualem relative permeability k_r(h) [-].

        k_r = S_e^ell * [1 - (1 - S_e^{1/m})^m]^2.

        Parameters
        ----------
        h : CF
            Pressure head [m].
        p : SoilParams

        Returns
        -------
        CF
        """
        n = p.n
        m = 1.0 - 1.0 / n
        Se = self.Se(h, p)
        Se_clamped = clamp_cf(Se, self.eps_Se, 1.0 - self.eps_Se)
        # S_e^{1/m}
        t = pow_cf(Se_clamped, 1.0 / m)
        a = clamp_cf(1.0 - t, self.eps_a, 1.0)
        term = 1.0 - pow_cf(a, m)
        return pow_cf(Se, p.ell) * (term * term)

    def K(self, h: CF, p: SoilParams) -> CF:
        """Hydraulic conductivity K(h) = K_s * k_r(h) [m/s].

        Parameters
        ----------
        h : CF
            Pressure head [m].
        p : SoilParams

        Returns
        -------
        CF
        """
        K_raw = p.K_s * self.kr(h, p)
        return IfPos(K_raw - self.K_min, K_raw, self.K_min)

    def Kinv(self, h: CF, p: SoilParams) -> CF:
        """Inverse hydraulic conductivity K^{-1}(h) [s/m].

        For the mixed formulation: the mass matrix uses K^{-1}.

        Parameters
        ----------
        h : CF
            Pressure head [m].
        p : SoilParams

        Returns
        -------
        CF
        """
        return 1.0 / self.K(h, p)

    def C(self, h: CF, p: SoilParams) -> CF:
        """Specific moisture capacity C(h) = d(theta)/dh [1/m].

        Used in the Newton Jacobian for transient problems.
        C = (theta_s - theta_r) * dS_e/dh.

        For h >= 0: C = 0 (saturated).
        For h < 0: analytical derivative of van Genuchten S_e.

        Parameters
        ----------
        h : CF
            Pressure head [m].
        p : SoilParams

        Returns
        -------
        CF
        """
        n = p.n
        m = 1.0 - 1.0 / n
        psi = IfPos(-h, -h, 0.0)
        base = p.alpha * psi
        u = 1.0 + pow_cf(base, n)
        # dS_e/dpsi = -m * n * alpha * (alpha*psi)^{n-1} * (1 + (alpha*psi)^n)^{-(m+1)}
        dSe_dpsi = m * n * p.alpha * pow_cf(base, n - 1.0) * pow_cf(u, -(m + 1.0))
        # dtheta/dh = (theta_s - theta_r) * dS_e/dpsi * (-dpsi/dh)
        # dpsi/dh = -1 for h < 0, so dtheta/dh = +(theta_s - theta_r) * dSe_dpsi
        active = IfPos(-h, 1.0, 0.0)  # 1 if h < 0
        return (p.theta_s - p.theta_r) * dSe_dpsi * active


# ── Module-level theta-based functions ──────────────────────────────────


def K_of_theta_cf(theta_cf, soil, *, eps_Se=1e-12, eps_a=1e-12, K_min=1e-16):
    """Hydraulic conductivity K(theta) [m/s] as CoefficientFunction.

    Computes K directly from theta (no h intermediate), using the Mualem
    formula. Works with proxy functions and ``.Other()`` for DG upwind.

    Validated at machine precision vs NumPy in NB 18.

    Parameters
    ----------
    theta_cf : CF or ProxyFunction
        Volumetric water content [-].
    soil : SoilParams
        Soil parameters.
    eps_Se : float
        Clamp bounds for effective saturation.
    eps_a : float
        Guard for 1 - S_e^{1/m}.
    K_min : float
        Floor for K [m/s].

    Returns
    -------
    CF
    """
    m = 1.0 - 1.0 / soil.n
    Se = clamp_cf(
        (theta_cf - soil.theta_r) / (soil.theta_s - soil.theta_r),
        eps_Se, 1.0 - eps_Se,
    )
    t = pow_cf(Se, 1.0 / m)
    a = clamp_cf(1.0 - t, eps_a, 1.0)
    bracket = 1.0 - pow_cf(a, m)
    kr = pow_cf(Se, soil.ell) * bracket * bracket
    K_raw = soil.K_s * kr
    return IfPos(K_raw - K_min, K_raw, CF(K_min))


def h_from_theta_cf(theta_cf, soil, *, eps_Se=1e-12):
    """Inverse van Genuchten: h(theta) [m] as CoefficientFunction.

    h = -(1/alpha) * (Se^{-1/m} - 1)^{1/n}.
    Returns 0 for theta >= theta_s (saturated).

    Parameters
    ----------
    theta_cf : CF or ProxyFunction
        Volumetric water content [-].
    soil : SoilParams
        Soil parameters.
    eps_Se : float
        Clamp bounds for effective saturation.

    Returns
    -------
    CF
        Pressure head h [m] (non-positive).
    """
    m = 1.0 - 1.0 / soil.n
    Se = clamp_cf(
        (theta_cf - soil.theta_r) / (soil.theta_s - soil.theta_r),
        eps_Se, 1.0 - eps_Se,
    )
    Se_inv_m = pow_cf(Se, -1.0 / m)
    bracket = clamp_cf(Se_inv_m - 1.0, 1e-20, 1e20)
    psi = pow_cf(bracket, 1.0 / soil.n) / soil.alpha
    return IfPos(theta_cf - (soil.theta_s - 1e-10), CF(0.0), -psi)


def h_from_theta_np(theta_np, soil):
    """Inverse van Genuchten: h(theta) [m], NumPy version for DOF arrays.

    Parameters
    ----------
    theta_np : ndarray
        Water content values [-].
    soil : SoilParams

    Returns
    -------
    ndarray
        Pressure head h [m] (non-positive).
    """
    m = 1.0 - 1.0 / soil.n
    Se = np.clip(
        (theta_np - soil.theta_r) / (soil.theta_s - soil.theta_r),
        1e-12, 1.0 - 1e-12,
    )
    psi = (Se ** (-1.0 / m) - 1.0) ** (1.0 / soil.n) / soil.alpha
    return np.where(theta_np >= soil.theta_s - 1e-10, 0.0, -psi)


def theta_from_h_np(h_np, soil):
    """Van Genuchten theta(h) [-], NumPy version for DOF arrays.

    Parameters
    ----------
    h_np : ndarray
        Pressure head values [m].
    soil : SoilParams

    Returns
    -------
    ndarray
        Volumetric water content [-].
    """
    m = 1.0 - 1.0 / soil.n
    psi = np.maximum(-h_np, 0.0)
    Se = (1.0 + (soil.alpha * psi) ** soil.n) ** (-m)
    return soil.theta_r + (soil.theta_s - soil.theta_r) * Se
