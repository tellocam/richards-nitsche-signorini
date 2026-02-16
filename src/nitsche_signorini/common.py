"""
Shared CoefficientFunction helpers for constitutive law implementations.

Ported and simplified from groundHeatRecovery common.py (hydraulic only).
"""

from ngsolve import CoefficientFunction as CF, IfPos, exp, log, y, z

__all__ = [
    "clamp_cf",
    "pow_cf",
    "ez_cf",
    "vertical_coord_cf",
    "neg_part",
    "pos_part",
]


def clamp_cf(val, lo, hi):
    """min(max(val, lo), hi) via IfPos.

    Parameters
    ----------
    val : CF or float
        Value to clamp.
    lo, hi : float
        Lower and upper bounds.

    Returns
    -------
    CF
    """
    return IfPos(val - lo, IfPos(hi - val, val, hi), lo)


def pow_cf(base, exponent, *, base_lo=1e-30, base_hi=1e30):
    """base^exponent for CoefficientFunctions with variable exponent.

    Uses exp(exponent * log(base)) with clamped base for safety.

    Parameters
    ----------
    base : CF or float
        Must be positive.
    exponent : CF or float
        Exponent (can be CF).
    base_lo, base_hi : float
        Clamp bounds for base before taking log.

    Returns
    -------
    CF
    """
    b = clamp_cf(base, base_lo, base_hi)
    return exp(exponent * log(b))


def ez_cf(mesh):
    """Vertical unit vector e_z as CoefficientFunction.

    Parameters
    ----------
    mesh : ngsolve.Mesh

    Returns
    -------
    CF
        (0, 1) in 2D, (0, 0, 1) in 3D.
    """
    if mesh.dim == 2:
        return CF((0, 1))
    if mesh.dim == 3:
        return CF((0, 0, 1))
    raise ValueError(f"Unsupported mesh dimension: {mesh.dim}")


def vertical_coord_cf(mesh):
    """Vertical coordinate CF (y in 2D, z in 3D).

    Parameters
    ----------
    mesh : ngsolve.Mesh

    Returns
    -------
    CF
    """
    return y if mesh.dim == 2 else z


def neg_part(x):
    """Negative part [x]_- = min(x, 0).

    Parameters
    ----------
    x : CF or float

    Returns
    -------
    CF
    """
    return IfPos(x, 0, x)


def pos_part(x):
    """Positive part [x]_+ = max(x, 0).

    Parameters
    ----------
    x : CF or float

    Returns
    -------
    CF
    """
    return IfPos(x, x, 0)
