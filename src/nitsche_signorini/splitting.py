"""
Lie-Trotter operator splitting for Richards equation.

Couples Godunov upwind convection (explicit, no K^{-1}) with mixed FEM
diffusion (implicit, HDiv x L2). Building blocks only; time loops live
in notebooks. Validated in NB 17 (NumPy), NB 18 (NGSolve DG), NB 19
(coupled splitting).

Sign convention for the split:
- Convective step: dtheta/dt + div(K(theta)*e_z) = 0  (Godunov upwind)
- Diffusion step:  K^{-1} q + grad h = 0,  div q + (1/dt)(theta(h) - theta*) = 0
  where q = -K grad h is the diffusive flux only (no gravity).
- Total Darcy flux: q_total = q_diff + q_conv = q - K*e_z
"""

import numpy as np
from ngsolve import (
    BilinearForm,
    CF,
    GridFunction,
    IfPos,
    InnerProduct,
    L2,
    HDiv,
    Mesh,
    div,
    dx,
    ds,
    specialcf,
    sqrt,
)
from netgen.geom2d import SplineGeometry

from .common import ez_cf
from .constitutive import ConstitutiveSharp, K_of_theta_cf
from .materials import SoilParams

__all__ = [
    "compute_cfl_dt",
    "build_mesh_and_spaces",
    "build_godunov_flux",
    "godunov_subcycle",
    "build_diffusion_form",
    "build_diffusion_form_signorini",
    "build_monolithic_form_signorini",
]

cons = ConstitutiveSharp()


def compute_cfl_dt(soil, dz, cfl=0.5):
    """CFL time step for Godunov convection based on max|dK/dtheta|.

    Parameters
    ----------
    soil : SoilParams
        Soil parameters.
    dz : float
        Mesh size [m].
    cfl : float
        CFL number (default 0.5).

    Returns
    -------
    dt : float
        CFL time step [s].
    max_speed : float
        Maximum characteristic speed |dK/dtheta| [m/s].
    """
    th = np.linspace(soil.theta_r + 1e-10, soil.theta_s, 10000)
    m = 1.0 - 1.0 / soil.n
    Se = np.clip(
        (th - soil.theta_r) / (soil.theta_s - soil.theta_r), 1e-12, 1 - 1e-12
    )
    t = Se ** (1.0 / m)
    a = np.clip(1.0 - t, 1e-12, 1.0)
    kr = Se ** soil.ell * (1.0 - a ** m) ** 2
    K_vals = np.maximum(soil.K_s * kr, 1e-16)
    dK = np.gradient(K_vals, th)
    max_speed = np.max(np.abs(dK))
    return cfl * dz / max_speed, max_speed


def build_mesh_and_spaces(H, soil, W=0.1, maxh=0.05):
    """Build 2D column mesh and FE spaces for operator splitting.

    Boundary labels: bottom, right, top, left.
    No-flow Dirichlet on left/right for HDiv.

    Parameters
    ----------
    H : float
        Column height [m].
    soil : SoilParams
        Soil parameters (unused, kept for interface consistency).
    W : float
        Column width [m] (default 0.1).
    maxh : float
        Maximum element size [m] (default 0.05).

    Returns
    -------
    mesh : Mesh
    Q_dg : L2
        Godunov space (P0, dgjumps=True).
    fes : ProductSpace
        HDiv(RT0) x L2(P0) for mixed diffusion.
    """
    geo = SplineGeometry()
    pts = [geo.AppendPoint(px, py) for px, py in
           [(0, 0), (W, 0), (W, H), (0, H)]]
    geo.Append(["line", pts[0], pts[1]], bc="bottom")
    geo.Append(["line", pts[1], pts[2]], bc="right")
    geo.Append(["line", pts[2], pts[3]], bc="top")
    geo.Append(["line", pts[3], pts[0]], bc="left")
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))

    Q_dg = L2(mesh, order=0, dgjumps=True)
    V = HDiv(mesh, order=0, RT=True, dirichlet="left|right")
    Q = L2(mesh, order=0)
    fes = V * Q

    assert Q_dg.ndof == Q.ndof, "DOF count mismatch Q_dg vs Q_mix"
    return mesh, Q_dg, fes


def build_godunov_flux(mesh, Q_dg, soil, theta_top=None):
    """Build upwind flux BilinearForm for Godunov convection.

    Flux: F(theta) = -K(theta) e_z. Upwind selection via IfPos on e_z * n.
    Boundary ghosts: theta_s at top/bottom (saturated), copy at sides.

    Parameters
    ----------
    mesh : Mesh
    Q_dg : L2
        Godunov space (P0, dgjumps=True).
    soil : SoilParams
    theta_top : CF or Parameter or None
        Top boundary ghost value. Default: soil.theta_s (seepage mode).
        Pass a Parameter for dynamic Signorini ghost.

    Returns
    -------
    a_flux : BilinearForm
        Non-assembled upwind flux form.
    """
    U, V_test = Q_dg.TnT()
    n = specialcf.normal(mesh.dim)
    ez = ez_cf(mesh)
    ez_n = InnerProduct(ez, n)

    if theta_top is None:
        theta_top = soil.theta_s

    theta_bnd = mesh.BoundaryCF({
        "top": theta_top,
        "bottom": soil.theta_s,
        "left": U, "right": U,
    })
    Uhat = U.Other(bnd=theta_bnd)

    # Upwind: ez_n > 0 => face above element => use neighbor
    theta_up = IfPos(ez_n, Uhat, U)
    K_up = K_of_theta_cf(theta_up, soil)
    Fhatn = (-K_up * ez_n).Compile()

    a_flux = BilinearForm(Q_dg, nonassemble=True)
    a_flux += Fhatn * V_test * dx(element_boundary=True)
    return a_flux


def godunov_subcycle(gf_theta, a_flux, Q_dg, soil, dt_macro, dt_cfl):
    """Run Godunov CFL sub-steps within one macro time step.

    Explicit Euler with in-place clipping to [theta_r, theta_s].

    Parameters
    ----------
    gf_theta : GridFunction
        Solution on Q_dg, updated in place.
    a_flux : BilinearForm
        Non-assembled upwind flux.
    Q_dg : L2
        Godunov space.
    soil : SoilParams
    dt_macro : float
        Macro time step [s].
    dt_cfl : float
        CFL time step [s].

    Returns
    -------
    n_sub : int
        Number of sub-steps taken.
    """
    n_sub = max(1, int(np.ceil(dt_macro / dt_cfl)))
    dt_sub = dt_macro / n_sub
    res = gf_theta.vec.CreateVector()

    for _ in range(n_sub):
        a_flux.Apply(gf_theta.vec, res)
        Q_dg.SolveM(rho=CF(1), vec=res)
        gf_theta.vec.data -= dt_sub * res
        theta_np = gf_theta.vec.FV().NumPy()
        np.clip(theta_np, soil.theta_r, soil.theta_s, out=theta_np)

    return n_sub


def build_diffusion_form(fes, soil, dt_param, gf_theta_star):
    """Diffusion-only mixed weak form (no gravity, no Signorini).

    Natural BC (omitting boundary integrals) enforces h=0 on all boundaries.
    Gravity is handled by the Godunov convective step.

    Parameters
    ----------
    fes : ProductSpace
        HDiv(RT0) x L2(P0).
    soil : SoilParams
    dt_param : Parameter
        Time step [s].
    gf_theta_star : GridFunction
        Godunov output theta* on L2(P0) component of fes.

    Returns
    -------
    R : BilinearForm
    """
    (q, h), (tau, v) = fes.TnT()

    R = BilinearForm(fes)
    # Darcy (diffusion only, NO gravity)
    R += cons.Kinv(h, soil) * InnerProduct(q, tau) * dx
    R += -h * div(tau) * dx

    # Mass conservation with storage
    R += div(q) * v * dx
    R += (1.0 / dt_param) * cons.theta(h, soil) * v * dx
    R += -(1.0 / dt_param) * gf_theta_star * v * dx

    return R


def _smooth_pos(x, eps):
    """Smooth approximation of [x]+ = max(x, 0).

    Uses (x + sqrt(x^2 + eps^2)) / 2. Equals [x]+ up to O(eps).
    Derivative is continuous everywhere (no kink at x = 0).
    """
    return (x + sqrt(x * x + eps * eps)) / 2.0


def build_diffusion_form_signorini(fes, mesh, soil, dt_param,
                                    gf_theta_star, p_rain=0.0):
    """Diffusion-only mixed weak form WITH Nitsche-Signorini on top.

    The Signorini switching function uses the total flux (diffusive q plus
    convective -K*e_z) to determine wet/dry mode on the top boundary.

    Total Darcy flux: q_total = q - K(h)*e_z, where q is the diffusive flux
    from the mixed system and -K*e_z is the convective (gravity) flux handled
    by the Godunov step.

    The complementarity [Q_h - gamma*h]_+ is regularized with a smooth
    approximation (x + sqrt(x^2 + eps^2))/2, removing the Jacobian
    discontinuity at the seepage/infiltration transition. The regularization
    parameter eps is set equal to gamma so that the smooth residual, after
    multiplication by 1/gamma, remains O(1) and does not produce spurious
    boundary contributions.

    Parameters
    ----------
    fes : ProductSpace
        HDiv(RT0) x L2(P0).
    mesh : Mesh
    soil : SoilParams
    dt_param : Parameter
        Time step [s].
    gf_theta_star : GridFunction
        Godunov output theta* on L2(P0) component of fes.
    p_rain : float
        Prescribed precipitation rate [m/s] (default 0).

    Returns
    -------
    R : BilinearForm
    """
    (q, h), (tau, v) = fes.TnT()
    nn = specialcf.normal(mesh.dim)
    h_mesh = specialcf.mesh_size
    ez = ez_cf(mesh)

    R = BilinearForm(fes)
    # Darcy (diffusion only, NO gravity)
    R += cons.Kinv(h, soil) * InnerProduct(q, tau) * dx
    R += -h * div(tau) * dx

    # Storage
    R += div(q) * v * dx
    R += (1.0 / dt_param) * cons.theta(h, soil) * v * dx
    R += -(1.0 / dt_param) * gf_theta_star * v * dx

    # Nitsche-Signorini on top boundary
    gamma_0 = 1e-10
    gamma = gamma_0 * h_mesh
    dS_top = ds(skeleton=True, definedon=mesh.Boundaries("top"))

    p_cf = -p_rain * ez  # prescribed precipitation flux (downward)
    K_h = cons.K(h, soil)
    # Total flux = diffusive q + convective (-K*e_z) = q - K*e_z
    q_total = q - K_h * ez
    Q_h = (p_cf - q_total) * nn  # flux excess
    nitsche_arg = Q_h - gamma * h
    # eps = gamma: smooth residual O(gamma), divided by gamma -> O(1)
    R += (-(1.0 / gamma)
          * _smooth_pos(nitsche_arg, gamma)
          * (tau * nn) * dS_top)

    return R


def build_monolithic_form_signorini(fes, mesh, soil, dt_param,
                                     gf_h_old, p_rain=0.0):
    """Monolithic mixed weak form with gravity + Nitsche-Signorini.

    Unlike the split form, gravity is implicit in the flux equation and
    the storage RHS uses theta(h_old) from the previous time step.
    Used for the base method (clay) where operator splitting is not needed.

    Adapts the NB 03 weak form into a reusable function with the same
    interface conventions as ``build_diffusion_form_signorini``.

    Parameters
    ----------
    fes : ProductSpace
        HDiv(RT0) x L2(P0).
    mesh : Mesh
    soil : SoilParams
    dt_param : Parameter
        Time step [s].
    gf_h_old : GridFunction
        Pressure head from previous time step on L2(P0).
    p_rain : float or Parameter
        Precipitation rate [m/s].

    Returns
    -------
    R : BilinearForm
    """
    (q, h), (tau, v) = fes.TnT()
    nn = specialcf.normal(mesh.dim)
    h_mesh = specialcf.mesh_size
    ez = ez_cf(mesh)

    R = BilinearForm(fes)
    # Darcy WITH gravity: K^{-1} q + grad(h) + e_z = 0
    R += cons.Kinv(h, soil) * InnerProduct(q, tau) * dx
    R += -h * div(tau) * dx
    R += InnerProduct(ez, tau) * dx

    # Storage: (1/dt)[theta(h^{n+1}) - theta(h^n)]
    R += div(q) * v * dx
    R += (1.0 / dt_param) * cons.theta(h, soil) * v * dx
    R += -(1.0 / dt_param) * cons.theta(gf_h_old, soil) * v * dx

    # Nitsche-Signorini on top (q is total flux, not split)
    gamma_0 = 1e-10
    gamma = gamma_0 * h_mesh
    dS_top = ds(skeleton=True, definedon=mesh.Boundaries("top"))
    p_cf = -p_rain * ez
    Q_h = (p_cf - q) * nn
    nitsche_arg = Q_h - gamma * h
    R += -(1.0 / gamma) * IfPos(nitsche_arg, nitsche_arg, 0.0) * (tau * nn) * dS_top

    return R
