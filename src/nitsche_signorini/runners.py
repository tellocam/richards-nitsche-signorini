"""
Time-stepping runners for reproducibility scripts.

Consolidated from research notebooks 03, 22--27. Each runner sets up
mesh, spaces, forms, and runs the time loop, returning profile snapshots
and diagnostic info.
"""

import numpy as np
from ngsolve import (
    BilinearForm,
    CF,
    GridFunction,
    Integrate,
    L2,
    HDiv,
    Mesh,
    Parameter,
    TaskManager,
    cos,
    div,
    sin,
    sqrt,
    y,
    dx,
    ds,
    specialcf,
    InnerProduct,
    IfPos,
)
from ngsolve.solvers import Newton
from netgen.geom2d import SplineGeometry
from scipy.integrate import solve_ivp

from .common import ez_cf
from .constitutive import ConstitutiveSharp, h_from_theta_np, theta_from_h_np
from .materials import SoilParams
from .splitting import (
    build_diffusion_form_signorini,
    build_godunov_flux,
    build_mesh_and_spaces,
    build_monolithic_form_signorini,
    compute_cfl_dt,
    godunov_subcycle,
)

__all__ = [
    "run_monolithic",
    "run_splitting",
    "run_splitting_ramp",
    "run_godunov_only",
    "find_front_position",
    "solve_manufactured",
    "run_temporal_convergence",
    "solve_rh_ode",
    "build_K_table",
    "theta_from_K_table",
    "make_piecewise_linear",
    "run_splitting_temporal_convergence",
]

cons = ConstitutiveSharp()


# ── Helpers ──────────────────────────────────────────────────────────────


def make_piecewise_linear(breakpoints):
    """Build a piecewise-linear function from (t, p) breakpoints.

    Uses np.interp: constant extrapolation beyond endpoints.

    Parameters
    ----------
    breakpoints : list of (float, float)
        [(t0, p0), (t1, p1), ...], sorted by t.

    Returns
    -------
    callable
        p(t) -> float.
    """
    ts = np.array([bp[0] for bp in breakpoints])
    ps = np.array([bp[1] for bp in breakpoints])
    return lambda t: float(np.interp(t, ts, ps))


def build_K_table(soil, n_pts=10000):
    """Precompute K(theta) lookup table for inverse Mualem.

    Parameters
    ----------
    soil : SoilParams
    n_pts : int

    Returns
    -------
    theta_tab, K_tab : ndarray
        K_tab is monotonically increasing.
    """
    theta_tab = np.linspace(
        soil.theta_r + 1e-10, soil.theta_s - 1e-10, n_pts)
    m = 1.0 - 1.0 / soil.n
    Se = np.clip(
        (theta_tab - soil.theta_r) / (soil.theta_s - soil.theta_r),
        1e-12, 1 - 1e-12)
    t = Se ** (1.0 / m)
    a = np.clip(1.0 - t, 1e-12, 1.0)
    kr = Se ** soil.ell * (1.0 - a ** m) ** 2
    K_tab = soil.K_s * kr
    return theta_tab, K_tab


def theta_from_K_table(K_target, soil, theta_tab, K_tab):
    """Inverse Mualem: find theta such that K(theta) = K_target.

    Uses precomputed table with np.interp. Cannot fail.

    Parameters
    ----------
    K_target : float
        Target hydraulic conductivity [m/s].
    soil : SoilParams
    theta_tab, K_tab : ndarray
        From build_K_table().

    Returns
    -------
    float
    """
    if K_target >= soil.K_s:
        return soil.theta_s
    if K_target <= 0:
        return soil.theta_r
    return float(np.interp(K_target, K_tab, theta_tab))


def find_front_position(th_prof, z_profile, soil):
    """Find wetting front position from theta profile.

    The front descends from z=H (top). Above the front, theta ~ theta_s.
    Below, theta follows the IC. The front is at the maximum positive
    gradient in the ascending z_profile.

    Parameters
    ----------
    th_prof : ndarray
    z_profile : ndarray
        Ascending vertical coordinates.
    soil : SoilParams

    Returns
    -------
    float
        Front position z [m].
    """
    if np.min(th_prof) > 0.8 * soil.theta_s:
        return 0.0

    dtheta = np.diff(th_prof)
    dz = np.diff(z_profile)
    grad = dtheta / dz

    max_grad_idx = np.argmax(grad)
    max_grad = grad[max_grad_idx]

    theta_range = soil.theta_s - soil.theta_r
    H = z_profile[-1] - z_profile[0]
    if max_grad < theta_range / H:
        return z_profile[-1]

    return 0.5 * (z_profile[max_grad_idx] + z_profile[max_grad_idx + 1])


def solve_rh_ode(soil, H, t_span, n_eval=500):
    """Solve the variable-speed Rankine-Hugoniot front ODE.

    dz/dt = -s_local(z), where s_local = [K_s - K(theta_IC(z))] / [theta_s - theta_IC(z)].

    Parameters
    ----------
    soil : SoilParams
    H : float
        Column height [m]. Front starts at z = H.
    t_span : tuple (t0, t1)
    n_eval : int

    Returns
    -------
    t_arr, z_arr : ndarray
    """
    theta_tab, K_tab = build_K_table(soil)

    def rhs(t, z_vec):
        z = z_vec[0]
        if z <= 0:
            return [0.0]
        theta_ic = float(theta_from_h_np(np.array([-z]), soil)[0])
        K_ic = float(np.interp(theta_ic, theta_tab, K_tab))
        denom = soil.theta_s - theta_ic
        if denom < 1e-12:
            return [-1e6]
        s_local = (soil.K_s - K_ic) / denom
        return [-s_local]

    def front_exit(t, z_vec):
        return z_vec[0]
    front_exit.terminal = True
    front_exit.direction = -1

    t_eval = np.linspace(t_span[0], t_span[1], n_eval)
    sol = solve_ivp(rhs, t_span, [float(H)], t_eval=t_eval,
                    events=front_exit, method='RK45',
                    rtol=1e-8, atol=1e-10)

    t_arr = sol.t
    z_arr = np.clip(sol.y[0], 0.0, H)

    if t_arr[-1] < t_span[1]:
        t_arr = np.append(t_arr, t_span[1])
        z_arr = np.append(z_arr, 0.0)

    return t_arr, z_arr


def _build_profile_tools(mesh, H, W, maxh):
    """Build profile extraction helpers.

    Returns
    -------
    z_profile, pts, pts_top
    """
    n_profile = max(int(H / maxh * 2), 50)
    z_profile = np.linspace(maxh / 2, H - maxh / 2, n_profile)
    pts = [mesh(W / 2, zi) for zi in z_profile]
    eps_top = 1e-6
    x_samp = np.linspace(eps_top, W - eps_top, 10)
    pts_top = [mesh(xi, H - eps_top) for xi in x_samp]
    return z_profile, pts, pts_top


# ── Monolithic runner (clay baseline) ────────────────────────────────────


def run_monolithic(soil, H, maxh, dt, t_final, p_rain, snap_times,
                   W=0.1, verbose=True):
    """Run monolithic mixed FEM (no splitting) with adaptive dt.

    Uses build_monolithic_form_signorini: gravity is implicit in the flux
    equation, storage uses theta(h_old). Newton failure triggers dt halving.

    Parameters
    ----------
    soil : SoilParams
    H : float
        Column height [m].
    maxh : float
    dt : float
        Time step [s].
    t_final : float
    p_rain : float
        Precipitation rate [m/s].
    snap_times : array-like
        Times for profile snapshots [s].
    W : float
    verbose : bool

    Returns
    -------
    z_profile : ndarray
    snapshots : list of (float, ndarray)
        (time, theta_profile) pairs.
    info : dict
    """
    snap_times = np.sort(np.asarray(snap_times, dtype=float))

    mesh, Q_dg, fes = build_mesh_and_spaces(H, soil, W=W, maxh=maxh)
    V, Q_mix = fes.components

    dt_param = Parameter(dt)
    gfu = GridFunction(fes)
    gf_q, gf_h = gfu.components
    gf_h_old = GridFunction(Q_mix)

    R = build_monolithic_form_signorini(
        fes, mesh, soil, dt_param, gf_h_old, p_rain=p_rain)

    gfu_backup = gfu.vec.CreateVector()

    # IC: hydrostatic h = -z
    gf_h.Set(-y)
    gf_h_old.Set(-y)

    z_profile, pts, pts_top = _build_profile_tools(mesh, H, W, maxh)

    def extract_theta():
        h_arr = np.array([gf_h(p) for p in pts])
        return theta_from_h_np(h_arr, soil)

    def extract_h():
        return np.array([gf_h(p) for p in pts])

    def extract_psi():
        return np.array([gf_h(p) for p in pts])

    def extract_qz():
        return np.array([gf_q(p)[1] for p in pts])

    t = 0.0
    step = 0
    snap_idx = 0
    dt_current = dt
    dt_macro = dt

    newton_its_log = []
    newton_failures = 0
    failed = False
    failure_time = None

    snapshots = [(0.0, extract_theta())]
    psi_snapshots = [(0.0, extract_psi())]
    qz_snapshots = [(0.0, extract_qz())]

    while t < t_final - 1e-12:
        dt_eff = min(dt_current, t_final - t)
        dt_param.Set(dt_eff)

        gfu_backup.data = gfu.vec
        gf_h_old.vec.FV().NumPy()[:] = gf_h.vec.FV().NumPy()

        with TaskManager():
            status, nits = Newton(R, gfu,
                                  freedofs=fes.FreeDofs(),
                                  maxerr=1e-6, maxit=50,
                                  inverse="pardiso",
                                  printing=False)

        if status != 0:
            gfu.vec.data = gfu_backup
            newton_failures += 1
            dt_current = dt_eff / 2.0
            if dt_current < 1e-4:
                failed = True
                failure_time = t
                break
            continue

        t += dt_eff
        step += 1
        newton_its_log.append(nits)

        if nits <= 5 and dt_current < dt_macro:
            dt_current = min(dt_current * 1.5, dt_macro)

        while (snap_idx < len(snap_times)
               and t >= snap_times[snap_idx] - 1e-10):
            snapshots.append((t, extract_theta()))
            psi_snapshots.append((t, extract_psi()))
            qz_snapshots.append((t, extract_qz()))
            snap_idx += 1

        if verbose and step % 20 == 0:
            print(f"  step {step:4d}, t={t:.1f}s "
                  f"({t / 3600:.2f}h), Newton={nits}")

    info = {
        "n_steps": step,
        "n_els": mesh.ne,
        "newton_its": np.array(newton_its_log),
        "newton_failures": newton_failures,
        "failed": failed,
        "failure_time": failure_time,
        "psi_snapshots": psi_snapshots,
        "qz_snapshots": qz_snapshots,
    }
    return z_profile, snapshots, info


# ── Splitting runner (constant rain) ────────────────────────────────────


def run_splitting(soil, H, maxh, dt_macro, t_final, p_rain, snap_times,
                  W=0.1, cfl=0.5, maxerr=1e-3, verbose=True,
                  uniform_dry_ic=False):
    """Run Lie-Trotter splitting with Nitsche-Signorini.

    Godunov convection (explicit subcycled) + mixed FEM diffusion
    (implicit Newton). Constant p_rain (seepage mode, ghost = theta_s).

    Parameters
    ----------
    soil : SoilParams
    H : float
    maxh : float
    dt_macro : float
    t_final : float
    p_rain : float
    snap_times : array-like
    W : float
    cfl : float
    maxerr : float
    verbose : bool
    uniform_dry_ic : bool
        If True, use h = -H everywhere. Default: hydrostatic h = -z.

    Returns
    -------
    z_profile : ndarray
    snapshots : list of (float, ndarray)
    info : dict
    """
    snap_times = np.sort(np.asarray(snap_times, dtype=float))

    mesh, Q_dg, fes = build_mesh_and_spaces(H, soil, W=W, maxh=maxh)
    V, Q_mix = fes.components

    a_flux = build_godunov_flux(mesh, Q_dg, soil)
    dt_cfl, _ = compute_cfl_dt(soil, maxh, cfl)
    gf_theta = GridFunction(Q_dg)

    dt_param = Parameter(dt_macro)
    gf_theta_star = GridFunction(Q_mix)
    gfu = GridFunction(fes)
    gf_q, gf_h = gfu.components
    R_diff = build_diffusion_form_signorini(
        fes, mesh, soil, dt_param, gf_theta_star, p_rain=p_rain)

    gf_theta_backup = gf_theta.vec.CreateVector()
    gfu_backup = gfu.vec.CreateVector()

    if uniform_dry_ic:
        gf_h.Set(CF(-H))
        gf_theta.Set(cons.theta(CF(-H), soil))
    else:
        gf_h.Set(-y)
        gf_theta.Set(cons.theta(-y, soil))

    z_profile, pts, pts_top = _build_profile_tools(mesh, H, W, maxh)

    def extract_theta():
        return np.array([gf_theta(p) for p in pts])

    def extract_h():
        return np.array([gf_h(p) for p in pts])

    def eval_h_top():
        return np.mean([gf_h(p) for p in pts_top])

    t = 0.0
    step = 0
    snap_idx = 0
    dt_current = dt_macro

    snapshots = [(0.0, extract_theta())]
    h_snapshots = [(0.0, extract_h())]
    newton_its_log = []
    front_log = [(0.0, float(H))]

    while t < t_final - 1e-12:
        dt_eff = min(dt_current, t_final - t)
        dt_param.Set(dt_eff)

        gf_theta_backup.data = gf_theta.vec
        gfu_backup.data = gfu.vec

        godunov_subcycle(gf_theta, a_flux, Q_dg, soil, dt_eff, dt_cfl)

        gf_theta_star.vec.FV().NumPy()[:] = gf_theta.vec.FV().NumPy()
        h_star = h_from_theta_np(gf_theta.vec.FV().NumPy(), soil)
        gf_h.vec.FV().NumPy()[:] = h_star

        with TaskManager():
            status, nits = Newton(R_diff, gfu,
                                  freedofs=fes.FreeDofs(),
                                  maxerr=maxerr, maxit=25,
                                  inverse="pardiso",
                                  printing=False)

        if status != 0:
            gf_theta.vec.data = gf_theta_backup
            gfu.vec.data = gfu_backup
            dt_current = dt_eff / 2.0
            if dt_current < 1e-4:
                print(f"ABORT: dt too small at t={t:.1f}s")
                break
            continue

        theta_new = theta_from_h_np(gf_h.vec.FV().NumPy(), soil)
        gf_theta.vec.FV().NumPy()[:] = theta_new

        t += dt_eff
        step += 1
        newton_its_log.append(nits)

        if nits <= 5 and dt_current < dt_macro:
            dt_current = min(dt_current * 1.5, dt_macro)

        while (snap_idx < len(snap_times)
               and t >= snap_times[snap_idx] - 1e-10):
            snapshots.append((t, extract_theta()))
            h_snapshots.append((t, extract_h()))
            th_prof = extract_theta()
            front_pos = find_front_position(th_prof, z_profile, soil)
            front_log.append((t, front_pos))
            snap_idx += 1

        if verbose and step % 30 == 0:
            h_top = eval_h_top()
            print(f"  step {step:4d}, t={t:.1f}s "
                  f"({t / 3600:.2f}h), Newton={nits}, "
                  f"h_top={h_top:.4f}")

    h_max_all = max(np.max(hs) for _, hs in h_snapshots)

    info = {
        "dt_cfl": dt_cfl,
        "dt_macro": dt_macro,
        "n_steps": step,
        "n_els": mesh.ne,
        "newton_its": np.array(newton_its_log),
        "h_snapshots": h_snapshots,
        "front_log": np.array(front_log),
        "h_max": h_max_all,
    }
    return z_profile, snapshots, info


# ── Splitting runner with time-varying rain ──────────────────────────────


def run_splitting_ramp(soil, H, maxh, dt_macro, t_final, p_rain_func,
                       snap_times, W=0.1, cfl=0.5, verbose=True):
    """Run Lie-Trotter splitting with dynamic ghost and time-varying rain.

    Key features: ghost cap (K(ghost) <= p_rain), IC-consistent ghost
    initialization, regularized Signorini, relaxed Newton tolerance.

    Parameters
    ----------
    soil : SoilParams
    H : float
    maxh : float
    dt_macro : float
    t_final : float
    p_rain_func : callable
        p_rain_func(t) -> float.
    snap_times : array-like
    W : float
    cfl : float
    verbose : bool

    Returns
    -------
    z_profile : ndarray
    snapshots : list of (float, ndarray)
    info : dict
    """
    snap_times = np.sort(np.asarray(snap_times, dtype=float))

    mesh, Q_dg, fes = build_mesh_and_spaces(H, soil, W=W, maxh=maxh)
    V, Q_mix = fes.components

    theta_ghost = Parameter(soil.theta_s)
    a_flux = build_godunov_flux(mesh, Q_dg, soil, theta_top=theta_ghost)
    dt_cfl, max_speed = compute_cfl_dt(soil, maxh, cfl)
    gf_theta = GridFunction(Q_dg)

    dt_param = Parameter(dt_macro)
    p_rain_param = Parameter(p_rain_func(0.0))
    gf_theta_star = GridFunction(Q_mix)
    gfu = GridFunction(fes)
    gf_q, gf_h = gfu.components
    R_diff = build_diffusion_form_signorini(
        fes, mesh, soil, dt_param, gf_theta_star, p_rain=p_rain_param)

    gf_theta_backup = gf_theta.vec.CreateVector()
    gfu_backup = gfu.vec.CreateVector()

    # IC: hydrostatic h = -z
    gf_h.Set(-y)
    gf_theta.Set(cons.theta(-y, soil))

    # Ghost initialization consistent with p_rain(0)
    theta_tab, K_tab = build_K_table(soil)
    p_init = p_rain_func(0.0)
    if p_init < soil.K_s:
        if p_init > 0:
            theta_ghost.Set(theta_from_K_table(p_init, soil, theta_tab, K_tab))
        else:
            h_top_ic = -H
            theta_ghost.Set(float(
                theta_from_h_np(np.array([h_top_ic]), soil)[0]))

    z_profile, pts, pts_top = _build_profile_tools(mesh, H, W, maxh)

    def extract_theta():
        return np.array([gf_theta(p) for p in pts])

    def extract_h():
        return np.array([gf_h(p) for p in pts])

    def eval_h_top():
        return np.mean([gf_h(p) for p in pts_top])

    t = 0.0
    step = 0
    snap_idx = 0
    dt_current = dt_macro

    snapshots = [(0.0, extract_theta())]
    h_snapshots = [(0.0, extract_h())]
    newton_its_log = []
    h_top_history = []
    t_history = [0.0]
    p_rain_history = [p_rain_func(0.0)]
    signorini_mode_history = []

    while t < t_final - 1e-12:
        dt_eff = min(dt_current, t_final - t)
        dt_param.Set(dt_eff)

        p_val = p_rain_func(t)
        p_rain_param.Set(p_val)

        # Ghost cap: K(ghost) <= p when 0 < p < K_s
        if 0 < p_val < soil.K_s:
            theta_cap = theta_from_K_table(p_val, soil, theta_tab, K_tab)
            if float(theta_ghost.Get()) > theta_cap:
                theta_ghost.Set(theta_cap)

        gf_theta_backup.data = gf_theta.vec
        gfu_backup.data = gfu.vec

        godunov_subcycle(gf_theta, a_flux, Q_dg, soil, dt_eff, dt_cfl)

        gf_theta_star.vec.FV().NumPy()[:] = gf_theta.vec.FV().NumPy()
        h_star = h_from_theta_np(gf_theta.vec.FV().NumPy(), soil)
        gf_h.vec.FV().NumPy()[:] = h_star

        with TaskManager():
            status, nits = Newton(R_diff, gfu,
                                  freedofs=fes.FreeDofs(),
                                  maxerr=5e-3, maxit=50,
                                  inverse="pardiso",
                                  printing=False)

        if status != 0:
            gf_theta.vec.data = gf_theta_backup
            gfu.vec.data = gfu_backup
            dt_current = dt_eff / 2.0
            if verbose:
                print(f"  Newton failed at t={t:.1f}s, "
                      f"halving dt to {dt_current:.2f}s")
            if dt_current < 1e-4:
                print(f"ABORT: dt too small at t={t:.1f}s")
                break
            continue

        theta_new = theta_from_h_np(gf_h.vec.FV().NumPy(), soil)
        gf_theta.vec.FV().NumPy()[:] = theta_new

        h_top = eval_h_top()
        theta_top_val = float(
            theta_from_h_np(np.array([h_top]), soil)[0])
        theta_ghost.Set(theta_top_val)

        t += dt_eff
        step += 1
        newton_its_log.append(nits)
        h_top_history.append(h_top)
        t_history.append(t)
        p_rain_history.append(p_val)
        mode = "seepage" if abs(h_top) < 1e-2 else "infiltration"
        signorini_mode_history.append(mode)

        if nits <= 5 and dt_current < dt_macro:
            dt_current = min(dt_current * 1.5, dt_macro)

        while (snap_idx < len(snap_times)
               and t >= snap_times[snap_idx] - 1e-10):
            snapshots.append((t, extract_theta()))
            h_snapshots.append((t, extract_h()))
            snap_idx += 1

        if verbose and step % 30 == 0:
            print(f"  step {step:4d}, t={t:.1f}s "
                  f"({t / 3600:.2f}h), Newton={nits}, "
                  f"h_top={h_top:.4f}, p={p_val:.2e}, "
                  f"mode={mode}")

    h_max_all = max(np.max(hs) for _, hs in h_snapshots)

    info = {
        "n_steps": step,
        "n_els": mesh.ne,
        "newton_its": np.array(newton_its_log),
        "h_snapshots": h_snapshots,
        "h_max": h_max_all,
        "h_top_history": np.array(h_top_history),
        "t_history": np.array(t_history),
        "p_rain_history": np.array(p_rain_history),
        "signorini_mode": signorini_mode_history,
    }
    return z_profile, snapshots, info


# ── Pure Godunov runner (no diffusion) ──────────────────────────────────


def run_godunov_only(soil, H, maxh, dt_macro, t_final, snap_times,
                     W=0.1, cfl=0.5):
    """Pure Godunov convection (no diffusion step).

    Uses uniform dry IC (h = -H) and theta_s ghost at top (seepage).

    Parameters
    ----------
    soil : SoilParams
    H : float
    maxh : float
    dt_macro : float
    t_final : float
    snap_times : array-like
    W : float
    cfl : float

    Returns
    -------
    z_profile : ndarray
    snapshots : list of (float, ndarray)
    front_log : ndarray
    """
    snap_times = np.sort(np.asarray(snap_times, dtype=float))

    mesh, Q_dg, fes = build_mesh_and_spaces(H, soil, W=W, maxh=maxh)
    a_flux = build_godunov_flux(mesh, Q_dg, soil)
    dt_cfl, _ = compute_cfl_dt(soil, maxh, cfl)
    gf_theta = GridFunction(Q_dg)

    gf_theta.Set(cons.theta(CF(-H), soil))

    z_profile, pts, _ = _build_profile_tools(mesh, H, W, maxh)

    def extract_theta():
        return np.array([gf_theta(p) for p in pts])

    t = 0.0
    snap_idx = 0

    snapshots = [(0.0, extract_theta())]
    front_log = [(0.0, float(H))]

    while t < t_final - 1e-12:
        dt_eff = min(dt_macro, t_final - t)
        godunov_subcycle(gf_theta, a_flux, Q_dg, soil, dt_eff, dt_cfl)
        th_vec = gf_theta.vec.FV().NumPy()
        np.clip(th_vec, soil.theta_r, soil.theta_s, out=th_vec)

        t += dt_eff

        while (snap_idx < len(snap_times)
               and t >= snap_times[snap_idx] - 1e-10):
            th_prof = extract_theta()
            snapshots.append((t, th_prof))
            front_pos = find_front_position(th_prof, z_profile, soil)
            front_log.append((t, front_pos))
            snap_idx += 1

    return z_profile, snapshots, np.array(front_log)


# ── Manufactured-solution convergence ────────────────────────────────────


def _build_column_geometry(W, H):
    """Build SplineGeometry for a W x H column."""
    geo = SplineGeometry()
    pts = [geo.AppendPoint(px, py) for px, py in
           [(0, 0), (W, 0), (W, H), (0, H)]]
    geo.Append(["line", pts[0], pts[1]], bc="bottom")
    geo.Append(["line", pts[1], pts[2]], bc="right")
    geo.Append(["line", pts[2], pts[3]], bc="top")
    geo.Append(["line", pts[3], pts[0]], bc="left")
    return geo


def solve_manufactured(mesh_obj, fes, soil, H=2.0, h_center=0.5,
                       h_amp=0.25):
    """Stationary manufactured-solution solve on pre-built mesh.

    Manufactured solution: h_ex(y) = -(h_center + h_amp * sin(pi*y/H)).

    Parameters
    ----------
    mesh_obj : Mesh
    fes : ProductSpace (HDiv x L2)
    soil : SoilParams
    H, h_center, h_amp : float

    Returns
    -------
    dict with err_h, err_q, err_theta_L1, err_theta_L2, nits, ne
    """
    (q, h), (tau, v) = fes.TnT()
    nn = specialcf.normal(mesh_obj.dim)
    ez = ez_cf(mesh_obj)

    h_ex = -(h_center + h_amp * sin(np.pi * y / H))
    dh_dy = -h_amp * np.pi / H * cos(np.pi * y / H)
    K_ex = cons.K(h_ex, soil)
    q_ex = CF((0, -K_ex * (dh_dy + 1.0)))
    theta_ex = cons.theta(h_ex, soil)

    R = BilinearForm(fes)
    R += cons.Kinv(h, soil) * InnerProduct(q, tau) * dx
    R += -h * div(tau) * dx
    R += InnerProduct(ez, tau) * dx
    R += h_ex * (tau * nn) * ds(skeleton=True,
                                definedon=mesh_obj.Boundaries("bottom|top"))
    R += div(q) * v * dx
    R += -InnerProduct(q_ex, nn) * v * dx(element_boundary=True)

    gfu = GridFunction(fes)
    gf_q, gf_h = gfu.components
    gf_h.Set(h_ex)

    with TaskManager():
        status, nits = Newton(R, gfu, freedofs=fes.FreeDofs(),
                              maxerr=1e-8, maxit=50,
                              inverse="pardiso", printing=False)
    assert status == 0, f"Newton failed: status={status}, nits={nits}"

    err_h = float(sqrt(Integrate((gf_h - h_ex) ** 2, mesh_obj)))
    err_q = float(sqrt(Integrate(InnerProduct(gf_q - q_ex, gf_q - q_ex),
                                 mesh_obj)))
    theta_h = cons.theta(gf_h, soil)
    diff_theta = theta_h - theta_ex
    err_theta_L2 = float(sqrt(Integrate(diff_theta ** 2, mesh_obj)))
    err_theta_L1 = float(Integrate(IfPos(diff_theta, diff_theta,
                                         -diff_theta), mesh_obj))

    return {
        "err_h": err_h, "err_q": err_q,
        "err_theta_L2": err_theta_L2, "err_theta_L1": err_theta_L1,
        "nits": nits, "ne": mesh_obj.ne,
    }


def run_spatial_convergence(surrogates, n_levels=5, H=2.0, W=0.1,
                            maxh_coarse=0.5, verbose=True):
    """Run spatial convergence study across surrogates.

    Parameters
    ----------
    surrogates : dict
        {label: SoilParams} pairs.
    n_levels : int
    H, W : float
    maxh_coarse : float
        Coarsest mesh size.
    verbose : bool

    Returns
    -------
    results : dict
        {label: {norm: array, ...}} for each surrogate.
    """
    norms = ["err_h", "err_q", "err_theta_L1", "err_theta_L2"]
    results = {label: {k: [] for k in norms + ["ne", "nits"]}
               for label in surrogates}

    geo = _build_column_geometry(W, H)
    mesh = Mesh(geo.GenerateMesh(maxh=maxh_coarse))

    for level in range(n_levels):
        if level > 0:
            mesh.Refine()
        V = HDiv(mesh, order=0, RT=True, dirichlet="left|right")
        Q = L2(mesh, order=0)
        fes = V * Q

        for label, soil in surrogates.items():
            r = solve_manufactured(mesh, fes, soil, H=H)
            for k in norms:
                results[label][k].append(r[k])
            results[label]["ne"].append(r["ne"])
            results[label]["nits"].append(r["nits"])

        if verbose:
            print(f"Level {level}: {mesh.ne} elements -- all surrogates OK")

    for label in surrogates:
        for k in results[label]:
            results[label][k] = np.array(results[label][k])

    return results


# ── Temporal convergence ─────────────────────────────────────────────────


def run_temporal_convergence(soil, H, maxh, dt_levels, dt_ref, t_final,
                             p_rain, W=0.1, verbose=True):
    """Run temporal self-convergence study.

    Monolithic backward Euler on fixed mesh, comparing each dt against
    a fine-dt reference.

    Parameters
    ----------
    soil : SoilParams
    H, maxh : float
    dt_levels : list of float
    dt_ref : float
        Reference (fine) dt.
    t_final : float
    p_rain : float
    W : float
    verbose : bool

    Returns
    -------
    results : list of dict
        Each dict: {dt, err_h, err_theta_L1}.
    """
    geo = _build_column_geometry(W, H)
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))
    V = HDiv(mesh, order=0, RT=True, dirichlet="left|right")
    Q = L2(mesh, order=0)
    fes = V * Q

    def _run_fixed_dt(dt_val):
        gfu = GridFunction(fes)
        gf_q, gf_h = gfu.components
        gf_h_old = GridFunction(Q)
        gf_h.Set(-y)
        gf_h_old.Set(-y)

        dt_param = Parameter(dt_val)
        R = build_monolithic_form_signorini(
            fes, mesh, soil, dt_param, gf_h_old, p_rain=p_rain)

        n_steps = int(round(t_final / dt_val))
        for step in range(n_steps):
            gf_h_old.vec.FV().NumPy()[:] = gf_h.vec.FV().NumPy()
            with TaskManager():
                status, nits = Newton(R, gfu, freedofs=fes.FreeDofs(),
                                      maxerr=1e-6, maxit=50,
                                      inverse="pardiso", printing=False)
            assert status == 0, f"Newton failed at step {step}"
        return gf_h.vec.FV().NumPy().copy()

    if verbose:
        print(f"Mesh: {mesh.ne} elements")
        print(f"Reference: dt={dt_ref}s ...", end=" ", flush=True)
    h_ref = _run_fixed_dt(dt_ref)
    theta_ref = theta_from_h_np(h_ref, soil)
    if verbose:
        print("done")

    gf_err = GridFunction(Q)
    results = []
    for dt_test in dt_levels:
        if verbose:
            n_steps = int(t_final / dt_test)
            print(f"dt={dt_test:.0f}s ({n_steps} steps) ...",
                  end=" ", flush=True)
        h_test = _run_fixed_dt(dt_test)
        theta_test = theta_from_h_np(h_test, soil)

        gf_err.vec.FV().NumPy()[:] = h_test - h_ref
        err_h = float(sqrt(Integrate(gf_err * gf_err, mesh)))

        gf_err.vec.FV().NumPy()[:] = np.abs(theta_test - theta_ref)
        err_theta_L1 = float(Integrate(gf_err, mesh))

        results.append({"dt": dt_test, "err_h": err_h,
                         "err_theta_L1": err_theta_L1})
        if verbose:
            print(f"err_h={err_h:.4e}, err_th={err_theta_L1:.4e}")

    return results


# ── Splitting temporal convergence ──────────────────────────────────────


def run_splitting_temporal_convergence(soil, H, maxh, dt_levels, dt_ref,
                                       t_final, p_rain, W=0.1, cfl=0.5,
                                       maxerr=1e-3, verbose=True):
    """Temporal self-convergence for Lie-Trotter splitting.

    Runs splitting at multiple dt_macro levels on a shared mesh, comparing
    each against a fine-dt_macro reference. Isolates the splitting temporal
    error: spatial discretization is identical across all runs.

    The Godunov substep is CFL-limited (independent of dt_macro), so only
    the diffusion-correction frequency varies.

    Parameters
    ----------
    soil : SoilParams
    H, maxh : float
    dt_levels : list of float
        Macro time steps to test [s].
    dt_ref : float
        Reference (fine) dt_macro [s].
    t_final : float
        Final time [s], should be within transient phase.
    p_rain : float
        Precipitation rate [m/s].
    W : float
    cfl : float
    maxerr : float
        Newton tolerance for diffusion step.
    verbose : bool

    Returns
    -------
    results : list of dict
        Each dict: {dt, err_theta_L1, err_h_L2, steps}.
    """
    mesh, Q_dg, fes = build_mesh_and_spaces(H, soil, W=W, maxh=maxh)
    V, Q_mix = fes.components

    a_flux = build_godunov_flux(mesh, Q_dg, soil)
    dt_cfl, _ = compute_cfl_dt(soil, maxh, cfl)

    def _run_fixed_dt(dt_macro_val):
        """Run splitting to t_final, return (theta_dofs, h_dofs, n_steps)."""
        dt_param = Parameter(dt_macro_val)
        gf_theta = GridFunction(Q_dg)
        gf_theta_star = GridFunction(Q_mix)
        gfu = GridFunction(fes)
        gf_q, gf_h = gfu.components

        gf_theta_backup = gf_theta.vec.CreateVector()
        gfu_backup = gfu.vec.CreateVector()

        # IC: hydrostatic h = -z
        gf_h.Set(-y)
        gf_theta.Set(cons.theta(-y, soil))

        R_diff = build_diffusion_form_signorini(
            fes, mesh, soil, dt_param, gf_theta_star, p_rain=p_rain)

        t = 0.0
        step = 0
        dt_current = dt_macro_val

        while t < t_final - 1e-12:
            dt_eff = min(dt_current, t_final - t)
            dt_param.Set(dt_eff)

            gf_theta_backup.data = gf_theta.vec
            gfu_backup.data = gfu.vec

            godunov_subcycle(gf_theta, a_flux, Q_dg, soil, dt_eff, dt_cfl)

            gf_theta_star.vec.FV().NumPy()[:] = (
                gf_theta.vec.FV().NumPy())
            h_star = h_from_theta_np(
                gf_theta.vec.FV().NumPy(), soil)
            gf_h.vec.FV().NumPy()[:] = h_star

            with TaskManager():
                status, nits = Newton(
                    R_diff, gfu,
                    freedofs=fes.FreeDofs(),
                    maxerr=maxerr, maxit=25,
                    inverse="pardiso", printing=False)

            if status != 0:
                gf_theta.vec.data = gf_theta_backup
                gfu.vec.data = gfu_backup
                dt_current = dt_eff / 2.0
                if dt_current < 1e-4:
                    raise RuntimeError(
                        f"dt too small at t={t:.1f}s")
                continue

            theta_new = theta_from_h_np(
                gf_h.vec.FV().NumPy(), soil)
            gf_theta.vec.FV().NumPy()[:] = theta_new

            t += dt_eff
            step += 1

            if nits <= 5 and dt_current < dt_macro_val:
                dt_current = min(dt_current * 1.5, dt_macro_val)

        return (gf_theta.vec.FV().NumPy().copy(),
                gf_h.vec.FV().NumPy().copy(), step)

    if verbose:
        Gr = H * soil.alpha
        print(f"H={H}m, Gr={Gr:.1f}, mesh: {mesh.ne} els, "
              f"dt_cfl={dt_cfl:.3f}s")
        print(f"Reference: dt_macro={dt_ref}s ...",
              end=" ", flush=True)

    theta_ref, h_ref, steps_ref = _run_fixed_dt(dt_ref)
    if verbose:
        print(f"done ({steps_ref} steps)")

    gf_err = GridFunction(Q_mix)
    results = []

    for dt_test in dt_levels:
        if verbose:
            n_exp = int(round(t_final / dt_test))
            print(f"dt_macro={dt_test:.0f}s (~{n_exp} steps) ...",
                  end=" ", flush=True)

        theta_test, h_test, steps_test = _run_fixed_dt(dt_test)

        # L1(theta)
        gf_err.vec.FV().NumPy()[:] = np.abs(theta_test - theta_ref)
        err_theta_L1 = float(Integrate(gf_err, mesh))

        # L2(h)
        gf_err.vec.FV().NumPy()[:] = h_test - h_ref
        err_h_L2 = float(sqrt(Integrate(gf_err * gf_err, mesh)))

        results.append({
            "dt": dt_test,
            "err_theta_L1": err_theta_L1,
            "err_h_L2": err_h_L2,
            "steps": steps_test,
        })
        if verbose:
            print(f"err_th={err_theta_L1:.4e}, "
                  f"err_h={err_h_L2:.4e} ({steps_test} steps)")

    return results
