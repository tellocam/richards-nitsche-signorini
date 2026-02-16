#!/usr/bin/env python
"""
NOT USED IN PAPER -- Diagnostic confirming that splitting temporal
convergence is unmeasurable. See scripts/README.md for context.

Cauchy convergence (no reference) between consecutive dt levels. Result:
differences are flat (~1.5e-4 for Gr=2), giving rates near zero. Confirms
the splitting temporal error is too small to measure for this soil.
Kept for reference only.

---

Original intent: Compute ||theta_{dt} - theta_{dt/2}||_L1 between
consecutive dt levels and extract convergence rates directly, avoiding
reference contamination issues found in fig_test2_splitting_temporal.py.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from ngsolve import (
    CF, GridFunction, Integrate, Parameter, TaskManager, y,
)
from ngsolve.solvers import Newton

from nitsche_signorini import SoilParams
from nitsche_signorini.constitutive import (
    ConstitutiveSharp, h_from_theta_np, theta_from_h_np,
)
from nitsche_signorini.splitting import (
    build_diffusion_form_signorini,
    build_godunov_flux,
    build_mesh_and_spaces,
    compute_cfl_dt,
    godunov_subcycle,
)

cons = ConstitutiveSharp()

soil = SoilParams(theta_s=0.40, theta_r=0.04, alpha=1.0,
                  n=1.5, ell=0.5, K_s=1e-6)

configs = [
    {"label": "Gr=2",  "H": 2.0,  "Gr": 2.0},
    {"label": "Gr=4",  "H": 4.0,  "Gr": 4.0},
    {"label": "Gr=10", "H": 10.0, "Gr": 10.0},
]

W = 0.1
maxh = 0.05
t_final = 36000.0
dt_levels = [7200.0, 3600.0, 1800.0, 900.0, 450.0]
p_rain = 10.0 * soil.K_s
maxerr = 1e-6

print("Cauchy convergence diagnostic")
print(f"  dt levels: {[f'{d:.0f}' for d in dt_levels]}s")
print()


def run_fixed_dt(mesh, Q_dg, fes, dt_macro_val):
    """Run splitting to t_final, return theta dof vector."""
    V, Q_mix = fes.components
    dt_param = Parameter(dt_macro_val)
    gf_theta = GridFunction(Q_dg)
    gf_theta_star = GridFunction(Q_mix)
    gfu = GridFunction(fes)
    gf_q, gf_h = gfu.components

    gf_theta_backup = gf_theta.vec.CreateVector()
    gfu_backup = gfu.vec.CreateVector()

    gf_h.Set(-y)
    gf_theta.Set(cons.theta(-y, soil))

    a_flux = build_godunov_flux(mesh, Q_dg, soil)
    dt_cfl, _ = compute_cfl_dt(soil, maxh, 0.5)

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

        gf_theta_star.vec.FV().NumPy()[:] = gf_theta.vec.FV().NumPy()
        h_star = h_from_theta_np(gf_theta.vec.FV().NumPy(), soil)
        gf_h.vec.FV().NumPy()[:] = h_star

        with TaskManager():
            status, nits = Newton(
                R_diff, gfu, freedofs=fes.FreeDofs(),
                maxerr=maxerr, maxit=25,
                inverse="pardiso", printing=False)

        if status != 0:
            gf_theta.vec.data = gf_theta_backup
            gfu.vec.data = gfu_backup
            dt_current = dt_eff / 2.0
            if dt_current < 1e-4:
                raise RuntimeError(f"dt too small at t={t:.1f}s")
            continue

        theta_new = theta_from_h_np(gf_h.vec.FV().NumPy(), soil)
        gf_theta.vec.FV().NumPy()[:] = theta_new

        t += dt_eff
        step += 1

        if nits <= 5 and dt_current < dt_macro_val:
            dt_current = min(dt_current * 1.5, dt_macro_val)

    return gf_theta.vec.FV().NumPy().copy(), step


for c in configs:
    H = c["H"]
    area = H * W
    print(f"--- {c['label']} (H={H:.0f}m) ---")

    mesh, Q_dg, fes = build_mesh_and_spaces(H, soil, W=W, maxh=maxh)
    V, Q_mix = fes.components
    gf_err = GridFunction(Q_mix)

    theta_vecs = []
    for dt_val in dt_levels:
        n_exp = int(round(t_final / dt_val))
        print(f"  dt={dt_val:.0f}s (~{n_exp} steps) ...", end=" ", flush=True)
        theta, steps = run_fixed_dt(mesh, Q_dg, fes, dt_val)
        theta_vecs.append(theta)
        print(f"done ({steps} steps)")

    # Cauchy differences
    diffs = []
    for i in range(len(dt_levels) - 1):
        gf_err.vec.FV().NumPy()[:] = np.abs(theta_vecs[i] - theta_vecs[i + 1])
        diff_L1 = float(Integrate(gf_err, mesh))
        diffs.append(diff_L1)

    # Cauchy rates
    print(f"\n  {'dt1':>6} {'dt2':>6} {'||th1-th2||_L1':>16} "
          f"{'normalized':>12} {'Cauchy rate':>12}")
    print("  " + "-" * 58)
    for i, diff in enumerate(diffs):
        dt1, dt2 = dt_levels[i], dt_levels[i + 1]
        rate = (np.log2(diffs[i - 1] / diffs[i])
                if i > 0 and diffs[i] > 0 else float("nan"))
        print(f"  {dt1:6.0f} {dt2:6.0f} {diff:16.4e} "
              f"{diff/area:12.4e} {rate:12.2f}")
    print()
