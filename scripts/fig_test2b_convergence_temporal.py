#!/usr/bin/env python
"""
Test 2b: Temporal convergence (monolithic backward Euler).

Self-convergence study on the D2 surrogate (beta=2.75, near-transition).
Fixed fine mesh, varying dt, reference at dt=15s.

Output: figures/fig_test2b_convergence_temporal.pdf + rates table to stdout
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from _style import DATA_DIR, FIGURES_DIR, apply_style

from nitsche_signorini import SoilParams

apply_style()

# ── Parameters ─────────────────────────────────────────────────────────
soil = SoilParams(theta_s=0.4, theta_r=0.05, alpha=1.0, n=1.8,
                  ell=0.5, K_s=1e-4)  # D2 surrogate
H = 2.0
maxh = 0.05
p_rain = 0.5 * soil.K_s
t_final = 3600.0
dt_levels = [240.0, 120.0, 60.0, 30.0]
dt_ref = 15.0

CACHE = DATA_DIR / "test2b_convergence_temporal.npz"
recompute = "--recompute" in sys.argv

print("Test 2b: Temporal convergence (D2 surrogate)")
print(f"  H={H}m, maxh={maxh}m, p/K_s=0.5, t_final={t_final/3600:.0f}h")
print(f"  dt levels: {dt_levels}, reference: {dt_ref}s")
print()

if not recompute and CACHE.exists():
    print(f"Loading cached data from {CACHE}")
    d = np.load(CACHE)
    dt_arr = d["dt_arr"]
    err_h_arr = d["err_h_arr"]
    err_theta_arr = d["err_theta_arr"]
else:
    from nitsche_signorini.runners import run_temporal_convergence

    results = run_temporal_convergence(
        soil, H, maxh, dt_levels, dt_ref, t_final, p_rain)

    dt_arr = np.array([r["dt"] for r in results])
    err_h_arr = np.array([r["err_h"] for r in results])
    err_theta_arr = np.array([r["err_theta_L1"] for r in results])

    CACHE.parent.mkdir(exist_ok=True)
    np.savez(CACHE, dt_arr=dt_arr, err_h_arr=err_h_arr,
             err_theta_arr=err_theta_arr)
    print(f"Saved: {CACHE}")

# ── Print convergence table ────────────────────────────────────────────
print(f"\n{'dt':>8} {'||h-href||':>12} {'rate':>6} "
      f"{'||th-thref||':>14} {'rate':>6}")
print("-" * 52)
for i in range(len(dt_arr)):
    rh = (np.log2(err_h_arr[i - 1] / err_h_arr[i])
          if i > 0 else float("nan"))
    rt = (np.log2(err_theta_arr[i - 1] / err_theta_arr[i])
          if i > 0 else float("nan"))
    print(f"{dt_arr[i]:>8.0f} {err_h_arr[i]:>12.4e} {rh:>6.2f} "
          f"{err_theta_arr[i]:>14.4e} {rt:>6.2f}")

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

ax = axes[0]
ax.loglog(dt_arr, err_h_arr, "o-", lw=1.5, ms=7)
dt_line = np.array([dt_arr[0], dt_arr[-1]])
ref = err_h_arr[0] * (dt_line / dt_line[0])
ax.loglog(dt_line, ref, "k:", lw=1, label="O(dt)")
ax.set_xlabel(r"$\Delta t$ [s]")
ax.set_ylabel(r"$\|h - h_{\mathrm{ref}}\|_{L^2}$")
ax.set_title(r"(a) $L^2(h)$ self-convergence")
ax.legend()

ax = axes[1]
ax.loglog(dt_arr, err_theta_arr, "s-", lw=1.5, ms=7)
ref = err_theta_arr[0] * (dt_line / dt_line[0])
ax.loglog(dt_line, ref, "k:", lw=1, label="O(dt)")
ax.set_xlabel(r"$\Delta t$ [s]")
ax.set_ylabel(r"$\|\theta - \theta_{\mathrm{ref}}\|_{L^1}$")
ax.set_title(r"(b) $L^1(\theta)$ self-convergence")
ax.legend()

fig.suptitle("Test 2b: Temporal convergence (D2, monolithic backward Euler)",
             fontsize=11)
fig.tight_layout()

out = Path(FIGURES_DIR) / "fig_test2b_convergence_temporal.pdf"
out.parent.mkdir(exist_ok=True)
fig.savefig(out)
print(f"\nSaved: {out}")
plt.close()
