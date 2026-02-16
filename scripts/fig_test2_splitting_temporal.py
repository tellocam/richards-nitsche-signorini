#!/usr/bin/env python
"""
NOT USED IN PAPER -- Splitting temporal convergence cannot be demonstrated
cleanly. See scripts/README.md for full explanation.

The splitting temporal error is unmeasurably small for convection-dominated
regimes (high Gr) and masked by backward Euler error for diffusion-dominated
regimes (low Gr). Both self-convergence and Cauchy convergence fail to
produce stable rates. Kept for reference only.

---

Original intent: Self-convergence study isolating the splitting temporal
error. Fixed soil (alpha=1.0, n=1.5, K_s=1e-6), varying column height H
to control Gr = H * alpha:

    Gr=2 (H=2m)    Gr=4 (H=4m)    Gr=10 (H=10m)

Was meant to validate Proposition 4.7:
    ||theta_dt - theta_ref||_{L1} ~ O(min(dt, sqrt(dt/Gr)))

Output: figures/fig_test2_splitting_temporal.pdf + rates table to stdout
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
soil = SoilParams(theta_s=0.40, theta_r=0.04, alpha=1.0,
                  n=1.5, ell=0.5, K_s=1e-6)

configs = [
    {"label": "Gr=2",  "H": 2.0,  "Gr": 2.0},
    {"label": "Gr=4",  "H": 4.0,  "Gr": 4.0},
    {"label": "Gr=10", "H": 10.0, "Gr": 10.0},
]

W = 0.1
maxh = 0.05
t_final = 36000.0  # 10 hours
dt_levels = [3600.0, 1800.0, 900.0, 450.0]  # 1h .. 7.5min
dt_ref = 112.5  # ~1.9 min (fine reference)
p_rain = 10.0 * soil.K_s  # 1e-5 m/s

CACHE = DATA_DIR / "test2_splitting_temporal.npz"
recompute = "--recompute" in sys.argv
all_results = {}

print("Test 2: Splitting temporal convergence (Lie-Trotter)")
print(f"  alpha={soil.alpha}, n={soil.n}, K_s={soil.K_s:.1e}")
print(f"  maxh={maxh}m, t_final={t_final/3600:.0f}h, "
      f"p_rain={p_rain:.1e} m/s")
print(f"  dt_macro levels: {[f'{d:.0f}' for d in dt_levels]}s, "
      f"reference: {dt_ref:.1f}s")
for c in configs:
    print(f"  {c['label']}: H={c['H']:.0f}m")
print()

if not recompute and CACHE.exists():
    print(f"Loading cached data from {CACHE}")
    d = np.load(CACHE)
    all_results = {}
    try:
        for c in configs:
            key = c["label"]
            all_results[key] = {
                "dt": d[f"{key}_dt"],
                "err_theta_L1": d[f"{key}_err_theta_L1"],
                "err_h_L2": d[f"{key}_err_h_L2"],
            }
    except KeyError:
        print("Cache keys mismatch (stale data), recomputing...")
        recompute = True

if recompute or not all_results:
    from nitsche_signorini.runners import run_splitting_temporal_convergence

    all_results = {}
    for c in configs:
        H = c["H"]
        print(f"\n--- {c['label']} (H={H:.0f}m, Gr={c['Gr']:.0f}) ---")
        results = run_splitting_temporal_convergence(
            soil, H, maxh, dt_levels, dt_ref, t_final, p_rain,
            maxerr=1e-6)
        all_results[c["label"]] = {
            "dt": np.array([r["dt"] for r in results]),
            "err_theta_L1": np.array([r["err_theta_L1"] for r in results]),
            "err_h_L2": np.array([r["err_h_L2"] for r in results]),
        }

    CACHE.parent.mkdir(exist_ok=True)
    save_data = {}
    for key, r in all_results.items():
        save_data[f"{key}_dt"] = r["dt"]
        save_data[f"{key}_err_theta_L1"] = r["err_theta_L1"]
        save_data[f"{key}_err_h_L2"] = r["err_h_L2"]
    np.savez(CACHE, **save_data)
    print(f"\nSaved: {CACHE}")

# ── Print convergence tables ──────────────────────────────────────────
for c in configs:
    key = c["label"]
    area = c["H"] * W
    r = all_results[key]
    dt_arr = r["dt"]
    eth = r["err_theta_L1"]
    ehL2 = r["err_h_L2"]

    print(f"\n{key} (H={c['H']:.0f}m, area={area:.1f}):")
    print(f"{'dt':>8} {'||th-ref||_L1':>14} {'rate':>6} "
          f"{'||h-ref||_L2':>14} {'rate':>6}")
    print("-" * 56)
    for i in range(len(dt_arr)):
        rt = (np.log2(eth[i - 1] / eth[i])
              if i > 0 and eth[i] > 0 else float("nan"))
        rh = (np.log2(ehL2[i - 1] / ehL2[i])
              if i > 0 and ehL2[i] > 0 else float("nan"))
        print(f"{dt_arr[i]:>8.0f} {eth[i]:>14.4e} {rt:>6.2f} "
              f"{ehL2[i]:>14.4e} {rh:>6.2f}")

# ── Figure ────────────────────────────────────────────────────────────
fig, ax = plt.subplots(1, 1, figsize=(5.5, 4.5))

colors = ["tab:blue", "tab:red", "tab:green"]
markers = ["o", "s", "D"]

for i, c in enumerate(configs):
    key = c["label"]
    area = c["H"] * W
    r = all_results[key]
    ax.loglog(r["dt"], r["err_theta_L1"] / area,
              marker=markers[i], color=colors[i],
              lw=1.5, ms=7,
              label=rf"$\mathrm{{Gr}}={c['Gr']:.0f}$"
                    rf" ($H={c['H']:.0f}\,\mathrm{{m}}$)")

# Reference slopes anchored at finest-dt point for each curve
dt_line = np.array([dt_levels[0], dt_levels[-1]])
ref_half_plotted = False
ref_one_plotted = False

for i, c in enumerate(configs):
    key = c["label"]
    area = c["H"] * W
    r = all_results[key]
    err_norm = r["err_theta_L1"] / area
    anchor_err = err_norm[-1]  # finest dt
    anchor_dt = r["dt"][-1]

    # O(sqrt(dt)) reference for all curves
    ref_half = anchor_err * np.sqrt(dt_line / anchor_dt)
    lbl_half = (r"$O(\sqrt{\Delta t})$" if not ref_half_plotted
                else None)
    ax.loglog(dt_line, ref_half, color=colors[i], ls="--", lw=1,
              label=lbl_half)
    ref_half_plotted = True

    # O(dt) reference for Euler-dominated curve (Gr=10)
    if c["Gr"] >= 10:
        ref_one = anchor_err * (dt_line / anchor_dt)
        lbl_one = r"$O(\Delta t)$" if not ref_one_plotted else None
        ax.loglog(dt_line, ref_one, color=colors[i], ls=":", lw=1,
                  label=lbl_one)
        ref_one_plotted = True

ax.set_xlabel(r"$\Delta t_{\mathrm{macro}}$ [s]")
ax.set_ylabel(r"$\|\theta_{\Delta t}"
              r" - \theta_{\mathrm{ref}}\|_{L^1}"
              r" / |\Omega|$")
ax.legend()
ax.set_title("Splitting temporal convergence (Lie-Trotter)")

fig.tight_layout()

out = Path(FIGURES_DIR) / "fig_test2_splitting_temporal.pdf"
out.parent.mkdir(exist_ok=True)
fig.savefig(out)
print(f"\nSaved: {out}")
plt.close()
