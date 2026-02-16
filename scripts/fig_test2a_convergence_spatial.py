#!/usr/bin/env python
"""
Test 2a: Spatial convergence via manufactured solution.

O(h) for RT0 x P0, four parameter surrogates spanning the degeneracy
transition. Stationary problem with mesh.Refine() for clean factor-2
refinement.

Output: figures/fig_test2a_convergence_spatial.pdf + rates table to stdout
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from _style import DATA_DIR, FIGURES_DIR, apply_style

from nitsche_signorini import SoilParams

apply_style()

# ── Surrogates spanning the degeneracy transition ──────────────────────
surrogates = {
    "D1": SoilParams(theta_s=0.4, theta_r=0.05, alpha=0.5, n=1.3,
                     ell=0.5, K_s=1e-4),
    "D2": SoilParams(theta_s=0.4, theta_r=0.05, alpha=1.0, n=1.8,
                     ell=0.5, K_s=1e-4),
    "C1": SoilParams(theta_s=0.4, theta_r=0.05, alpha=2.0, n=2.5,
                     ell=0.5, K_s=1e-4),
    "C2": SoilParams(theta_s=0.4, theta_r=0.05, alpha=2.0, n=3.5,
                     ell=0.5, K_s=1e-4),
}

n_levels = 5
norms = ["err_h", "err_q", "err_theta_L1", "err_theta_L2"]
CACHE = DATA_DIR / "test2a_convergence_spatial.npz"
recompute = "--recompute" in sys.argv

print("Test 2a: Spatial convergence (manufactured solution)")
print(f"{'Label':>5} {'n':>5} {'alpha':>6} {'beta':>6} {'Regime'}")
print("-" * 45)
for label, s in surrogates.items():
    beta = s.ell + s.n / (s.n - 1)
    regime = "diffusion" if beta > 2.5 else "convection"
    print(f"{label:>5} {s.n:>5.1f} {s.alpha:>6.1f} {beta:>6.2f} {regime}")
print()

if not recompute and CACHE.exists():
    print(f"Loading cached data from {CACHE}")
    d = np.load(CACHE)
    results = {}
    for label in surrogates:
        results[label] = {k: d[f"{label}_{k}"] for k in norms + ["ne"]}
else:
    from nitsche_signorini.runners import run_spatial_convergence

    results = run_spatial_convergence(surrogates, n_levels=n_levels)

    CACHE.parent.mkdir(exist_ok=True)
    data = {}
    for label, r in results.items():
        for k in norms + ["ne"]:
            data[f"{label}_{k}"] = np.array(r[k])
    np.savez(CACHE, **data)
    print(f"Saved: {CACHE}")

# ── Print convergence tables ───────────────────────────────────────────
for label in surrogates:
    r = results[label]
    beta = surrogates[label].ell + surrogates[label].n / (surrogates[label].n - 1)
    print(f"\n{label} (beta={beta:.2f}):")
    print(f"{'lvl':>4} {'ne':>6} {'||h||':>10} {'r':>5} "
          f"{'||q||':>10} {'r':>5} "
          f"{'||th||_L1':>10} {'r':>5} "
          f"{'||th||_L2':>10} {'r':>5}")
    print("-" * 75)
    for i in range(n_levels):
        rates = []
        for k in norms:
            if i > 0 and r[k][i] > 0 and r[k][i - 1] > 0:
                rates.append(np.log2(r[k][i - 1] / r[k][i]))
            else:
                rates.append(float("nan"))
        print(f"{i:>4} {r['ne'][i]:>6.0f} "
              f"{r['err_h'][i]:>10.2e} {rates[0]:>5.2f} "
              f"{r['err_q'][i]:>10.2e} {rates[1]:>5.2f} "
              f"{r['err_theta_L1'][i]:>10.2e} {rates[2]:>5.2f} "
              f"{r['err_theta_L2'][i]:>10.2e} {rates[3]:>5.2f}")

# ── Figure: log-log convergence ────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(10, 8))

levels = np.arange(n_levels)
markers = {"D1": "o-", "D2": "s-", "C1": "^--", "C2": "v--"}
colors = {"D1": "tab:blue", "D2": "tab:orange",
          "C1": "tab:green", "C2": "tab:red"}
norm_labels = {
    "err_h": r"$\|h - h_{\mathrm{ex}}\|_{L^2}$",
    "err_q": r"$\|q - q_{\mathrm{ex}}\|_{L^2}$",
    "err_theta_L1": r"$\|\theta - \theta_{\mathrm{ex}}\|_{L^1}$",
    "err_theta_L2": r"$\|\theta - \theta_{\mathrm{ex}}\|_{L^2}$",
}
titles = ["(a)", "(b)", "(c)", "(d)"]

for idx, (norm_key, ax) in enumerate(zip(norms, axes.flat)):
    for label in surrogates:
        ax.semilogy(levels, results[label][norm_key],
                    markers[label], color=colors[label],
                    label=label, lw=1.5, ms=6)
    ref_anchor = results["D2"][norm_key][0]
    ref_line = ref_anchor * 2.0 ** (-levels)
    ax.semilogy(levels, ref_line, "k:", lw=1, label="O(h)")
    ax.set_xlabel("Refinement level")
    ax.set_ylabel(norm_labels[norm_key])
    ax.set_title(f"{titles[idx]} {norm_labels[norm_key]}")
    ax.legend(fontsize=7)
    ax.set_xticks(levels)

fig.suptitle("Test 2a: Spatial convergence (manufactured solution)",
             fontsize=12)
fig.tight_layout()

out = Path(FIGURES_DIR) / "fig_test2a_convergence_spatial.pdf"
out.parent.mkdir(exist_ok=True)
fig.savefig(out)
print(f"\nSaved: {out}")
plt.close()
