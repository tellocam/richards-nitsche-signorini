#!/usr/bin/env python
"""
Test 1: Clay baseline -- reproduces Gatti et al. 2024 Figure 6.

Monolithic mixed FEM (RT0 x P0) with Nitsche-Signorini on top boundary.
Clay soil, H=5m, p/K_s=10, dt=5h, t_final=50h.

Output: figures/fig_test1_clay_baseline.pdf
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from _style import DATA_DIR, FIGURES_DIR, apply_style

apply_style()

# ── Parameters (match Gatti et al. Figure 6) ───────────────────────────
H = 5.0          # m
maxh = 0.05      # m
dt = 5 * 3600    # 5 hours [s]
t_final = 50 * 3600  # 50 hours [s]

CACHE = DATA_DIR / "test1_clay_baseline.npz"
recompute = "--recompute" in sys.argv

print(f"Test 1: Clay baseline (Gatti et al. Figure 6)")
print(f"  H={H}m, maxh={maxh}m, dt={dt/3600:.0f}h, "
      f"t_final={t_final/3600:.0f}h, p/K_s=10")
print()

if not recompute and CACHE.exists():
    print(f"Loading cached data from {CACHE}")
    d = np.load(CACHE)
    z = d["z"]
    psi_snaps = list(zip(d["psi_times"], d["psi_profiles"]))
    qz_snaps = list(zip(d["qz_times"], d["qz_profiles"]))
    nits = d["newton_its"]
    n_prof = len(psi_snaps)
    print(f"  {n_prof} snapshots, Newton mean={nits.mean():.1f}, "
          f"max={nits.max()}")
else:
    from nitsche_signorini import CLAY
    from nitsche_signorini.runners import run_monolithic

    soil = CLAY
    p_rain = 10 * soil.K_s
    snap_times = np.arange(dt, t_final + 1, dt)

    z, snaps, info = run_monolithic(
        soil, H, maxh, dt, t_final, p_rain, snap_times, verbose=True)

    nits = info["newton_its"]
    psi_snaps = info["psi_snapshots"]
    qz_snaps = info["qz_snapshots"]
    n_prof = len(snaps)

    print(f"\nCompleted: {info['n_steps']} steps, {info['n_els']} elements")
    print(f"Newton: mean={nits.mean():.1f}, max={nits.max()}")

    CACHE.parent.mkdir(exist_ok=True)
    np.savez(CACHE,
             z=z,
             psi_times=np.array([t for t, _ in psi_snaps]),
             psi_profiles=np.array([p for _, p in psi_snaps]),
             qz_times=np.array([t for t, _ in qz_snaps]),
             qz_profiles=np.array([q for _, q in qz_snaps]),
             newton_its=nits)
    print(f"Saved: {CACHE}")

# ── Figure: psi(z) and q_z(z) profiles ─────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(10, 6))
cmap = plt.cm.viridis

# (a) Pressure head psi(z)
ax = axes[0]
for i, (t_s, psi_s) in enumerate(psi_snaps):
    c = cmap(i / max(len(psi_snaps) - 1, 1))
    label = f"{t_s/3600:.0f}h"
    ax.plot(psi_s, z, color=c, lw=1.2, label=label)
ax.axvline(0, color="grey", ls=":", lw=0.5)
ax.set_xlabel(r"$\psi$ [m]")
ax.set_ylabel("z [m]")
ax.set_title("(a)")
ax.legend(fontsize=6, ncol=2, loc="lower left")

# (b) Vertical Darcy flux q_z(z)
ax = axes[1]
for i, (t_s, qz_s) in enumerate(qz_snaps):
    c = cmap(i / max(len(qz_snaps) - 1, 1))
    label = f"{t_s/3600:.0f}h"
    ax.plot(qz_s, z, color=c, lw=1.2, label=label)
ax.set_xlabel(r"$q_z$ [m/s]")
ax.set_ylabel("z [m]")
ax.set_title("(b)")
ax.legend(fontsize=6, ncol=2, loc="lower left")
fig.tight_layout()

out = Path(FIGURES_DIR) / "fig_test1_clay_baseline.pdf"
out.parent.mkdir(exist_ok=True)
fig.savefig(out)
print(f"\nSaved: {out}")
plt.close()
