#!/usr/bin/env python
"""
Test 3b: Splitting on Sand H=5m (was STAGNATED regime).

The monolithic method stagnates because linearized K^{-1} kills the
dynamics. The splitting method resolves the front completely.

Output: figures/fig_test3b_splitting_5m.pdf
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from _style import DATA_DIR, FIGURES_DIR, apply_style

apply_style()

# ── Parameters ─────────────────────────────────────────────────────────
H = 5.0
maxh = 0.1
dt_macro = 120.0
t_final = 36000.0  # 10 hours

CACHE = DATA_DIR / "test3b_splitting_5m.npz"
recompute = "--recompute" in sys.argv

print("Test 3b: Splitting, Sand H=5m")
print(f"  maxh={maxh}m, dt={dt_macro}s, t_final={t_final/3600:.0f}h, p/K_s=10")
print()

if not recompute and CACHE.exists():
    print(f"Loading cached data from {CACHE}")
    d = np.load(CACHE)
    z = d["z"]
    snap_times_arr = d["snap_times"]
    snap_theta = d["snap_theta"]
    snaps = list(zip(snap_times_arr, snap_theta))
    h_snap_times = d["h_snap_times"]
    h_snap_profiles = d["h_snap_profiles"]
    h_snaps = list(zip(h_snap_times, h_snap_profiles))
    nits = d["newton_its"]
    h_max = float(d["h_max"])
    n_steps = int(d["n_steps"])
else:
    from nitsche_signorini import SAND
    from nitsche_signorini.runners import run_splitting

    soil = SAND
    p_rain = 10 * soil.K_s
    snap_times = np.arange(3600, t_final + 1, 3600)

    z, snaps, info = run_splitting(
        soil, H, maxh, dt_macro, t_final, p_rain, snap_times)

    nits = info["newton_its"]
    h_max = info["h_max"]
    n_steps = info["n_steps"]
    h_snaps = info["h_snapshots"]

    snap_times_arr = np.array([t for t, _ in snaps])
    snap_theta = np.array([th for _, th in snaps])
    h_snap_times = np.array([t for t, _ in h_snaps])
    h_snap_profiles = np.array([h for _, h in h_snaps])

    CACHE.parent.mkdir(exist_ok=True)
    np.savez(CACHE, z=z,
             snap_times=snap_times_arr, snap_theta=snap_theta,
             h_snap_times=h_snap_times, h_snap_profiles=h_snap_profiles,
             newton_its=nits, h_max=h_max, n_steps=n_steps)
    print(f"Saved: {CACHE}")

print(f"\nCompleted: {n_steps} steps, Newton mean={nits.mean():.1f}, "
      f"max={nits.max()}, h_max={h_max:.6f}")
h_ok = "PASS" if h_max <= 1e-6 else "FAIL"
print(f"h <= 0 check: {h_ok}")

# ── Figure ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
cmap = plt.cm.viridis
n_snaps = len(snaps)

# (a) theta(z) profiles
ax = axes[0]
for i, (t_s, th_s) in enumerate(snaps):
    c = cmap(i / max(n_snaps - 1, 1))
    label = f"t={t_s/3600:.0f}h" if t_s > 0 else "IC"
    ax.plot(th_s, z, color=c, lw=1.2, label=label)
ax.set_xlabel(r"$\theta$ [-]")
ax.set_ylabel("z [m]")
ax.set_title("(a)")
ax.legend(fontsize=6, ncol=2)

# (b) h(z) profiles
ax = axes[1]
for i, (t_s, h_s) in enumerate(h_snaps):
    c = cmap(i / max(len(h_snaps) - 1, 1))
    label = f"t={t_s/3600:.0f}h" if t_s > 0 else "IC"
    ax.plot(h_s, z, color=c, lw=1.2, label=label)
ax.axvline(0, color="grey", ls=":", lw=0.5, label="h=0")
ax.set_xlabel("h [m]")
ax.set_ylabel("z [m]")
ax.set_title("(b)")
ax.legend(fontsize=6, ncol=2)
fig.tight_layout()

out = Path(FIGURES_DIR) / "fig_test3b_splitting_5m.pdf"
out.parent.mkdir(exist_ok=True)
fig.savefig(out)
print(f"\nSaved: {out}")
plt.close()
