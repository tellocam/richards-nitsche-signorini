#!/usr/bin/env python
"""
Test 4: Bidirectional Signorini switching under time-varying rain.

Full rain cycle: ramp up [0,0.5h], hold [0.5h,1.5h], ramp down [1.5h,2h],
drainage [2h,5h]. Demonstrates infiltration -> seepage -> infiltration
switching with dynamic Godunov ghost.

Output: figures/fig_test4_rain_switching.pdf
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from _style import DATA_DIR, FIGURES_DIR, apply_style

apply_style()

# ── Parameters ─────────────────────────────────────────────────────────
H = 2.0
maxh = 0.05
dt_macro = 60.0
t_final = 18000.0  # 5 hours

CACHE = DATA_DIR / "test4_rain_switching.npz"
recompute = "--recompute" in sys.argv

print("Test 4: Bidirectional Signorini switching")
print(f"  H={H}m, maxh={maxh}m, dt={dt_macro}s, t_final={t_final/3600:.0f}h")
print(f"  Rain schedule: ramp up [0,0.5h], hold [0.5h,1.5h], "
      f"ramp down [1.5h,2h], off [2h,5h]")
print()

if not recompute and CACHE.exists():
    print(f"Loading cached data from {CACHE}")
    d = np.load(CACHE)
    z = d["z"]
    snap_times_arr = d["snap_times"]
    snap_theta = d["snap_theta"]
    snaps = list(zip(snap_times_arr, snap_theta))
    nits = d["newton_its"]
    t_history = d["t_history"]
    p_rain_history = d["p_rain_history"]
    h_top_history = d["h_top_history"]
    mode_int = d["signorini_mode_int"]
    modes = ["seepage" if m == 0 else "infiltration" for m in mode_int]
    n_steps = int(d["n_steps"])
    h_max = float(d["h_max"])
    K_s = float(d["K_s"])
else:
    from nitsche_signorini import SAND
    from nitsche_signorini.runners import make_piecewise_linear, run_splitting_ramp

    soil = SAND
    K_s = soil.K_s
    p_high = 10.0 * K_s
    p_func = make_piecewise_linear([
        (0, 0), (1800, p_high), (5400, p_high), (7200, 0), (99999, 0)])

    snap_times = np.array(
        [900, 1800, 3600, 5400, 6300, 7200, 10800, 14400, 18000],
        dtype=float)

    z, snaps, info = run_splitting_ramp(
        soil, H, maxh, dt_macro, t_final, p_func, snap_times)

    nits = info["newton_its"]
    n_steps = info["n_steps"]
    h_max = info["h_max"]
    t_history = np.array(info["t_history"])
    p_rain_history = np.array(info["p_rain_history"])
    h_top_history = np.array(info["h_top_history"])
    modes = info["signorini_mode"]
    mode_int = np.array([0 if m == "seepage" else 1 for m in modes])

    snap_times_arr = np.array([t for t, _ in snaps])
    snap_theta = np.array([th for _, th in snaps])

    CACHE.parent.mkdir(exist_ok=True)
    np.savez(CACHE, z=z,
             snap_times=snap_times_arr, snap_theta=snap_theta,
             newton_its=nits, t_history=t_history,
             p_rain_history=p_rain_history,
             h_top_history=h_top_history,
             signorini_mode_int=mode_int,
             n_steps=n_steps, h_max=h_max, K_s=K_s)
    print(f"Saved: {CACHE}")

print(f"\nCompleted: {n_steps} steps, "
      f"Newton mean={nits.mean():.1f}, max={nits.max()}, "
      f"h_max={h_max:.6f}")

n_seepage = sum(1 for m in modes if m == "seepage")
n_infilt = sum(1 for m in modes if m == "infiltration")
print(f"Mode: {n_seepage} seepage, {n_infilt} infiltration steps")
h_ok = "PASS" if h_max <= 1e-6 else "FAIL"
print(f"h <= 0 check: {h_ok}")

# ── Figure: 2x2 panel ─────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
cmap = plt.cm.viridis
ramp_up = (0.0, 0.5)
ramp_dn = (1.5, 2.0)

# (a) Rain schedule + Newton iterations
ax = axes[0, 0]
t_hist = t_history / 3600
p_hist = p_rain_history
ax.plot(t_hist, p_hist / K_s, "b-", lw=1.5, label="p/K_s")
ax.axvspan(*ramp_up, alpha=0.15, color="green", label="ramp up")
ax.axvspan(*ramp_dn, alpha=0.15, color="red", label="ramp down")
ax.set_xlabel("Time [h]")
ax.set_ylabel("p / K_s")
ax2 = ax.twinx()
step_t = t_history[1:] / 3600
ax2.plot(step_t, nits, ".", ms=2, color="orange", alpha=0.5)
ax2.set_ylabel("Newton its", color="orange")
ax.set_title("(a)")
ax.legend(fontsize=7, loc="upper right")

# (b) theta(z) profiles
ax = axes[0, 1]
n_snaps = len(snaps)
for i, (t_s, th_s) in enumerate(snaps):
    c = cmap(i / max(n_snaps - 1, 1))
    label = f"t={t_s/3600:.2f}h" if t_s > 0 else "IC"
    ax.plot(th_s, z, color=c, lw=1.2, label=label)
ax.set_xlabel(r"$\theta$ [-]")
ax.set_ylabel("z [m]")
ax.set_title("(b)")
ax.legend(fontsize=5)

# (c) h_top(t)
ax = axes[1, 0]
t_steps = t_history[1:] / 3600
ax.plot(t_steps, h_top_history, lw=1)
ax.axhline(0, color="grey", ls=":", lw=0.5)
ax.axvspan(*ramp_up, alpha=0.15, color="green", label="ramp up")
ax.axvspan(*ramp_dn, alpha=0.15, color="red", label="ramp down")
ax.set_xlabel("Time [h]")
ax.set_ylabel("h_top [m]")
ax.set_title("(c)")
ax.legend(fontsize=7)

# (d) Signorini mode timeline
ax = axes[1, 1]
mode_num = np.array([0 if m == "seepage" else 1 for m in modes])
changes = np.where(np.diff(mode_num) != 0)[0] + 1
boundaries = np.concatenate([[0], changes, [len(mode_num)]])
for k in range(len(boundaries) - 1):
    si = boundaries[k]
    ei = boundaries[k + 1] - 1
    t0 = 0.0 if si == 0 else t_steps[si]
    t1 = t_steps[min(ei, len(t_steps) - 1)]
    color = "steelblue" if mode_num[si] == 0 else "sandybrown"
    ax.axvspan(t0, t1, alpha=0.4, color=color)
ax.axvline(ramp_up[1], color="green", ls="--", lw=1.5)
ax.axvline(ramp_dn[0], color="red", ls="--", lw=1.5)
ax.axvline(ramp_dn[1], color="red", ls="--", lw=1.5)
ax.set_xlim(t_steps[0], t_steps[-1])
ax.set_ylim(0, 1)
ax.set_yticks([])
ax.set_xlabel("Time [h]")
ax.set_title("(d)")
legend_elements = [
    Patch(facecolor="steelblue", alpha=0.4, label="seepage"),
    Patch(facecolor="sandybrown", alpha=0.4, label="infiltration"),
    plt.Line2D([0], [0], color="green", ls="--", lw=1.5,
               label="ramp up end"),
    plt.Line2D([0], [0], color="red", ls="--", lw=1.5,
               label="ramp down"),
]
ax.legend(handles=legend_elements, fontsize=7, loc="center right")

fig.tight_layout()

out = Path(FIGURES_DIR) / "fig_test4_rain_switching.pdf"
out.parent.mkdir(exist_ok=True)
fig.savefig(out)
print(f"\nSaved: {out}")
plt.close()
