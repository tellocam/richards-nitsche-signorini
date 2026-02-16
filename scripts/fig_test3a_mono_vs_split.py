#!/usr/bin/env python
"""
Test 3a: Monolithic vs Splitting on Sand H=2m (BLOWUP regime).

The monolithic method struggles with K^{-1} variation (~10 orders);
the Lie-Trotter splitting resolves the wetting front with ~1 Newton
iteration per step.

Output: figures/fig_test3a_mono_vs_split.pdf + comparison table to stdout
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from _style import DATA_DIR, FIGURES_DIR, apply_style

apply_style()

# ── Parameters ─────────────────────────────────────────────────────────
H = 2.0
maxh = 0.05
dt_macro = 60.0
t_final = 10800.0  # 3 hours

CACHE = DATA_DIR / "test3a_mono_vs_split.npz"
recompute = "--recompute" in sys.argv

print("Test 3a: Monolithic vs Splitting, Sand H=2m")
print(f"  maxh={maxh}m, dt={dt_macro}s, t_final={t_final/3600:.0f}h, p/K_s=10")
print()

if not recompute and CACHE.exists():
    print(f"Loading cached data from {CACHE}")
    d = np.load(CACHE)
    z_split = d["z_split"]
    snap_split_times = d["snap_split_times"]
    snap_split_theta = d["snap_split_theta"]
    snaps_split = list(zip(snap_split_times, snap_split_theta))
    nits_m = d["nits_m"]
    nits_s = d["nits_s"]
    mono_n_steps = int(d["mono_n_steps"])
    mono_failed = bool(d["mono_failed"])
    mono_newton_failures = int(d["mono_newton_failures"])
    split_n_steps = int(d["split_n_steps"])
else:
    from nitsche_signorini import SAND
    from nitsche_signorini.runners import run_monolithic, run_splitting

    soil = SAND
    p_rain = 10 * soil.K_s
    snap_times = np.sort(np.unique(np.concatenate([
        np.arange(300, 5401, 300),
        np.arange(6300, t_final + 1, 900),
    ])))

    print("=== Monolithic ===")
    z_mono, snaps_mono, info_mono = run_monolithic(
        soil, H, maxh, dt_macro, t_final, p_rain, snap_times)

    nits_m = info_mono["newton_its"]
    mono_failed = info_mono["failed"]
    mono_n_steps = info_mono["n_steps"]
    mono_newton_failures = info_mono["newton_failures"]
    if mono_failed:
        print(f"FAILED at t={info_mono['failure_time']:.1f}s "
              f"({mono_n_steps} steps, {mono_newton_failures} failures)")
    else:
        total_mono = int(nits_m.sum()) + mono_newton_failures * 50
        print(f"Completed: {mono_n_steps} steps, "
              f"Newton mean={nits_m.mean():.1f}, max={nits_m.max()}, "
              f"{mono_newton_failures} failures, "
              f"total Newton ~ {total_mono}")
    print()

    print("=== Splitting ===")
    z_split, snaps_split, info_split = run_splitting(
        soil, H, maxh, dt_macro, t_final, p_rain, snap_times)

    nits_s = info_split["newton_its"]
    split_n_steps = info_split["n_steps"]
    total_split = int(nits_s.sum())
    print(f"Completed: {split_n_steps} steps, "
          f"Newton mean={nits_s.mean():.1f}, max={nits_s.max()}, "
          f"total Newton ~ {total_split}, h_max={info_split['h_max']:.6f}")

    snap_split_times = np.array([t for t, _ in snaps_split])
    snap_split_theta = np.array([th for _, th in snaps_split])

    CACHE.parent.mkdir(exist_ok=True)
    np.savez(CACHE,
             z_split=z_split,
             snap_split_times=snap_split_times,
             snap_split_theta=snap_split_theta,
             nits_m=nits_m, nits_s=nits_s,
             mono_n_steps=mono_n_steps,
             mono_failed=mono_failed,
             mono_newton_failures=mono_newton_failures,
             split_n_steps=split_n_steps)
    print(f"Saved: {CACHE}")

# ── Comparison table ───────────────────────────────────────────────────
total_split = int(nits_s.sum())
print("\n" + "=" * 70)
print(f"{'Method':<12} {'Steps':>6} {'Newton mean':>12} {'Newton max':>11} "
      f"{'Failures':>9} {'Total':>6}")
print("-" * 70)
if len(nits_m) > 0:
    total_m = int(nits_m.sum()) + mono_newton_failures * 50
    print(f"{'Monolithic':<12} {mono_n_steps:>6} "
          f"{nits_m.mean():>12.1f} {nits_m.max():>11} "
          f"{mono_newton_failures:>9} {total_m:>6}")
print(f"{'Splitting':<12} {split_n_steps:>6} "
      f"{nits_s.mean():>12.1f} {nits_s.max():>11} "
      f"{'0':>9} {total_split:>6}")
if len(nits_m) > 0 and total_split > 0:
    print(f"Cost ratio: {total_m / total_split:.0f}x")
print("=" * 70)

# ── Figure ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
cmap = plt.cm.viridis

# (a) Splitting theta(z) profiles every 10 min up to 90 min
ax = axes[0]
show_times_s = [0] + list(range(600, 5401, 600))
n_show = len(show_times_s)
for i, t_target in enumerate(show_times_s):
    best = min(snaps_split, key=lambda s: abs(s[0] - t_target))
    t_s, th_s = best
    if abs(t_s - t_target) > 120:
        continue
    c = cmap(i / (n_show - 1))
    label = "IC" if t_s == 0 else f"t={t_s / 60:.0f}min"
    ls = ":" if t_s == 0 else "-"
    ax.plot(th_s, z_split, ls, color=c, lw=1.5, label=label)
ax.set_xlabel(r"$\theta$ [-]")
ax.set_ylabel("z [m]")
ax.set_title("(a)")
ax.legend(fontsize=6, ncol=2)

# (b) Newton iterations per step
ax = axes[1]
if len(nits_m) > 0:
    mono_steps = np.arange(1, len(nits_m) + 1)
    ax.plot(mono_steps, nits_m, "r.-", ms=3, lw=0.8,
            label=f"Monolithic ({mono_newton_failures} failures)")
split_steps = np.arange(1, len(nits_s) + 1)
ax.plot(split_steps, nits_s, "b.-", ms=2, lw=0.5, alpha=0.7,
        label=f"Splitting (mean={nits_s.mean():.1f})")
ax.set_xlabel("Time step")
ax.set_ylabel("Newton iterations")
ax.set_title("(b)")
ax.legend(fontsize=8)
fig.tight_layout()

out = Path(FIGURES_DIR) / "fig_test3a_mono_vs_split.pdf"
out.parent.mkdir(exist_ok=True)
fig.savefig(out)
print(f"\nSaved: {out}")
plt.close()
