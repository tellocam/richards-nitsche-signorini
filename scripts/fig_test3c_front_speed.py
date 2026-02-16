#!/usr/bin/env python
"""
Test 3c: Front speed validation.

Pure Godunov (convection only) should track the Rankine-Hugoniot prediction.
Full splitting (Godunov + diffusion) produces faster fronts due to the
diffusion step's water table BC.

Output: figures/fig_test3c_front_speed.pdf
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from _style import DATA_DIR, FIGURES_DIR, apply_style

from nitsche_signorini import SAND
from nitsche_signorini.constitutive import theta_from_h_np
from nitsche_signorini.runners import build_K_table

apply_style()

soil = SAND
p_rain = 10 * soil.K_s

# ── R-H speeds (cheap -- always recompute) ────────────────────────────
theta_tab, K_tab = build_K_table(soil, n_pts=50000)

configs = [
    {"H": 2.0, "maxh": 0.05, "dt_macro": 60.0, "t_final": 10800.0},
    {"H": 5.0, "maxh": 0.1, "dt_macro": 120.0, "t_final": 36000.0},
]

for cfg in configs:
    H = cfg["H"]
    theta_ic = float(theta_from_h_np(np.array([-H]), soil)[0])
    K_ic = float(np.interp(theta_ic, theta_tab, K_tab))
    s_rh = (soil.K_s - K_ic) / (soil.theta_s - theta_ic)
    cfg["s_rh"] = s_rh
    cfg["snaps"] = np.arange(
        cfg["dt_macro"] * 15, cfg["t_final"] + 1,
        cfg["dt_macro"] * 15)

CACHE = DATA_DIR / "test3c_front_speed.npz"
recompute = "--recompute" in sys.argv

print("Test 3c: Front speed validation")
for i, cfg in enumerate(configs):
    print(f"  H={cfg['H']}m: s_RH={cfg['s_rh']:.4e} m/s")
print()

if not recompute and CACHE.exists():
    print(f"Loading cached data from {CACHE}")
    d = np.load(CACHE)
    results = []
    for j, cfg in enumerate(configs):
        results.append({
            "cfg": cfg,
            "fl_g": d[f"fl_g_{j}"],
            "fl_s": d[f"fl_s_{j}"],
        })
else:
    from nitsche_signorini.runners import run_godunov_only, run_splitting

    results = []
    for cfg in configs:
        H = cfg["H"]
        print(f"--- H={H}m ---")

        print("  Godunov only ...", end=" ", flush=True)
        z_g, snaps_g, fl_g = run_godunov_only(
            soil, H, cfg["maxh"], cfg["dt_macro"], cfg["t_final"],
            cfg["snaps"])
        print(f"{len(fl_g)-1} front samples")

        print("  Splitting (dry IC) ...", end=" ", flush=True)
        z_s, snaps_s, info_s = run_splitting(
            soil, H, cfg["maxh"], cfg["dt_macro"], cfg["t_final"],
            p_rain, cfg["snaps"], uniform_dry_ic=True, verbose=False)
        fl_s = info_s["front_log"]
        print(f"{info_s['n_steps']} steps, "
              f"Newton mean={info_s['newton_its'].mean():.1f}")

        results.append({"cfg": cfg, "fl_g": fl_g, "fl_s": fl_s})

    CACHE.parent.mkdir(exist_ok=True)
    data = {}
    for j, r in enumerate(results):
        data[f"fl_g_{j}"] = np.asarray(r["fl_g"])
        data[f"fl_s_{j}"] = np.asarray(r["fl_s"])
    np.savez(CACHE, **data)
    print(f"Saved: {CACHE}")

# ── Deviation table ────────────────────────────────────────────────────
print("\nFront deviations from R-H prediction")
print("=" * 50)
print(f"{'':<10} {'Godunov':<12} {'Splitting':<12}")
print("-" * 50)
for r in results:
    H = r["cfg"]["H"]
    s_rh = r["cfg"]["s_rh"]
    dev_g = max(abs(z - max(H - s_rh * t, 0.0))
                for t, z in r["fl_g"])
    dev_s = max(abs(z - max(H - s_rh * t, 0.0))
                for t, z in r["fl_s"])
    print(f"H={H:.0f}m     {dev_g:<12.3f}m {dev_s:<12.3f}m")

# ── Figure ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

for idx, (r, ax) in enumerate(zip(results, axes)):
    H = r["cfg"]["H"]
    s_rh = r["cfg"]["s_rh"]
    t_final = r["cfg"]["t_final"]

    t_rh = np.linspace(0, t_final, 200)
    ax.plot(t_rh / 3600, np.maximum(H - s_rh * t_rh, 0.0),
            "k-", lw=2, label="R-H prediction")
    ax.plot(r["fl_g"][:, 0] / 3600, r["fl_g"][:, 1],
            "g^-", ms=5, lw=1, label="Godunov only")
    ax.plot(r["fl_s"][:, 0] / 3600, r["fl_s"][:, 1],
            "bo-", ms=4, lw=1, label="Splitting")
    ax.set_xlabel("Time [h]")
    ax.set_ylabel("Front position z [m]")
    ax.set_title(f"({'a' if idx == 0 else 'b'}) H = {H:.0f} m")
    ax.legend(fontsize=7)

fig.suptitle("Test 3c: Front speed (uniform dry IC)", fontsize=11)
fig.tight_layout()

out = Path(FIGURES_DIR) / "fig_test3c_front_speed.pdf"
out.parent.mkdir(exist_ok=True)
fig.savefig(out)
print(f"\nSaved: {out}")
plt.close()
