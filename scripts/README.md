# Reproducibility scripts

Figure-generation scripts for the CMAME paper. Each `fig_*.py` produces one
PDF in `figures/` and prints convergence tables to stdout.

## Running

```bash
uv run python scripts/fig_testN_name.py           # use cached data
uv run python scripts/fig_testN_name.py --recompute  # recompute from scratch
```

Cached `.npz` files are stored in `data/`.

## File structure

```
scripts/
  fig_test1_clay_baseline.py          # Test 1: clay + monolithic Nitsche-Signorini
  fig_test2a_convergence_spatial.py   # Test 2a: spatial convergence (manufactured soln)
  fig_test2b_convergence_temporal.py  # Test 2b: monolithic temporal convergence (BE)
  fig_test3a_mono_vs_split.py         # Test 3a: monolithic vs splitting (sand)
  fig_test3b_splitting_5m.py          # Test 3b: splitting on 5m sand column
  fig_test3c_front_speed.py           # Test 3c: front speed vs Rankine-Hugoniot ODE
  fig_test4_rain_switching.py         # Test 4: time-varying rain, Signorini switching
  _style.py                           # Shared matplotlib rcParams (helper)
  fig_test2_splitting_temporal.py     # NOT USED -- splitting temporal convergence (failed)
  _diag_cauchy.py                     # NOT USED -- Cauchy diagnostic (confirmed failure)
```

## Script inventory

### Active (used in paper)

| Script | Figure | Description |
|--------|--------|-------------|
| `fig_test1_clay_baseline.py` | `fig_test1_clay_baseline.pdf` | Clay + monolithic Nitsche-Signorini, reproduces Gatti et al. |
| `fig_test2a_convergence_spatial.py` | `fig_test2a_convergence_spatial.pdf` | Spatial convergence via manufactured solution |
| `fig_test2b_convergence_temporal.py` | `fig_test2b_convergence_temporal.pdf` | Monolithic backward Euler temporal convergence |
| `fig_test3a_mono_vs_split.py` | `fig_test3a_mono_vs_split.pdf` | Monolithic fails for sand; splitting resolves it |
| `fig_test3b_splitting_5m.py` | `fig_test3b_splitting_5m.pdf` | Splitting on 5m sand column (was STAGNATED) |
| `fig_test3c_front_speed.py` | `fig_test3c_front_speed.pdf` | Numerical front speed vs Rankine-Hugoniot ODE |
| `fig_test4_rain_switching.py` | `fig_test4_rain_switching.pdf` | Time-varying rain with bidirectional Signorini switching |

### Not used (retired experiments)

| Script | Why retired |
|--------|------------|
| `fig_test2_splitting_temporal.py` | Splitting temporal convergence rate unmeasurable (see below) |
| `_diag_cauchy.py` | Diagnostic confirming the above; Cauchy differences are flat |

---

## Retired: splitting temporal convergence

**Do not invest time trying to make this experiment work.** We conducted a
thorough investigation (Feb 2026) and found a fundamental catch-22 that
prevents clean demonstration of the splitting temporal convergence rate
from Proposition 4.7.

### What we tried

Self-convergence study: fixed mesh (maxh=0.05), single soil
(alpha=1.0, n=1.5, K\_s=1e-6), three column heights H=2,4,10m giving
Gr=2,4,10. dt\_macro in [450, 3600]s, reference at dt\_ref=112.5s.
Goal: show O(sqrt(dt/Gr)) for splitting-dominated regime and
O(dt) for Euler-dominated regime.

### What we found

1. **Two-term error fit (A\*sqrt(dt) + B\*dt) gives B < 0 for all
   configs.** The backward Euler term is negligible. No Euler-dominated
   regime is visible -- all configs are splitting-dominated.

2. **Self-convergence rates are unreliable.** With dt\_ref/dt\_finest = 0.25,
   the reference error contaminates the fine-dt measurements by ~50%
   (for rate 0.5). Observed rates for Gr=2 were 0.30, 0.48, 0.80 --
   drifting upward due to reference contamination, not reflecting a true
   power law.

3. **Cauchy convergence (no reference) reveals flat differences.** The
   `_diag_cauchy.py` diagnostic computes ||theta\_dt - theta\_{dt/2}||\_L1
   between consecutive dt levels. For Gr=2, these differences are
   approximately constant (~1.5e-4) across all dt pairs. Cauchy rates:
   0.40, -0.07, -0.06. The splitting temporal error is too small to
   measure.

### Root cause: catch-22

The proposition bounds the splitting error by O(min(dt, sqrt(dt/Gr))).

- **High Gr (convection-dominated):** The diffusion correction is nearly
  identity. The Godunov step does all the work. Changing dt\_macro barely
  affects the solution. The splitting temporal error is unmeasurably small.
  This is *why the method works so well* (Newton ~2 its, h <= 0), but it
  also means we cannot measure the convergence rate.

- **Low Gr (diffusion-dominated):** The splitting error sqrt(dt/Gr)
  exceeds dt, so the bound reduces to O(dt) -- just backward Euler.
  The splitting-specific rate (0.5) is invisible because the Euler
  error dominates.

A narrow sweet spot may exist, but finding it requires knowing
C\_split/C\_BE a priori, which depends on the solution profile. We could
not determine this ratio theoretically or fit it reliably from data
(B < 0 in all cases).

### Conclusion

The splitting's value is **robustness** (converges where monolithic fails),
not a specific temporal convergence rate. This is demonstrated by the
monolithic-vs-splitting comparison (Test 3a) and the sand failure phase
diagram, not by a convergence rate study.
