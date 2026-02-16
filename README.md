# richards-nitsche-signorini

Reproducibility code for:

> **Operator-split mixed finite elements for Richards equation with Nitsche-Signorini seepage conditions**
>
> Camilo Tello Fachin, Federico Gatti, Wouter Tonnon
>
> *Computer Methods in Applied Mechanics and Engineering* (submitted)

## Quick start

### Prerequisites

Install [uv](https://docs.astral.sh/uv/getting-started/installation/) (Python package manager):

```bash
# macOS
brew install uv

# Linux
curl -LsSf https://astral.sh/uv/install.sh | sh
```

### Install & run

```bash
git clone https://github.com/tellocam/richards-nitsche-signorini.git
cd richards-nitsche-signorini
uv sync                    # install dependencies
uv run python scripts/fig_test1_clay_baseline.py   # generate one figure
```

The generated figure is saved to `figures/`. Pre-computed data in `data/`
is used by default; pass `--recompute` to re-run the simulation from scratch.

## Repository structure

```
richards-nitsche-signorini/
├── pyproject.toml                          # Project metadata & dependencies
├── uv.lock                                 # Locked dependency versions
├── Makefile                                # Convenience targets (figures, notebooks, clean)
├── LICENSE                                 # MIT license
│
├── src/nitsche_signorini/                  # Python package (installed via uv sync)
│   ├── __init__.py                         # Public API exports
│   ├── common.py                           # CoefficientFunction helpers
│   ├── constitutive.py                     # Sharp van Genuchten-Mualem model
│   ├── materials.py                        # Soil parameter presets (CLAY, SAND, ...)
│   ├── splitting.py                        # Lie-Trotter operator splitting core
│   └── runners.py                          # Time-stepping runners for all test cases
│
├── scripts/                                # Figure-generation scripts (one per paper figure)
│   ├── README.md                           # Script inventory & notes on retired experiments
│   ├── _style.py                           # Shared matplotlib rcParams & path constants
│   ├── fig_test1_clay_baseline.py          # Fig. 1: clay + monolithic Nitsche-Signorini
│   ├── fig_test2a_convergence_spatial.py   # Fig. 2: spatial convergence (manufactured soln)
│   ├── fig_test2b_convergence_temporal.py  # Table:  temporal convergence (backward Euler)
│   ├── fig_test3a_mono_vs_split.py         # Fig. 3: monolithic fails, splitting fixes it
│   ├── fig_test3b_splitting_5m.py          # Fig. 4: splitting on 5 m sand column
│   ├── fig_test3c_front_speed.py           # Fig. 5: front speed vs Rankine-Hugoniot ODE
│   ├── fig_test4_rain_switching.py         # Fig. 6: bidirectional Signorini switching
│   ├── fig_test2_splitting_temporal.py     # (retired) splitting temporal convergence
│   └── _diag_cauchy.py                     # (retired) Cauchy diagnostic
│
├── notebooks/                              # Interactive Jupyter notebooks
│   ├── 01_clay_baseline.ipynb              # Clay baseline walkthrough
│   ├── 02_convergence_study.ipynb          # Convergence rate exploration
│   ├── 03_sand_splitting.ipynb             # Operator splitting on sand
│   ├── 04_signorini_switching.ipynb        # Dynamic boundary switching
│   └── 05_failed_approaches.ipynb          # Why standard methods fail for sand
│
├── data/                                   # Pre-computed results (.npz cache, 132 KB)
│   └── *.npz                               # One file per test case
│
└── figures/                                # Generated PDF figures (not committed)
    └── *.pdf                               # One file per script
```

## Figure index

| Script | Paper element | Description | ~Runtime |
|--------|--------------|-------------|----------|
| `fig_test1_clay_baseline.py` | Fig. 1 | Clay baseline: psi(z) and q_z(z) profiles | 1 min |
| `fig_test2a_convergence_spatial.py` | Fig. 2 + Table | Spatial convergence (manufactured solution) | 2 min |
| `fig_test2b_convergence_temporal.py` | Table | Temporal convergence (backward Euler) | 3 min |
| `fig_test3a_mono_vs_split.py` | Fig. 3 + Table | Monolithic vs splitting, sand H=2m | 5 min |
| `fig_test3b_splitting_5m.py` | Fig. 4 | Splitting on sand H=5m (was STAGNATED) | 5 min |
| `fig_test3c_front_speed.py` | Fig. 5 | Front speed validation vs Rankine-Hugoniot | 8 min |
| `fig_test4_rain_switching.py` | Fig. 6 | Bidirectional Signorini switching | 5 min |

## Generate all figures

```bash
make figures
# or: for f in scripts/fig_test*.py; do uv run python "$f"; done
```

Total runtime: approximately 30 minutes on a modern workstation (instant with cached data).

## Notebooks

Interactive Jupyter notebooks for exploration and education.

| Notebook | Description |
|----------|-------------|
| `01_clay_baseline.ipynb` | Clay baseline (monolithic method), reproduces Gatti et al. Figure 6 |
| `02_convergence_study.ipynb` | Spatial and temporal convergence rates |
| `03_sand_splitting.ipynb` | Operator splitting on sand (the stabilised method) |
| `04_signorini_switching.ipynb` | Dynamic Signorini switching under time-varying rain |
| `05_failed_approaches.ipynb` | **Failed approaches**: why standard methods fail for sand |

Notebook 05 documents the systematic investigation of alternative methods
(broken RT0, HDG, TVB limiter, linearised schemes) and demonstrates that
all fail due to the same root cause: $K^{-1}$ blow-up spanning ~7 orders
of magnitude for sand. This motivates the operator splitting approach.

```bash
uv sync --extra notebooks   # install Jupyter (first time only)
uv run jupyter lab notebooks/
# or: make notebooks
```

## Library overview

The `src/nitsche_signorini/` package provides:

| Module | Description |
|--------|-------------|
| `common.py` | CoefficientFunction helpers (`clamp_cf`, `pow_cf`, `ez_cf`) |
| `constitutive.py` | Sharp van Genuchten-Mualem constitutive law (no regularization) |
| `materials.py` | Soil parameter presets (CLAY, SAND, SANDY_LOAM) |
| `splitting.py` | Lie-Trotter operator splitting: Godunov upwind + mixed FEM diffusion + Nitsche-Signorini BCs |
| `runners.py` | Time-stepping runners for all test cases |

## Solver parameters

| Parameter | Clay (Test 1) | Sand splitting (Tests 3-4) |
|-----------|--------------|---------------------------|
| FE space | RT0 x P0 | RT0 x P0 |
| Solver | Newton (monolithic) | Newton (diffusion step only) |
| Convection | implicit (gravity in flux eq.) | explicit Godunov subcycled |
| Nitsche gamma_0 | 1e-10 | 1e-10 |
| Newton tol | 1e-6 | 1e-3 (splitting), 5e-3 (rain switching) |
| Linear solver | PARDISO | PARDISO |

## Dependencies

- Python >= 3.12
- [NGSolve](https://ngsolve.org/) (finite element library, includes PARDISO solver)
- NumPy, SciPy, matplotlib

All dependencies are managed via `pyproject.toml` and installable with `uv sync`.
Works on Linux (x86_64) and macOS (Intel & Apple Silicon).

## License

MIT. See [LICENSE](LICENSE).
