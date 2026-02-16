"""
Shared matplotlib styling for CMAME paper figures.
"""

from pathlib import Path

import matplotlib.pyplot as plt

# Resolve to repo root / figures regardless of CWD
REPO_ROOT = Path(__file__).resolve().parent.parent
FIGURES_DIR = str(REPO_ROOT / "figures")
DATA_DIR = REPO_ROOT / "data"


def apply_style():
    """Set consistent matplotlib rcParams for paper-quality figures."""
    plt.rcParams.update({
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 8,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "lines.linewidth": 1.5,
        "lines.markersize": 5,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "font.family": "serif",
        "mathtext.fontset": "cm",
    })
