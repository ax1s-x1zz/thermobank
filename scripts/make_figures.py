"""Generate figures for the thermobank report.

Outputs PNG files into results/figures/.
"""

from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from thermobank.data import REACTIONS          # noqa: E402
from thermobank import thermo, screen           # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "figures")
os.makedirs(OUT, exist_ok=True)

WINDOWS = [
    ("data center\n60-90 C", 333.15, 363.15, "#2a9d8f"),
    ("industrial low\n150-200 C", 423.15, 473.15, "#e9c46a"),
    ("industrial high\n300-500 C", 573.15, 773.15, "#e76f51"),
]


def fig_dg_curves():
    T = np.linspace(280, 1200, 400)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for rxn in REACTIONS:
        g = np.array([thermo.delta_g(rxn, t) / 1000.0 for t in T])
        ax.plot(T - 273.15, g, lw=1.8, label=rxn["id"])
    for _, lo, hi, c in WINDOWS:
        ax.axvspan(lo - 273.15, hi - 273.15, color=c, alpha=0.18)
    ax.axhline(0, color="k", lw=0.8, ls="--")
    ax.set_xlabel("Temperature (C)")
    ax.set_ylabel("dG (kJ/mol)")
    ax.set_title("Gibbs energy of dehydration/decomposition vs temperature")
    ax.set_ylim(-120, 260)
    ax.set_xlim(0, 950)
    ax.legend(fontsize=8, ncol=3)
    fig.tight_layout()
    p = os.path.join(OUT, "fig1_dg_vs_T.png")
    fig.savefig(p, dpi=200)
    plt.close(fig)
    print("wrote", p)


def fig_map():
    rows = screen.run_screen()
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for r in rows:
        tstar = r.T_star - 273.15
        color = "#264653" if r.quality == "ref" else "#8ab17d"
        ax.scatter(tstar, r.energy_density, s=60, color=color, zorder=3)
        ax.annotate(r.id, (tstar, r.energy_density),
                    textcoords="offset points", xytext=(6, 4), fontsize=8)
    for label, lo, hi, c in WINDOWS:
        ax.axvspan(lo - 273.15, hi - 273.15, color=c, alpha=0.18)
        ax.text((lo + hi) / 2 - 273.15, ax.get_ylim()[1] * 0.96,
                label, ha="center", va="top", fontsize=8)
    ax.set_xlabel("Threshold temperature T* (C)")
    ax.set_ylabel("Energy density (kJ/kg)")
    ax.set_title("Storage reactions: temperature requirement vs energy density")
    ax.set_xlim(0, 950)
    fig.tight_layout()
    p = os.path.join(OUT, "fig2_map.png")
    fig.savefig(p, dpi=200)
    plt.close(fig)
    print("wrote", p)


if __name__ == "__main__":
    fig_dg_curves()
    fig_map()
