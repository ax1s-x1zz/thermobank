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

from thermobank.data import REACTIONS              # noqa: E402
from thermobank import thermo, screen, uncertainty  # noqa: E402
from thermobank.equilibrium import peq, AMBIENT     # noqa: E402

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
    fig.savefig(os.path.join(OUT, "fig1_dg_vs_T.png"), dpi=200)
    plt.close(fig)
    print("wrote fig1_dg_vs_T.png")


def fig_map():
    rows = screen.run_screen()
    fig, ax = plt.subplots(figsize=(10, 5.8))
    offsets = {
        "glauber": (8, -14), "cacl2_hydrate": (8, 8), "gypsum": (8, -2),
        "cuso4_hydrate": (8, 4), "mgoh2": (8, 4), "mgh2": (8, -14),
        "caoh2": (8, 4), "baoh2": (8, 4), "caco3": (-52, 4),
    }
    for r in rows:
        tstar = r.T_star - 273.15
        color = "#264653" if r.quality == "ref" else "#8ab17d"
        ax.scatter(tstar, r.energy_density, s=60, color=color, zorder=3)
        ax.annotate(r.id, (tstar, r.energy_density),
                    textcoords="offset points", xytext=offsets.get(r.id, (6, 4)),
                    fontsize=8)
    ymax = max(r.energy_density for r in rows)
    for label, lo, hi, c in WINDOWS:
        ax.axvspan(lo - 273.15, hi - 273.15, color=c, alpha=0.18)
        ax.text((lo + hi) / 2 - 273.15, ymax * 1.08, label,
                ha="center", va="top", fontsize=8)
    ax.set_xlabel("Threshold temperature T* (C)")
    ax.set_ylabel("Energy density (kJ/kg)")
    ax.set_title("Storage reactions: temperature requirement vs energy density")
    ax.set_xlim(0, 950)
    ax.set_ylim(0, ymax * 1.18)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig2_map.png"), dpi=200)
    plt.close(fig)
    print("wrote fig2_map.png")


def fig_peq():
    """Equilibrium water-vapour pressure of each dehydration vs ambient humidity."""
    T = np.linspace(275, 420, 300)
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    for rxn in REACTIONS:
        if "H2O(g)" not in rxn["products"]:
            continue
        p = np.array([peq(rxn, t) for t in T])
        ax.plot(T - 273.15, p, lw=1.8, label=rxn["id"])
    for name, p_amb in AMBIENT.items():
        ax.axhline(p_amb, ls=":", lw=1.1, color="gray")
        ax.text(400, p_amb * 1.08, name, va="bottom", ha="right",
                fontsize=7, color="gray")
    ax.axhline(1.0, ls="--", lw=1.0, color="k")
    ax.text(280, 1.05, "p = 1 bar  ->  T*", fontsize=7)
    ax.axvspan(60, 90, color="#2a9d8f", alpha=0.15)
    ax.set_yscale("log")
    ax.set_xlabel("Temperature (C)")
    ax.set_ylabel("Equilibrium water-vapour pressure (bar)")
    ax.set_title("Dehydration proceeds whenever p_eq(T) > ambient humidity")
    ax.set_xlim(0, 420)
    ax.set_ylim(2e-3, 3)
    ax.legend(fontsize=8, ncol=2, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig3_peq.png"), dpi=200)
    plt.close(fig)
    print("wrote fig3_peq.png")


def fig_uncertainty():
    rows = screen.run_screen()
    by_id = {r.id: r for r in rows}
    dists = {d["id"]: d for d in uncertainty.all_distributions()}
    fig, ax = plt.subplots(figsize=(10, 5.8))
    for rid, d in dists.items():
        r = by_id[rid]
        x = d["median"] - 273.15
        xerr = np.array([[x - (d["p05"] - 273.15)], [(d["p95"] - 273.15) - x]])
        ax.errorbar(x, r.energy_density, xerr=xerr, fmt="o", capsize=4,
                    color="#264653" if r.quality == "ref" else "#8ab17d")
        ax.annotate(rid, (x, r.energy_density), textcoords="offset points",
                    xytext=(6, 6), fontsize=8)
    ax.set_xlabel("Threshold temperature T* (C), median with 5-95% interval")
    ax.set_ylabel("Energy density (kJ/kg)")
    ax.set_title("Screening with propagated data uncertainty")
    ax.set_xlim(0, 1000)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig4_uncertainty.png"), dpi=200)
    plt.close(fig)
    print("wrote fig4_uncertainty.png")


if __name__ == "__main__":
    fig_dg_curves()
    fig_map()
    fig_peq()
    fig_uncertainty()
