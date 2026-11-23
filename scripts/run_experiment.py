"""Run the computational kinetics experiment and write results + figure."""

from __future__ import annotations

import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from thermobank import simulation as sim          # noqa: E402
from thermobank.data import REACTIONS             # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT_FIG = os.path.join(ROOT, "results", "figures")
OUT_CSV = os.path.join(ROOT, "results")
os.makedirs(OUT_FIG, exist_ok=True)

EA = sim.EA_DEFAULT
T_REF = sim.T_REF_DEFAULT
T90_REF = sim.T90_REF_DEFAULT
A = sim.pre_exponential(EA, T_REF, T90_REF)

TEMPS_C = [60, 70, 80, 90, 100]


def main() -> int:
    print(f"model: Ea={EA} kJ/mol, t90(100C)={T90_REF/3600:.1f} h  ->  A={A:.3e} /s")

    # --- table -------------------------------------------------------------
    rows = []
    T_list, k_list = [], []
    for tc in TEMPS_C:
        T = tc + 273.15
        k = sim.rate_constant(A, EA, T)
        t90 = sim.time_to_conversion(k, 0.9)
        rows.append({
            "T_C": tc,
            "k_per_s": f"{k:.3e}",
            "t90_s": f"{t90:.0f}",
            "t90_h": f"{t90/3600:.2f}",
            "t90_days": f"{t90/86400:.2f}",
            "alpha_1h": f"{1-np.exp(-k*3600):.4f}",
        })
        T_list.append(T)
        k_list.append(k)

    with open(os.path.join(OUT_CSV, "kinetics_simulation.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote results/kinetics_simulation.csv")
    for r in rows:
        print("  ", r)

    # --- validation: RK4 vs analytic --------------------------------------
    T = 80 + 273.15
    k = sim.rate_constant(A, EA, T)
    ts, alphas = sim.integrate(A, EA, T, t_end=sim.time_to_conversion(k, 0.9) * 1.3,
                               dt=max(1.0, sim.time_to_conversion(k, 0.9) / 2000))
    analytic = 1 - np.exp(-k * np.array(ts))
    max_err = float(np.max(np.abs(np.array(alphas) - analytic)))
    print(f"RK4 vs analytic max abs error (80 C): {max_err:.2e}")

    # --- recover Ea --------------------------------------------------------
    Ea_fit, A_fit = sim.fit_arrhenius(T_list, k_list)
    print(f"Arrhenius fit recovers Ea = {Ea_fit:.1f} kJ/mol (input {EA})")

    # --- kinetic bottleneck ratio -----------------------------------------
    t90_60 = sim.time_to_conversion(sim.rate_constant(A, EA, 333.15))
    t90_90 = sim.time_to_conversion(sim.rate_constant(A, EA, 363.15))
    print(f"t90(60C)/t90(90C) = {t90_60/t90_90:.0f}x")

    # --- cycle energy ------------------------------------------------------
    for rid in ("gypsum", "glauber"):
        rxn = next(r for r in REACTIONS if r["id"] == rid)
        ce = sim.cycle_energy(rxn, alpha=0.9)
        print(f"cycle[{rid}]: stored {ce['stored_kJ_kg']:.0f} kJ/kg at 90% conversion")

    # --- figure ------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.8))
    tt = np.logspace(1, 7, 300)
    for tc in TEMPS_C:
        k = sim.rate_constant(A, EA, tc + 273.15)
        ax1.plot(tt, 1 - np.exp(-k * tt), lw=1.8, label=f"{tc} C")
    ax1.axhline(0.9, color="k", ls="--", lw=0.8)
    ax1.set_xscale("log")
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Conversion alpha")
    ax1.set_title("(a) Dehydration conversion vs time")
    ax1.legend(fontsize=8, loc="upper left")
    ax1.set_ylim(0, 1.02)

    invT = 1.0 / np.array(T_list)
    ax2.plot(invT * 1000, np.log(k_list), "o-", color="#264653")
    ax2.set_xlabel("1000/T  (1/K)")
    ax2.set_ylabel("ln k")
    ax2.set_title(f"(b) Arrhenius fit: Ea = {Ea_fit:.0f} kJ/mol")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_FIG, "fig5_kinetics.png"), dpi=200)
    print("wrote fig5_kinetics.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
