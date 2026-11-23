"""Uncertainty propagation for the threshold temperature T*.

T* = dH / dS inherits the uncertainty of the formation data. Gas and simple
solids are good to ~1 kJ/mol and ~2 J/(mol*K); salt-hydrate values are much
looser. This module runs a Monte Carlo over those uncertainties and reports a
central value with a 5-95% interval, so the screening table carries error bars
instead of false precision.
"""

from __future__ import annotations

import random

from .data import SPECIES, REACTIONS
from .thermo import delta_h, delta_s

# (sigma_dHf kJ/mol, sigma_S J/mol/K) per quality class
SIGMA = {
    "ref": (1.0, 2.0),
    "approx": (6.0, 12.0),
}


def _perturbed(rxn: dict, rng: random.Random) -> tuple[float, float]:
    dH = 0.0
    dS = 0.0
    for side, sign in (("products", +1), ("reactants", -1)):
        for sp, coeff in rxn[side].items():
            sH, sS = SIGMA[SPECIES[sp]["quality"]]
            dH += sign * coeff * rng.gauss(SPECIES[sp]["dHf"], sH)
            dS += sign * coeff * rng.gauss(SPECIES[sp]["S"], sS)
    return dH, dS


def tstar_distribution(rxn: dict, n: int = 4000, seed: int = 20261123) -> dict:
    rng = random.Random(seed)
    vals = []
    for _ in range(n):
        dH, dS = _perturbed(rxn, rng)
        if abs(dS) < 1e-9:
            continue
        vals.append(dH * 1000.0 / dS)
    vals.sort()
    k = len(vals)
    return {
        "id": rxn["id"],
        "median": vals[k // 2],
        "p05": vals[int(0.05 * k)],
        "p95": vals[int(0.95 * k)],
        "n": k,
    }


def all_distributions(n: int = 4000) -> list[dict]:
    return [tstar_distribution(r, n) for r in REACTIONS]
