"""Equilibrium water-vapour pressure of a dehydration reaction.

For a reaction that releases ``n`` moles of gas,

    M.nH2O(s)  <=>  M(s) + n H2O(g)

the equilibrium constant is K = (p/p0)^n and dG(T) = -R T ln K, so

    ln(p_eq/p0) = dS/(n R) - dH/(n R T)

This is the van 't Hoff form. It matters because a dehydration does *not* need
to reach T* (where p_eq = 1 bar) to proceed: it only needs the ambient vapour
pressure to fall below p_eq(T). In ordinary air (p ~ 0.02 bar) that happens at a
much lower temperature than T*, which shifts the practical limit from
thermodynamics to kinetics.
"""

from __future__ import annotations

import math

from .data import SPECIES
from .thermo import delta_h, delta_s, R

P0 = 1.0  # bar, reference pressure

# Typical ambient water-vapour partial pressures (bar).
AMBIENT = {
    "dry_air_20RH": 0.005,
    "typical_50RH": 0.016,
    "humid_80RH": 0.025,
}


def gas_moles(rxn: dict) -> float:
    """Moles of gaseous product per reaction as written."""
    n = 0.0
    for sp, coeff in rxn["products"].items():
        if "(g)" in sp:
            n += coeff
    return n


def ln_peq_ratio(rxn: dict, T: float) -> float:
    n = gas_moles(rxn)
    if n == 0:
        raise ValueError("reaction releases no gas; p_eq is undefined")
    dH = delta_h(rxn) * 1000.0
    dS = delta_s(rxn)
    return dS / (n * R) - dH / (n * R * T)


def peq(rxn: dict, T: float) -> float:
    """Equilibrium gas partial pressure at temperature T (K), in bar."""
    return P0 * math.exp(ln_peq_ratio(rxn, T))


def effective_threshold(rxn: dict, p_ambient: float) -> float | None:
    """Temperature (K) above which p_eq(T) > p_ambient, i.e. dehydration proceeds.

    Solves ln(p_ambient/p0) = dS/(nR) - dH/(nRT) for T.
    """
    n = gas_moles(rxn)
    if n == 0 or p_ambient <= 0:
        return None
    dH = delta_h(rxn) * 1000.0
    dS = delta_s(rxn)
    target = math.log(p_ambient / P0)
    denom = dS / (n * R) - target
    if abs(denom) < 1e-12:
        return None
    return dH / (n * R) / denom
