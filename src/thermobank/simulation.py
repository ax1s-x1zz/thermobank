"""Computational experiment: kinetics of a thermochemical storage cycle.

This is a *computational* experiment, not a wet-lab measurement. It integrates a
lumped kinetic model of the dehydration reaction

    d(alpha)/dt = k(T) * (1 - alpha),        k(T) = A * exp(-Ea / (R T))

to obtain the conversion alpha(t) at a given heat-source temperature, and a
simple energy balance for a charge/discharge cycle. The point is to quantify the
kinetic bottleneck: thermodynamics says the dehydration is allowed at 60-90 C,
but the model shows how the *rate* collapses at low temperature.

The activation energy and the reference time are model inputs taken from the
literature range for salt-hydrate dehydration; they are stated explicitly so the
result can be recomputed with other values.
"""

from __future__ import annotations

import math

from .thermo import R, delta_h
from .equilibrium import gas_moles

# Literature-range kinetic inputs for gypsum-type dehydration (assumptions).
EA_DEFAULT = 150.0        # kJ/mol
T_REF_DEFAULT = 373.15    # K   (100 C)
T90_REF_DEFAULT = 3600.0  # s   (time to 90% conversion at T_ref)


def pre_exponential(Ea_kJ: float, T_ref: float, t90_ref: float) -> float:
    """A calibrated so that the model reaches 90% conversion in t90_ref at T_ref.

    First-order: alpha = 1 - exp(-k t)  =>  k = -ln(0.1)/t90 = ln(10)/t90.
    """
    k_ref = math.log(10.0) / t90_ref
    return k_ref * math.exp(Ea_kJ * 1000.0 / (R * T_ref))


def rate_constant(A: float, Ea_kJ: float, T: float) -> float:
    return A * math.exp(-Ea_kJ * 1000.0 / (R * T))


def time_to_conversion(k: float, alpha: float = 0.9) -> float:
    """First-order time to reach ``alpha`` (analytic)."""
    return -math.log(1.0 - alpha) / k


def integrate(A: float, Ea_kJ: float, T: float, t_end: float, dt: float = 1.0):
    """RK4 integration of d(alpha)/dt = k(1-alpha). Returns (t[], alpha[])."""
    k = rate_constant(A, Ea_kJ, T)
    f = lambda a: k * (1.0 - a)  # noqa: E731
    t, a = 0.0, 0.0
    ts, alphas = [t], [a]
    while t < t_end and a < 0.999999:
        k1 = f(a)
        k2 = f(a + 0.5 * dt * k1)
        k3 = f(a + 0.5 * dt * k2)
        k4 = f(a + dt * k3)
        a += dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt
        ts.append(t)
        alphas.append(a)
    return ts, alphas


def fit_arrhenius(T_list, k_list):
    """Recover (Ea [kJ/mol], A) from ln k vs 1/T by least squares."""
    import numpy as np
    x = np.array([1.0 / t for t in T_list])
    y = np.array([math.log(k) for k in k_list])
    slope, intercept = np.polyfit(x, y, 1)
    return -slope * R / 1000.0, math.exp(intercept)


def cycle_energy(rxn: dict, alpha: float = 0.9):
    """Energy stored/released per mole and per kg for a given conversion."""
    from .thermo import reactant_mass
    dH = delta_h(rxn)              # kJ/mol, endothermic (stored)
    m_kg = reactant_mass(rxn) / 1000.0
    return {
        "stored_kJ_mol": alpha * dH,
        "stored_kJ_kg": alpha * dH / m_kg,
        "gas_mol": gas_moles(rxn),
    }
