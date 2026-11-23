"""Simple kinetics helpers.

Thermodynamic spontaneity (dG < 0) does not guarantee a usable reaction rate.
This module quantifies that gap with an Arrhenius model so the screening can
flag reactions that are thermodynamically allowed but kinetically slow.
"""

from __future__ import annotations

from .thermo import R


def rate_constant(A: float, Ea_kJ: float, T: float) -> float:
    """Arrhenius rate constant k = A * exp(-Ea / (R T)).

    ``Ea_kJ`` in kJ/mol, ``T`` in K.
    """
    return A * __import__("math").exp(-Ea_kJ * 1000.0 / (R * T))


def relative_rate(Ea_kJ: float, T: float, T_ref: float) -> float:
    """k(T) / k(T_ref) with the pre-exponential factor cancelled out.

    A convenient dimensionless measure of how much faster the reaction runs at
    ``T`` than at a reference temperature.
    """
    import math
    return math.exp(-Ea_kJ * 1000.0 / R * (1.0 / T - 1.0 / T_ref))


def time_to_fraction(Ea_kJ: float, T: float, T_ref: float, tau_ref: float) -> float:
    """Characteristic time to reach a fixed conversion, scaled from a reference.

    ``tau_ref`` is the measured (or assumed) time to that conversion at
    ``T_ref``; the result is tau_ref * k(T_ref)/k(T).
    """
    return tau_ref / relative_rate(Ea_kJ, T, T_ref)
