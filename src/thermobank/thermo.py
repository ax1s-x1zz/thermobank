"""Core thermochemistry: reaction enthalpy, entropy, Gibbs energy, threshold T."""

from __future__ import annotations

from .data import SPECIES

R = 8.314462618  # J/(mol*K)


def _sum(rxn: dict, side: str, key: str) -> float:
    total = 0.0
    for sp, coeff in rxn[side].items():
        total += coeff * SPECIES[sp][key]
    return total


def delta_h(rxn: dict) -> float:
    """Standard reaction enthalpy at 298 K, kJ/mol (per reaction as written)."""
    return _sum(rxn, "products", "dHf") - _sum(rxn, "reactants", "dHf")


def delta_s(rxn: dict) -> float:
    """Standard reaction entropy at 298 K, J/(mol*K)."""
    return _sum(rxn, "products", "S") - _sum(rxn, "reactants", "S")


def delta_g(rxn: dict, T: float) -> float:
    """Gibbs energy of reaction at temperature ``T`` (K), J/mol.

    Uses the constant-dH/dS (van 't Hoff style) approximation:
        dG(T) ~= dH(298) - T * dS(298)
    This is the approximation used in the Korean high-school chemistry
    curriculum for judging spontaneity, and is accurate enough to locate the
    threshold temperature to within tens of degrees.
    """
    return delta_h(rxn) * 1000.0 - T * delta_s(rxn)


def threshold_temperature(rxn: dict) -> float | None:
    """Temperature (K) at which dG = 0, i.e. T* = dH/dS.

    For an endothermic dehydration (dH > 0, dS > 0) the reaction is
    spontaneous for T > T*. Returns ``None`` if dS is (near) zero.
    """
    ds = delta_s(rxn)
    if abs(ds) < 1e-9:
        return None
    return delta_h(rxn) * 1000.0 / ds


def reactant_mass(rxn: dict) -> float:
    """Molar mass of the storage material (reactants), g/mol."""
    return _sum(rxn, "reactants", "M")


def energy_density(rxn: dict) -> float:
    """Gravimetric energy density of the storage material, kJ/kg.

    The stored energy is |dH| per mole of reactant consumed.
    """
    m_kg = reactant_mass(rxn) / 1000.0
    return abs(delta_h(rxn)) / m_kg


def is_spontaneous(rxn: dict, T: float) -> bool:
    return delta_g(rxn, T) < 0.0
