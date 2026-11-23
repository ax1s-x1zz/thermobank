"""Calorimetry analysis for heat-of-solution experiments.

A student-level but rigorous calorimetry routine: measure the temperature rise
of a known mass of water when a known mass of salt dissolves, then convert to a
molar enthalpy of solution and compare with the literature value.

The heat balance used is
    q_water = m_water * c_water * dT
    dH_sol  = -q_water / n_solute          (kJ/mol)
"""

from __future__ import annotations

C_WATER = 4.184  # J/(g*K)


def heat_released(m_water_g: float, dT: float, c: float = C_WATER) -> float:
    """Heat absorbed by the water, J."""
    return m_water_g * c * dT


def enthalpy_of_solution(m_salt_g: float, M_salt: float, m_water_g: float,
                         dT: float, c: float = C_WATER) -> float:
    """Molar enthalpy of solution, kJ/mol (negative = exothermic)."""
    n = m_salt_g / M_salt
    q = heat_released(m_water_g, dT, c)  # J absorbed by water
    return -q / n / 1000.0


def percent_error(measured: float, literature: float) -> float:
    return (measured - literature) / literature * 100.0


def heat_capacity_of_calorimeter(m_water_g: float, m_salt_g: float, M_salt: float,
                                 dT: float, dH_lit_kJ: float) -> float:
    """Estimate the calorimeter heat capacity (J/K) from a known reaction.

    Solving  m*c*dT + C_cal*dT = -n*dH_lit  for C_cal lets the calorimeter
    constant be calibrated from a salt whose dH_sol is well known (e.g. KCl).
    """
    n = m_salt_g / M_salt
    q_rxn = -n * dH_lit_kJ * 1000.0          # J released by reaction (<0 exothermic)
    q_water = m_water_g * C_WATER * dT       # J absorbed by water
    return (q_rxn - q_water) / dT
