"""Unit tests for thermobank."""

from __future__ import annotations

import math

from thermobank import thermo, kinetics, calorimetry, screen
from thermobank.data import REACTIONS


def test_caco3_enthalpy():
    rxn = next(r for r in REACTIONS if r["id"] == "caco3")
    # dH ~ +178 kJ/mol, dS ~ +160 J/mol/K, T* ~ 1110 K
    assert abs(thermo.delta_h(rxn) - 178.3) < 1.0
    assert abs(thermo.delta_s(rxn) - 158.9) < 2.0
    assert abs(thermo.threshold_temperature(rxn) - 1122) < 30


def test_gypsum_low_temperature():
    rxn = next(r for r in REACTIONS if r["id"] == "gypsum")
    # gypsum dehydration is thermodynamically allowed below 100 C
    assert thermo.threshold_temperature(rxn) < 373.15


def test_delta_g_sign_changes_at_threshold():
    rxn = next(r for r in REACTIONS if r["id"] == "caoh2")
    ts = thermo.threshold_temperature(rxn)
    assert thermo.delta_g(rxn, ts - 20) > 0
    assert thermo.delta_g(rxn, ts + 20) < 0


def test_energy_density_positive():
    for r in REACTIONS:
        assert thermo.energy_density(r) > 0


def test_arrhenius_ratio_one_at_reference():
    assert abs(kinetics.relative_rate(100.0, 350.0, 350.0) - 1.0) < 1e-12


def test_heat_of_solution_exothermic_sign():
    # dissolving 5 g CaCl2 (M=110.98) in 50 g water raising it by 8 K
    dH = calorimetry.enthalpy_of_solution(5.0, 110.98, 50.0, 8.0)
    assert dH < 0  # exothermic


def test_screen_runs():
    rows = screen.run_screen()
    assert len(rows) == len(REACTIONS)
    assert all(r.T_star is not None for r in rows)
