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


def test_equilibrium_threshold_below_boiling():
    from thermobank import equilibrium as eq
    rxn = next(r for r in REACTIONS if r["id"] == "gypsum")
    t_eff = eq.effective_threshold(rxn, eq.AMBIENT["typical_50RH"])
    # in ordinary air the effective threshold is far below the 1-bar T*
    assert t_eff < 300.0  # < 27 C


def test_simulation_arrhenius_roundtrip():
    from thermobank import simulation as sim
    A = sim.pre_exponential(sim.EA_DEFAULT, sim.T_REF_DEFAULT, sim.T90_REF_DEFAULT)
    Ts = [333.15, 353.15, 373.15]
    ks = [sim.rate_constant(A, sim.EA_DEFAULT, T) for T in Ts]
    Ea_fit, _ = sim.fit_arrhenius(Ts, ks)
    assert abs(Ea_fit - sim.EA_DEFAULT) < 0.5


def test_simulation_rk4_matches_analytic():
    from thermobank import simulation as sim
    A = sim.pre_exponential(sim.EA_DEFAULT, sim.T_REF_DEFAULT, sim.T90_REF_DEFAULT)
    T = 353.15
    k = sim.rate_constant(A, sim.EA_DEFAULT, T)
    ts, alphas = sim.integrate(A, sim.EA_DEFAULT, T,
                               t_end=sim.time_to_conversion(k, 0.9) * 1.2,
                               dt=sim.time_to_conversion(k, 0.9) / 1000)
    err = max(abs(a - (1 - math.exp(-k * t))) for t, a in zip(ts, alphas))
    assert err < 1e-6
