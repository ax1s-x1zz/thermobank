"""Command-line interface for thermobank."""

from __future__ import annotations

import argparse
import sys

from . import thermo
from .data import REACTIONS
from .screen import run_screen, to_dataframe, WASTE_HEAT
from . import calorimetry as calo
from . import kinetics


def _cmd_screen(args: argparse.Namespace) -> int:
    df = to_dataframe(run_screen())
    if args.csv:
        df.to_csv(args.csv, index=False)
        print(f"wrote {args.csv} ({len(df)} reactions)")
    with __import__("pandas").option_context("display.width", 200,
                                             "display.max_columns", 20):
        print(df.to_string(index=False))
    return 0


def _cmd_reaction(args: argparse.Namespace) -> int:
    rxn = next((r for r in REACTIONS if r["id"] == args.id), None)
    if rxn is None:
        print(f"unknown reaction id: {args.id}", file=sys.stderr)
        return 2
    T = args.temp + 273.15
    print(f"{rxn['label']}")
    print(f"  dH(298)          = {thermo.delta_h(rxn):+.1f} kJ/mol")
    print(f"  dS(298)          = {thermo.delta_s(rxn):+.1f} J/(mol*K)")
    print(f"  dG({args.temp:.0f} C)      = {thermo.delta_g(rxn, T)/1000:+.1f} kJ/mol "
          f"({'spontaneous' if thermo.is_spontaneous(rxn, T) else 'not spontaneous'})")
    ts = thermo.threshold_temperature(rxn)
    if ts:
        print(f"  T*               = {ts:.0f} K ({ts-273.15:.0f} C)")
    print(f"  energy density   = {thermo.energy_density(rxn):.0f} kJ/kg")
    return 0


def _cmd_calo(args: argparse.Namespace) -> int:
    dH = calo.enthalpy_of_solution(args.mass, args.M, args.water, args.dT)
    print(f"dH_solution = {dH:+.1f} kJ/mol  ({'exothermic' if dH < 0 else 'endothermic'})")
    if args.lit is not None:
        err = calo.percent_error(dH, args.lit)
        print(f"literature  = {args.lit:+.1f} kJ/mol  ->  error {err:+.1f}%")
    return 0


def _cmd_kinetics(args: argparse.Namespace) -> int:
    r = kinetics.relative_rate(args.Ea, args.temp + 273.15, args.tref + 273.15)
    print(f"relative rate k({args.temp:.0f}C)/k({args.tref:.0f}C) = {r:.3g}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="thermobank",
        description="Thermochemical heat-storage screening for low-grade waste heat.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("screen", help="screen all candidate reactions")
    s.add_argument("--csv", help="write results to this CSV path")
    s.set_defaults(func=_cmd_screen)

    r = sub.add_parser("reaction", help="inspect one reaction")
    r.add_argument("id")
    r.add_argument("--temp", type=float, default=80.0, help="temperature in C")
    r.set_defaults(func=_cmd_reaction)

    c = sub.add_parser("calo", help="heat of solution from calorimetry data")
    c.add_argument("--mass", type=float, required=True, help="solute mass, g")
    c.add_argument("--M", type=float, required=True, help="solute molar mass, g/mol")
    c.add_argument("--water", type=float, required=True, help="water mass, g")
    c.add_argument("--dT", type=float, required=True, help="temperature change, K")
    c.add_argument("--lit", type=float, default=None, help="literature dH_sol, kJ/mol")
    c.set_defaults(func=_cmd_calo)

    k = sub.add_parser("kinetics", help="relative Arrhenius rate")
    k.add_argument("--Ea", type=float, required=True, help="activation energy, kJ/mol")
    k.add_argument("--temp", type=float, required=True, help="temperature, C")
    k.add_argument("--tref", type=float, default=25.0, help="reference temp, C")
    k.set_defaults(func=_cmd_kinetics)

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
