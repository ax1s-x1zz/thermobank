"""Screening of candidate reactions against waste-heat temperature windows."""

from __future__ import annotations

from dataclasses import dataclass

from .data import REACTIONS, SPECIES
from . import thermo

# Temperature windows of interest (K).
WASTE_HEAT = {
    "data_center": (333.15, 363.15),      # 60-90 C, liquid-cooled racks
    "industrial_low": (423.15, 473.15),   # 150-200 C, low-grade process heat
    "industrial_high": (573.15, 773.15),  # 300-500 C
}


@dataclass
class ScreenRow:
    id: str
    label: str
    cls: str
    dH: float              # kJ/mol
    dS: float              # J/mol/K
    T_star: float          # K
    energy_density: float  # kJ/kg
    quality: str           # "ref" or "approx"

    def dG_at(self, T: float) -> float:
        rxn = next(r for r in REACTIONS if r["id"] == self.id)
        return thermo.delta_g(rxn, T)

    def fits(self, window: str) -> bool:
        lo, hi = WASTE_HEAT[window]
        return self.T_star is not None and lo <= self.T_star <= hi


def _quality(rxn: dict) -> str:
    for side in ("reactants", "products"):
        for sp in rxn[side]:
            if SPECIES[sp]["quality"] == "approx":
                return "approx"
    return "ref"


def run_screen() -> list[ScreenRow]:
    rows: list[ScreenRow] = []
    for rxn in REACTIONS:
        rows.append(ScreenRow(
            id=rxn["id"],
            label=rxn["label"],
            cls=rxn["class"],
            dH=thermo.delta_h(rxn),
            dS=thermo.delta_s(rxn),
            T_star=thermo.threshold_temperature(rxn),
            energy_density=thermo.energy_density(rxn),
            quality=_quality(rxn),
        ))
    return rows


def to_dataframe(rows: list[ScreenRow]):
    import pandas as pd
    data = []
    for r in rows:
        data.append({
            "id": r.id,
            "reaction": r.label,
            "class": r.cls,
            "dH_kJ_mol": round(r.dH, 1),
            "dS_J_molK": round(r.dS, 1),
            "T_star_K": round(r.T_star, 0) if r.T_star else None,
            "T_star_C": round(r.T_star - 273.15, 0) if r.T_star else None,
            "energy_density_kJ_kg": round(r.energy_density, 0),
            "data_center": r.fits("data_center"),
            "industrial_low": r.fits("industrial_low"),
            "industrial_high": r.fits("industrial_high"),
            "data_quality": r.quality,
        })
    df = pd.DataFrame(data)
    return df.sort_values("T_star_K").reset_index(drop=True)
