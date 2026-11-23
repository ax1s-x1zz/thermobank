"""Standard thermochemical data used by thermobank.

All values are at 298.15 K.

* ``dHf``  standard molar enthalpy of formation, kJ/mol
* ``S``    standard molar entropy, J/(mol*K)
* ``M``    molar mass, g/mol

Primary sources
---------------
* CRC Handbook of Chemistry and Physics, 97th ed. (2016-2017).
* NIST Chemistry WebBook, SRD 69.
* Chase, M. W. (1998), NIST-JANAF Thermochemical Tables, 4th ed.

Hydrate entries are marked ``approx``. Their enthalpies and entropies carry
larger uncertainty than the gas/solid reference species, so any result that
depends on them is reported with that caveat (see ``docs/methods.md``).
"""

from __future__ import annotations

# name: {dHf (kJ/mol), S (J/mol/K), M (g/mol), quality, source}
SPECIES: dict[str, dict] = {
    # --- reference gases / water -------------------------------------------
    "H2O(g)":   {"dHf": -241.82, "S": 188.83, "M": 18.015, "quality": "ref", "source": "NIST"},
    "H2O(l)":   {"dHf": -285.83, "S":  69.95, "M": 18.015, "quality": "ref", "source": "NIST"},
    "CO2(g)":   {"dHf": -393.51, "S": 213.74, "M": 44.010, "quality": "ref", "source": "NIST"},
    "H2(g)":    {"dHf":    0.00, "S": 130.68, "M":  2.016, "quality": "ref", "source": "NIST"},

    # --- calcium system -----------------------------------------------------
    "CaO(s)":       {"dHf":  -635.09, "S":  38.10, "M":  56.077, "quality": "ref",    "source": "CRC"},
    "Ca(OH)2(s)":   {"dHf":  -986.09, "S":  83.40, "M":  74.093, "quality": "ref",    "source": "CRC"},
    "CaCO3(s)":     {"dHf": -1206.92, "S":  92.90, "M": 100.086, "quality": "ref",    "source": "CRC"},
    "CaSO4(s)":     {"dHf": -1434.11, "S": 106.70, "M": 136.140, "quality": "ref",    "source": "CRC"},
    "CaSO4.2H2O(s)":{"dHf": -2022.63, "S": 194.10, "M": 172.170, "quality": "approx", "source": "CRC"},
    "CaCl2(s)":     {"dHf":  -795.80, "S": 104.60, "M": 110.980, "quality": "ref",    "source": "CRC"},
    "CaCl2.6H2O(s)":{"dHf": -2607.90, "S": 300.00, "M": 219.080, "quality": "approx", "source": "CRC"},

    # --- magnesium system ---------------------------------------------------
    "MgO(s)":       {"dHf":  -601.70, "S":  26.90, "M":  40.304, "quality": "ref",    "source": "CRC"},
    "Mg(OH)2(s)":   {"dHf":  -924.54, "S":  63.20, "M":  58.320, "quality": "ref",    "source": "CRC"},
    "MgH2(s)":      {"dHf":   -75.30, "S":  31.10, "M":  26.321, "quality": "ref",    "source": "CRC"},
    "Mg(s)":        {"dHf":     0.00, "S":  32.70, "M":  24.305, "quality": "ref",    "source": "CRC"},

    # --- copper / sodium / barium ------------------------------------------
    "CuSO4(s)":       {"dHf":  -771.36, "S": 109.20, "M": 159.610, "quality": "ref",    "source": "CRC"},
    "CuSO4.5H2O(s)":  {"dHf": -2279.70, "S": 300.40, "M": 249.680, "quality": "approx", "source": "CRC"},
    "Na2SO4(s)":      {"dHf": -1387.10, "S": 149.60, "M": 142.040, "quality": "ref",    "source": "CRC"},
    "Na2SO4.10H2O(s)":{"dHf": -4327.00, "S": 590.00, "M": 322.190, "quality": "approx", "source": "CRC"},
    "BaO(s)":         {"dHf":  -553.50, "S":  70.40, "M": 153.330, "quality": "ref",    "source": "CRC"},
    "Ba(OH)2(s)":     {"dHf":  -944.70, "S": 107.10, "M": 171.340, "quality": "ref",    "source": "CRC"},
}


# Candidate thermochemical energy-storage reactions.
# Dehydration / decomposition reactions are written forward (heat in, endothermic);
# discharging is the reverse (heat out, exothermic).
REACTIONS: list[dict] = [
    {
        "id": "gypsum",
        "label": "CaSO4.2H2O -> CaSO4 + 2H2O(g)",
        "class": "salt hydrate",
        "reactants": {"CaSO4.2H2O(s)": 1},
        "products": {"CaSO4(s)": 1, "H2O(g)": 2},
        "note": "gypsum/bassanite; classic low-temperature dehydration",
    },
    {
        "id": "glauber",
        "label": "Na2SO4.10H2O -> Na2SO4 + 10H2O(g)",
        "class": "salt hydrate",
        "reactants": {"Na2SO4.10H2O(s)": 1},
        "products": {"Na2SO4(s)": 1, "H2O(g)": 10},
        "note": "Glauber's salt; high water uptake, large entropy term",
    },
    {
        "id": "cacl2_hydrate",
        "label": "CaCl2.6H2O -> CaCl2 + 6H2O(g)",
        "class": "salt hydrate",
        "reactants": {"CaCl2.6H2O(s)": 1},
        "products": {"CaCl2(s)": 1, "H2O(g)": 6},
        "note": "high energy density, deliquescent",
    },
    {
        "id": "cuso4_hydrate",
        "label": "CuSO4.5H2O -> CuSO4 + 5H2O(g)",
        "class": "salt hydrate",
        "reactants": {"CuSO4.5H2O(s)": 1},
        "products": {"CuSO4(s)": 1, "H2O(g)": 5},
        "note": "textbook hydrate dehydration",
    },
    {
        "id": "mgoh2",
        "label": "Mg(OH)2 -> MgO + H2O(g)",
        "class": "hydroxide",
        "reactants": {"Mg(OH)2(s)": 1},
        "products": {"MgO(s)": 1, "H2O(g)": 1},
        "note": "medium-temperature hydroxide",
    },
    {
        "id": "caoh2",
        "label": "Ca(OH)2 -> CaO + H2O(g)",
        "class": "hydroxide",
        "reactants": {"Ca(OH)2(s)": 1},
        "products": {"CaO(s)": 1, "H2O(g)": 1},
        "note": "CaO/Ca(OH)2 loop; widely studied",
    },
    {
        "id": "baoh2",
        "label": "Ba(OH)2 -> BaO + H2O(g)",
        "class": "hydroxide",
        "reactants": {"Ba(OH)2(s)": 1},
        "products": {"BaO(s)": 1, "H2O(g)": 1},
        "note": "high-temperature hydroxide",
    },
    {
        "id": "mgh2",
        "label": "MgH2 -> Mg + H2(g)",
        "class": "hydride",
        "reactants": {"MgH2(s)": 1},
        "products": {"Mg(s)": 1, "H2(g)": 1},
        "note": "chemical hydrogen storage; no water required",
    },
    {
        "id": "caco3",
        "label": "CaCO3 -> CaO + CO2(g)",
        "class": "carbonate",
        "reactants": {"CaCO3(s)": 1},
        "products": {"CaO(s)": 1, "CO2(g)": 1},
        "note": "calcium looping; high temperature",
    },
]
