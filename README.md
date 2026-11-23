# thermobank

**Thermochemical heat-storage screening for low-grade waste heat.**

thermobank asks a single question: *given a waste-heat stream at some
temperature, which reversible chemical reactions can it drive, and how much
energy can they store?* It screens candidate dehydration / decomposition
reactions using the same thermodynamics taught in high-school chemistry
(`dG = dH - T dS`) and reports the threshold temperature `T* = dH / dS` at which
each reaction becomes spontaneous.

The motivation is concrete: data centres reject large amounts of low-grade heat
(60-90 C from liquid-cooled racks), and industry rejects more at 150-200 C.
Storing that heat chemically instead of throwing it away is an active research
area (salt-hydrate and hydroxide loops, zeolite sorption batteries).

## Quick start

```bash
pip install -e .
thermobank screen --csv results/screening.csv
thermobank reaction gypsum --temp 80
```

```
           id                          reaction        class  dH_kJ_mol  dS_J_molK  T_star_C  energy_density_kJ_kg  data_center
      glauber Na2SO4.10H2O -> Na2SO4 + 10H2O(g) salt hydrate      521.7     1447.9      87.0                1619.0         True
       gypsum     CaSO4.2H2O -> CaSO4 + 2H2O(g) salt hydrate      104.9      290.3      88.0                 609.0         True
 cacl2_hydrate     CaCl2.6H2O -> CaCl2 + 6H2O(g) salt hydrate      361.2      937.6     112.0                1649.0        False
        mgoh2           Mg(OH)2 -> MgO + H2O(g)    hydroxide       81.0      152.5     258.0                1389.0        False
        caco3             CaCO3 -> CaO + CO2(g)    carbonate      178.3      158.9     849.0                1782.0        False
```

## What it computes

| module | quantity |
|---|---|
| `thermo` | reaction `dH`, `dS`, `dG(T)`, threshold temperature `T*`, energy density |
| `screen` | ranks every candidate against data-centre / industrial heat windows |
| `kinetics` | Arrhenius relative rate — spontaneity is not the same as speed |
| `calorimetry` | molar heat of solution from a simple calorimetry measurement |

## Method (one line)

`dG(T) = dH(298) - T dS(298)`, with `dH` and `dS` from standard formation
enthalpies and entropies. A dehydration with `dH > 0`, `dS > 0` is spontaneous
above `T* = dH/dS`. The approximation is the one used in the curriculum and is
accurate enough to place `T*` within tens of degrees. See `docs/methods.md`.

## Data provenance

Standard formation data come from the CRC Handbook and the NIST Chemistry
WebBook. Salt-hydrate entries are marked `approx` and carry larger uncertainty;
results that depend on them are flagged. See `src/thermobank/data.py`.

## Limitations

* `dH` and `dS` are treated as temperature-independent.
* Hydrate thermodynamic data are less well established than gas/solid data.
* Thermodynamic feasibility does not include kinetics, mass transport, or
  cycle stability — those need experiment.

## License

Apache-2.0. See `LICENSE`.

## Author

Jangseowoo, `ax1s-x1zz`, [ax1s@x1zz.com](mailto:ax1s@x1zz.com).
