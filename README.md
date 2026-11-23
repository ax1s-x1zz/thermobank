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
| `equilibrium` | equilibrium vapour pressure `p_eq(T)` and the *effective* threshold in ambient air |
| `kinetics` | Arrhenius relative rate — spontaneity is not the same as speed |
| `uncertainty` | Monte Carlo propagation of formation-data errors into `T*` |
| `calorimetry` | molar heat of solution from a simple calorimetry measurement |

## Key findings

1. **The 1-bar threshold is misleading.** `T*` (where `p_eq = 1 bar`) is 87-88 C
   for gypsum and Glauber's salt, which looks like a hard floor for low-grade
   heat. But a dehydration only needs `p_eq(T) > ambient humidity`. In ordinary
   50% RH air (~0.016 bar) the effective threshold drops to **~19 C** for both —
   far below data-centre waste heat.
2. **So the real bottleneck is kinetics, not thermodynamics.** At 60-90 C the
   dehydration is thermodynamically favoured; whether it happens fast enough is
   set by the activation energy. This reframes the design question.
3. **Uncertainty is large for hydrates.** Propagating formation-data errors
   gives gypsum `T* = 87 C` with a 5-95% interval of **[45, 134] C**, versus a
   tight [72, 103] C for Glauber's salt. Screening tables need error bars.


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
