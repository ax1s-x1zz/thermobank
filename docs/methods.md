# Methods

## 1. Thermodynamic model

For a reaction written forward (heat in, endothermic),

```
dH_rxn(298) = sum dHf(products) - sum dHf(reactants)      [kJ/mol]
dS_rxn(298) = sum S(products)   - sum S(reactants)        [J/(mol*K)]
dG_rxn(T)   = dH_rxn(298) - T * dS_rxn(298)               [J/mol]
```

A reaction is spontaneous when `dG < 0`. For a dehydration that releases a gas
(`dH > 0`, `dS > 0`) the sign of `dG` flips at

```
T* = dH / dS        (Kelvin)
```

Below `T*` the reaction is not spontaneous; above `T*` it is. This is exactly
the criterion used in the Korean high-school chemistry curriculum for judging
spontaneity, so the calculation stays inside the subject.

## 2. Assumptions

1. **Constant dH and dS.** Both are taken at 298 K and assumed independent of
   temperature. The neglected heat-capacity term is typically a few percent of
   `dH` over a 300 K range; `T*` shifts by tens of degrees, not hundreds.
2. **Standard-state, pure condensed phases.** Activities of solids are 1;
   released water is treated as an ideal gas at 1 bar.
3. **Single-step reaction.** Intermediate hydrates (e.g. bassanite,
   `CaSO4.0.5H2O`) are not modelled; only the full dehydration is screened.

## 3. Uncertainty

* Gas and simple-solid data (NIST, CRC) are good to ~1 kJ/mol, giving `T*` to
  roughly +-10 K.
* Salt-hydrate data (`quality: approx`) can be off by several kJ/mol, so `T*`
  for those rows should be read as +-20-40 K. Rows are flagged in the output.

## 4. Kinetics

Thermodynamic feasibility is necessary but not sufficient. The Arrhenius form

```
k(T) = A * exp(-Ea / (R T))
```

is used only in *relative* terms, `k(T)/k(T_ref)`, so the unknown prefactor `A`
cancels. This lets the report state how much faster a reaction runs when the
heat source is hotter, without claiming absolute rates.

## 5. Calorimetry (experiment)

The heat of solution of a salt is measured in a simple insulated calorimeter:

```
q_water = m_water * c_water * dT          (c_water = 4.184 J/(g*K))
dH_sol  = -q_water / n_solute             (kJ/mol)
```

The calorimeter constant `C_cal` is calibrated with a salt of known `dH_sol`
(e.g. KCl, `+17.2 kJ/mol`). See `docs/protocol.md`.

## 6. References

* CRC Handbook of Chemistry and Physics, 97th ed.
* NIST Chemistry WebBook, SRD 69.
* Chase, M. W. (1998), NIST-JANAF Thermochemical Tables.
