# Experiment protocol

Two experiments are described. **A** is the core measurement (heat of solution,
calorimetry). **B** is an optional dehydration check for the salt hydrates the
screening flags.

## A. Heat of solution (calorimetry)

**Purpose.** Measure `dH_sol` of several salts and connect an endothermic,
spontaneous process to `dG = dH - T dS`.

**Materials.** Polystyrene (styrofoam) cup + lid, thermometer (0.1 C resolution
or a digital probe), balance (0.01 g), graduated cylinder, stirring rod.
Salts: KCl (calibration), NH4NO3 (endothermic), CaCl2 (exothermic),
NaOH (exothermic, handle with care).

**Procedure (per salt).**
1. Add 50.0 g of distilled water to the cup. Record the initial temperature `T_i`
   every 15 s for 1 min to establish a stable baseline.
2. Weigh 5.00 g of the salt.
3. Add the salt, stir gently, and record the maximum (or minimum) temperature
   `T_f`. Use the extreme value, not the reading right after mixing.
4. Compute `dT = T_f - T_i` and feed it to the CLI:

```bash
thermobank calo --mass 5.00 --M 74.55 --water 50.0 --dT -3.2 --lit 17.2   # KCl
```

**Calibration.** Run KCl first. The measured `dH_sol` should be within ~15% of
`+17.2 kJ/mol`; if not, estimate `C_cal` (see `calorimetry.heat_capacity_of_calorimeter`)
and correct the others.

**Controls.** Same water mass, same starting temperature, same stirring, run each
salt three times, report the median and the spread.

**Safety.** NaOH and CaCl2 are corrosive / strongly exothermic; wear goggles and
add salt in small portions.

## B. Dehydration check (optional)

**Purpose.** Confirm that a screened salt hydrate loses water at the predicted
threshold temperature.

1. Weigh ~5 g of `CaSO4.2H2O` (gypsum) into a crucible; record mass.
2. Heat in an oven or on a hotplate at 80, 100, 130 C for 30 min each.
3. Cool in a desiccator and reweigh. Mass loss = released water.
4. Compare the onset temperature with `thermobank reaction gypsum`.
   Expect the mass loss to become appreciable near/above `T* = 88 C`; the gap
   between theory and observation is the kinetic barrier (see `kinetics`).

## Data to record

| salt | m_water (g) | m_salt (g) | T_i (C) | T_f (C) | dT (K) | dH_meas (kJ/mol) | dH_lit | error |
|---|---|---|---|---|---|---|---|---|
| KCl |  |  |  |  |  |  | +17.2 |  |
| NH4NO3 |  |  |  |  |  |  | +25.7 |  |
| CaCl2 |  |  |  |  |  |  | -82.8 |  |
