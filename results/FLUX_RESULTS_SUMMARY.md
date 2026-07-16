# Stem CH₄/CO₂ flux — results so far

*Working summary of the 2026 campaign. Fluxes from LI-7810 + goFlux (best.flux).
CH₄ reported as nmol m⁻² s⁻¹ and, in parentheses, µg CH₄ m⁻² h⁻¹ (×57.744);
CO₂ as µmol m⁻² s⁻¹.*

## Campaigns

| dataset | what | trees | measurements |
|---|---|---|---|
| **Climbing (intensive)** | flooded-forest stems, base → canopy (0.4–10 m) | 8 | 40 |
| **Basic (ground)** | flooded-forest stems, base only (0.4–1.6 m) | 51 | 128 |
| **BCI (terra firme)** | Barro Colorado I., 2 species | 2 spp | 12 |

## Headline: stems are strong CH₄ sources concentrated at the base

CH₄ emission is large at the waterline/base and falls steeply with height, reaching
near-ambient by ~1.6–2 m. Median CH₄ by height:

| height | climbing — nmol m⁻² s⁻¹ (µg m⁻² h⁻¹) | basic — nmol m⁻² s⁻¹ (µg m⁻² h⁻¹) |
|---|---|---|
| ~40 cm | **39.0 (2,254)** | **35.9 (2,073)** |
| ~80 cm | 6.7 (387) | 9.9 (569) |
| 160 cm | 0.4 (22) | 2.3 (132) |
| 500 cm | 0.0 (1) | — |
| 8–10 m (top) | 0.0 (2) | — |

- **Essentially all measurements are emission** (climbing 100 % positive; basic 96 %).
- Peak CH₄ ≈ 169 (climbing) / 445 (basic) nmol m⁻² s⁻¹ ≈ 9,800 / 25,700 µg m⁻² h⁻¹ at the base.
- Above ~1.6 m the stem is a negligible CH₄ source.

## Climbing and ground methods agree at the base

The two independent methods give near-identical base medians (~40 cm: 39 vs 36;
~80 cm: 6.7 vs 9.9 nmol m⁻² s⁻¹), cross-validating the ground-based approach against
the harnessed climbs on the 6 trees measured both ways.

## What controls the flux

Mixed model, `asinh(CH₄) ~ height + DBH + (1 | species/tree)` (basic, n = 128):

| term | effect (per SD) | t |
|---|---|---|
| height | **−0.61** | −5.2 |
| DBH | **+0.61** | +2.1 |

- **Height** is the dominant control — steep decline (above).
- **DBH** is a significant positive effect: larger-diameter stems emit more.
- Time-of-day and transect direction are not robust (transect is confounded with
  measurement date).

## Model form and the upscaling problem

Fitting CH₄ vs height on the raw scale (basic):

| model | AIC |
|---|---|
| power law  a·hᵇ | **1497** |
| exponential a·e^(kh) | 1497 |
| linear | 1501 |

- Power law and exponential are indistinguishable and both beat linear.
- Best-fit **power-law exponent b ≈ −1.32** (flux ∝ h⁻¹·³).
- **Consequence for upscaling:** with b < −1 the power law diverges as height → 0, so
  the whole-stem flux is dominated by the base/waterline — the one value we cannot
  currently measure (chamber is 20 cm tall; the forest is flooded, so nothing below
  the climber/boat). The base flux is therefore the key unconstrained quantity, and a
  direct waterline measurement is the highest-value next step.

## Contrast: BCI terra firme is CH₄-neutral

Non-flooded stems at BCI (Simarouba amara, Heisteria concinna) are **not** a CH₄
source: median CH₄ ≈ −0.006 nmol m⁻² s⁻¹ (≈ 0), range −0.02 to +0.01, only 25 %
positive — while still clearly respiring (CO₂ median 2.4 µmol m⁻² s⁻¹). This supports
a flooding/anaerobic-sediment origin for the transect emissions (CH₄ produced in
saturated soil, transported up and released through the lower stem) rather than
in-stem production.

## Isotope sampling (status)

δ¹³C vials collected at t0 (ambient) / t1 / t2 during the flux measurements:
**170 vials, 64 measurements, 15 trees.** All reconciled against the field sampling
sheet; concentrations estimated from the trace at each sampling time are compiled in
`results/sample_inventory_isotopes.csv` for the analyzer runs.

## Caveats / open items

- Base flux is extrapolated, not measured — the dominant term for upscaling is the
  least constrained. Recommend a waterline/near-surface measurement to confirm.
- Transect-direction differences are confounded with date; not interpreted.
- Flux estimates use goFlux `best.flux`; low-flux points near/below the detection
  floor (mostly high on the stem) contribute little to totals.
