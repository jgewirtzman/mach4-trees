# Transect Trees 2026 — goFlux processing pipeline

Calculates CO2 and CH4 fluxes from LI-COR **LI-7810** stem-chamber data using
the **goFlux** R package, following the workflow conventions in
`whole_tree_flux/data processing/goFlux_reprocessing` (import → auxfile →
interactive `click.peak2` → `goFlux` → `best.flux` → compile) and the MDF/LOD
methods in `ch4-data-filtering/scripts/03_mdf_computation.R`.

## Inputs
- LI-7810 1 Hz data: `~/Downloads/TG10-01503-2026-05-01T000000.data`
- Field notes:       `~/Downloads/Transect Trees - 052026.csv`
  (8 measured trees × 5 stem heights = 40 measurements; `fluxStartNotes_2026_05`
  gives each measurement's **start** time, `heights_2026_05` the stem height in cm.
  Times are UTC — same clock as the LI-7810, **no offset**.)

## Scripts (run in order)
| # | script | interactive? | what it does |
|---|--------|:---:|---|
| 00 | `00_setup.R` | — | packages, paths, **constants you must set** (chamber, T, P) |
| 01 | `01_import.R` | no (`Rscript`) | subset .data to campaign dates, `import.LI7810()` |
| 02 | `02_build_auxfile.R` | no | parse field notes → goFlux auxfile (40 windows) |
| 03 | `03_manual_id.R` | **YES — RStudio** | stacked **CO2 (top) + CH4 (bottom)** click window (`helpers/click_peak2_stacked.R`): click true start/end on the CO2 panel; same window applied to CH4 (two batches of 20) |
| 04 | `04_flux_calculation.R` | no | `goFlux()` CO2 & CH4 + `best.flux()` (LM/HM) |
| 05 | `05_compile_results.R` | no | merge fluxes + field metadata → tidy table |
| 06 | `06_mdf_lod.R` | no | MDF/LOD: Manufacturer/goFlux, Wassmann 2018, Christiansen 2015 × 99/95/90% |

Steps 01, 02, 04, 05, 06 run with `Rscript`. **Step 03 must run in an
interactive R session** (RStudio/R GUI) because it needs a clickable plot.

## >>> Additional data you need to supply (in `00_setup.R`) <<<
None of these are in the LI-7810 file or the field notes. Flux **and** MDF scale
with the flux term `V·P/(A·R·T)`, so they directly change the numbers:

1. **Chamber head-space volume** `chamber_vol_L` (L)
2. **Analyzer + tubing dead volume** `vtot_addition` (L) — LI-7810 optical cell
   (~0.016 L) + your intake tubing volume → together give `Vtot`
3. **Enclosed bark surface area** `chamber_area_cm2` (cm²)
4. **Chamber/air temperature** `Tcham_default_C` (°C) — a site value or a
   per-measurement series
5. **Barometric pressure** `Pcham_default_kPa` (kPa) — site value or from elevation

If a different chamber was used at different heights/trees, replace the single
constants with a per-measurement lookup joined into the auxfile in `02`.

Already handled (you do **not** need to provide):
- **Dry mole fractions** — `CO2dry_ppm` / `CH4dry_ppb` are the LI-7810's native
  *dry* (water-corrected) outputs; the instrument applies the water-vapour
  dilution + band-broadening correction internally. `goFlux()` additionally uses
  the `H2O_ppm` column (~2.9% vol here) for the head-space dilution term — so the
  water correction is applied once, not double-counted.
- **Instrument precision** for MDF — uses the LI-7810 values from your
  `ch4-data-filtering` project (CO2 3.5 ppm, CH4 0.6 ppb), = goFlux defaults.
- **Measurement end times** — you set them by clicking in step 03.

## Outputs (`results/`)
- `auxfile.csv`, `field_metadata.csv`
- `CO2_best.xlsx`, `CH4_best.xlsx`
- `transect_trees_2026_fluxes.csv` / `.xlsx` — compiled fluxes + goFlux MDF/g.fact
- `transect_trees_2026_fluxes_with_mdf.csv` / `.xlsx` — + all MDF approaches & below-MDF flags
- `mdf_lod_table.csv`
- `plots/click_peak/` — saved click.peak2 QC plots
