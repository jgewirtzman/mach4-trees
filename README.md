# mach4-trees

Tree-stem CH₄ and CO₂ flux processing (goFlux / LI-COR LI-7810), consolidated.
Three parallel pipelines share the same structure (`import → auxfile →
click.peak2 → goFlux → best.flux → compile → MDF/LOD`):

| folder | what | chamber | licor file |
|---|---|---|---|
| `transect-trees-2026/` | Mamirauá **climbing** trees (5 heights: 40/80/160/500/1000 cm) | 20×30×2 cm | `TG10-01503` |
| `transect-basic-2026/`  | Mamirauá **basic** trees (3 heights: ~40/80/160 cm), from transcribed datasheets | 20×30×2 cm | `TG10-01503` |
| `bci-fluxes-2026/`       | BCI (Panama) trees (2 heights: Base + 2 m) | Blueflux class C | `TG10-01128` |

## Layout

```
mach4-trees/
├── data/
│   ├── licor/                LI-7810 .data files            (git-ignored, large)
│   ├── field_notes/          field CSVs + Transect Trees.xlsx
│   └── datasheets/
│       ├── photos/           basic-sheet photos             (git-ignored, large)
│       └── parsed/           transcribed datasheet CSVs (fixed/ = verified;
│                             all_basic_sheets.csv = master; paste_ready_* = workbook merge)
├── transect-trees-2026/      scripts/ results/ plots/ (RData/ git-ignored)
├── transect-basic-2026/
└── bci-fluxes-2026/
```

## Running a pipeline

Each project's `scripts/` is numbered. From R/RStudio:

1. `01_import.R`  — subset the LI-7810 file to campaign dates + `import.LI7810()`
2. `02_build_auxfile.R` — build the goFlux auxfile from the field notes
3. `03_manual_id.R` — **interactive** (RStudio/R GUI): click start/end on stacked CO₂/CH₄ traces
4. `04_flux_calculation.R` → `05_compile_results.R` → `06_mdf_lod.R` (Rscript-able)

All paths are set in each `scripts/00_setup.R` (sourced by every step). Fluxes:
CO₂ in µmol m⁻² s⁻¹, CH₄ in nmol m⁻² s⁻¹.

## Notes

- **Git inside Google Drive:** this repo lives in My Drive. Drive's file-by-file
  sync can occasionally disturb the `.git` directory — if git acts up, pause
  Drive sync during git operations, or push to a remote as the real backup.
- **Large files are git-ignored** (the LI-7810 `.data` files and datasheet
  photos). They're present locally under `data/` but not in version control.
- `RData/` and the date-subset `.data` are regenerable (re-run `01`–`04`).
