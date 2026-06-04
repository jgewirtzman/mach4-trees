# =============================================================================
# 00_setup.R
# Packages, paths, and constants for the Transect Trees 2026 goFlux pipeline.
# Source this script at the top of every other script.
#
# Structure mirrors:
#   whole_tree_flux/data processing/goFlux_reprocessing  (import -> auxfile ->
#       click.peak2 -> goFlux -> best.flux -> compile)
#   ch4-data-filtering/scripts                            (MDF / LOD step)
#
# Instrument: LI-COR LI-7810 (CO2 / CH4 / H2O) -> goFlux::import.LI7810()
# Manual identification (click.peak2) is INTERACTIVE and done by the user.
# =============================================================================

# --- Locale ------------------------------------------------------------------
# Rscript can default to the C locale, which mangles accented species names
# (Macucú, Ingá, ...) into "<U+00FA>" escapes when writing CSVs. Force UTF-8.
try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)

# --- Packages ----------------------------------------------------------------

if (!require("remotes", quietly = TRUE)) install.packages("remotes")
if (!require("goFlux", quietly = TRUE)) remotes::install_github("Qepanna/goFlux")
suppressMessages(library(goFlux))

pkgs <- c("dplyr", "tidyr", "readr", "stringr", "lubridate", "purrr", "openxlsx")
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) install.packages(pkg)
  suppressMessages(library(pkg, character.only = TRUE))
}

# --- Paths -------------------------------------------------------------------

proj_dir    <- "/Users/jongewirtzman/My Drive/Research/mach4-trees/transect-trees-2026"
scripts_dir <- file.path(proj_dir, "scripts")
data_dir    <- file.path(proj_dir, "data")
rdata_dir   <- file.path(proj_dir, "RData")
results_dir <- file.path(proj_dir, "results")
plots_dir   <- file.path(proj_dir, "plots")

# Raw inputs (read-only)
licor_raw <- "/Users/jongewirtzman/My Drive/Research/mach4-trees/data/licor/TG10-01503-2026-05-01T000000.data"
field_csv <- "/Users/jongewirtzman/My Drive/Research/mach4-trees/data/field_notes/Transect Trees - 052026.csv"

# Staged date-subset of the LI-7810 file (campaign dates only; built in 01).
# This is date-windowing for import speed only -- NO gas values are altered.
licor_subset <- file.path(data_dir, "li7810_measurement_dates.data")

for (d in c(data_dir, rdata_dir, results_dir, plots_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# --- Timezone ----------------------------------------------------------------
# LI-7810 logs in UTC and the field-note times use the same clock: NO offset.
tz_data <- "UTC"

# --- Observation window (initial guess for click.peak2) ----------------------
# Only START times were recorded. The CONTEXT shown around each trace is
#   start - shoulder  ...  start + obs.length + shoulder
# and you click the true start/end anywhere in that range (length is whatever
# you select -- per measurement). Set to show 3 min BEFORE and 12 min AFTER:
#   before = shoulder              = 180 s = 3 min
#   after  = obs.length + shoulder = 720 s = 12 min  (=> obs.length = 540 s)
obs_length    <- 540   # seconds; with shoulder below -> context ends at start + 12 min
shoulder_secs <- 180   # seconds before the start (and added after obs.length) = 3 min

# Typical measurement length, drawn as a light dashed GUIDE line (not a limit) so
# you have a visual anchor for "~3 min" while clicking. Measurements vary.
typical_obs_secs <- 180   # 3 min guide

# =============================================================================
# >>> REQUIRED INPUTS NOT PRESENT IN THE TWO RAW FILES -- SET THESE <<<
# These are not in the LI-7810 file or the field notes. Placeholders below let
# the pipeline run; replace with the real campaign values. Flux scales linearly
# with Vtot/Area, so these directly affect the numbers.
# =============================================================================

# Chamber head-space + analyzer-cell + tubing dead volume -> Vtot (LITERS).
# Chamber inner dimensions: 20 x 30 x 2 cm.
chamber_vol_cm3  <- 20 * 30 * 2     # 1200 cm^3 head-space
chamber_area_cm2 <- 20 * 30         # 600 cm^2 enclosed bark footprint

# Dead volume = LI-7810 optical cell + closed-loop tubing.
li7810_cell_cm3 <- 16               # LI-7810 optical cell (~16 cm^3; confirm if exact)
tubing_len_cm   <- 2 * 300          # 2 x 3 m intake + return = 600 cm
tubing_id_cm    <- 0.3175           # 1/8 inch inner diameter (3.175 mm)
tubing_vol_cm3  <- pi * (tubing_id_cm / 2)^2 * tubing_len_cm
dead_vol_cm3    <- li7810_cell_cm3 + tubing_vol_cm3

Vtot_L <- (chamber_vol_cm3 + dead_vol_cm3) / 1000

# Chamber air temperature (C) and pressure (kPa) -> Tcham / Pcham
# (whole_tree_flux pulled these from Harvard Forest met data; supply your own
#  site values or a per-measurement series.)
Tcham_default_C   <- 28.0   # <-- TODO site/chamber air temperature (C)
Pcham_default_kPa <- 101.3  # <-- TODO site barometric pressure (kPa)

# --- Instrument precision (for goFlux / MDF) ---------------------------------
# c(CO2dry_ppm, CH4dry_ppb, H2O_ppm). import.LI7810 default is the LI-7810
# spec below; replace with your measured 1-sigma precision if you have it.
# (whole_tree_flux used UGGA: ugga_prec <- c(0.2, 1.4, 50).)
li7810_prec <- c(3.5, 0.6, 45)   # <-- confirm / replace with measured precision

# --- best.flux model-selection criteria --------------------------------------
# (same family as whole_tree_flux / ch4-data-filtering)
flux_criteria <- c("MAE", "AICc", "g.factor", "MDF")

# --- Campaign dates (used to subset the big .data file in 01_import.R) --------
measurement_dates <- c("2026-05-23", "2026-05-24", "2026-05-25")

message("Setup complete. Project: ", proj_dir)
message(sprintf("  Chamber INPUTS (edit 00_setup.R): Vtot=%.3f L, Area=%.1f cm2, Tcham=%.1f C, Pcham=%.1f kPa",
                Vtot_L, chamber_area_cm2, Tcham_default_C, Pcham_default_kPa))
