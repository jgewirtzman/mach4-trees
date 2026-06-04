# =============================================================================
# 00_setup.R  (Transect BASIC trees 2026)
# Packages, paths, constants for the basic-tree goFlux pipeline.
# Source this at the top of every other script.
#
# Same campaign/instrument/chamber as transect-trees-2026 (the CLIMBING trees);
# this project processes the BASIC trees (3 heights: ~40/80/160 cm) from the
# hand-entered datasheets (convert/parsed/all_basic_sheets.csv).
#
# Instrument: LI-COR LI-7810 (CO2 / CH4 / H2O) -> goFlux::import.LI7810()
# Manual identification (click.peak2) is INTERACTIVE and done by the user.
# =============================================================================

# --- Locale ------------------------------------------------------------------
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
proj_dir    <- "/Users/jongewirtzman/My Drive/Research/mach4-trees/transect-basic-2026"
scripts_dir <- file.path(proj_dir, "scripts")
data_dir    <- file.path(proj_dir, "data")
rdata_dir   <- file.path(proj_dir, "RData")
results_dir <- file.path(proj_dir, "results")
plots_dir   <- file.path(proj_dir, "plots")

# Raw inputs (read-only) -- SAME licor file as transect-trees-2026
licor_raw <- "/Users/jongewirtzman/My Drive/Research/mach4-trees/data/licor/TG10-01503-2026-05-01T000000.data"
# Basic-tree field data = the transcribed datasheet master
field_csv <- "/Users/jongewirtzman/My Drive/Research/mach4-trees/data/datasheets/parsed/all_basic_sheets.csv"

licor_subset <- file.path(data_dir, "li7810_measurement_dates.data")
for (d in c(data_dir, rdata_dir, results_dir, plots_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# --- Timezone (same convention as transect-trees-2026: NO offset) ------------
tz_data <- "UTC"

# --- Observation window (initial guess for click.peak2) ----------------------
# Show 3 min BEFORE and 12 min AFTER the recorded start time (same as transect).
obs_length       <- 540   # s
shoulder_secs    <- 180   # s before start (and added after obs.length)
typical_obs_secs <- 180   # dashed ~3 min length guide

# --- Chamber geometry: SAME as transect-trees-2026 ---------------------------
# Chamber inner dimensions 20 x 30 x 2 cm.
chamber_vol_cm3  <- 20 * 30 * 2     # 1200 cm^3 head-space
chamber_area_cm2 <- 20 * 30         # 600 cm^2 enclosed bark footprint
li7810_cell_cm3 <- 16               # LI-7810 optical cell
tubing_len_cm   <- 2 * 300          # 2 x 3 m intake + return
tubing_id_cm    <- 0.3175           # 1/8 inch ID
tubing_vol_cm3  <- pi * (tubing_id_cm / 2)^2 * tubing_len_cm
dead_vol_cm3    <- li7810_cell_cm3 + tubing_vol_cm3
Vtot_L <- (chamber_vol_cm3 + dead_vol_cm3) / 1000   # ~1.2635 L

Tcham_default_C   <- 28.0    # chamber air temperature (C)
Pcham_default_kPa <- 101.3   # barometric pressure (kPa)

# --- Instrument precision (LI-7810) ------------------------------------------
li7810_prec <- c(3.5, 0.6, 45)   # c(CO2dry_ppm, CH4dry_ppb, H2O_ppm)

# --- best.flux model-selection criteria --------------------------------------
flux_criteria <- c("MAE", "AICc", "g.factor", "MDF")

# --- Campaign dates (subset the big .data file in 01_import.R) ----------------
measurement_dates <- c("2026-05-22","2026-05-23","2026-05-24",
                       "2026-05-25","2026-05-26","2026-05-27")

message("Setup complete (BASIC trees). Project: ", proj_dir)
message(sprintf("  Chamber: Vtot=%.3f L, Area=%.1f cm2, Tcham=%.1f C, Pcham=%.1f kPa",
                Vtot_L, chamber_area_cm2, Tcham_default_C, Pcham_default_kPa))
