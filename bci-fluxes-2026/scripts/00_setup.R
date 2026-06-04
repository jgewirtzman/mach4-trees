# =============================================================================
# 00_setup.R  (BCI fluxes 2026)
# Packages, paths, and constants for the BCI goFlux pipeline.
# Source this script at the top of every other script.
#
# Same processing as transect-trees-2026 (import.LI7810 -> auxfile ->
# click.peak2 stacked -> goFlux -> best.flux -> compile -> MDF/LOD).
# Differences for BCI:
#   * field notes have TWO heights per tree: "2m" and "Base"
#   * instrument logged in America/New_York; site is Panama (see Timezone below)
#   * chamber is a Blueflux-STYLE class B chamber (dimensions only; NOT a
#     Blueflux project) -- SET ITS GEOMETRY BELOW.
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
proj_dir    <- "/Users/jongewirtzman/My Drive/Research/mach4-trees/bci-fluxes-2026"
scripts_dir <- file.path(proj_dir, "scripts")
data_dir    <- file.path(proj_dir, "data")
rdata_dir   <- file.path(proj_dir, "RData")
results_dir <- file.path(proj_dir, "results")
plots_dir   <- file.path(proj_dir, "plots")

# Raw inputs (read-only)
licor_raw <- "/Users/jongewirtzman/My Drive/Research/mach4-trees/data/licor/TG10-01128-2026-05-01T000000.data"
field_csv <- "/Users/jongewirtzman/My Drive/Research/mach4-trees/data/field_notes/BCI Fluxes - Sheet1.csv"

# Staged date-subset of the LI-7810 file (campaign dates only; built in 01).
licor_subset <- file.path(data_dir, "li7810_measurement_dates.data")

for (d in c(data_dir, rdata_dir, results_dir, plots_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# --- Timezone ----------------------------------------------------------------
# The LI-7810 header reports Timezone: America/New_York, so the DATE/TIME columns
# are in that clock. We import AND build the auxfile in the same tz so the
# recorded HH:MM start times line up with the data. If the field notes were
# instead written in Panama local time (UTC-5, i.e. 1 h behind New_York in May),
# set field_offset_hours <- -1 (you'll also see the trace sit ~1 h off the blue
# start line in click.peak -- easy to confirm interactively).
tz_data            <- "America/New_York"
field_offset_hours <- 0    # add to field-note times to match the data clock

# --- Observation window (initial guess for click.peak2) ----------------------
# Show 3 min BEFORE and 12 min AFTER the recorded start time (same as transect):
#   before = shoulder              = 180 s = 3 min
#   after  = obs.length + shoulder = 720 s = 12 min  (=> obs.length = 540 s)
obs_length       <- 540   # seconds
shoulder_secs    <- 180   # seconds before start (and added after obs.length)
typical_obs_secs <- 180   # dashed ~3 min length guide

# =============================================================================
# CLASS C CHAMBER GEOMETRY -- from the Blueflux project (dimensions only; this
# is NOT a Blueflux project). Source files:
#   Blueflux/blueflux-ground/data/field_notes/dimension_csvs/
#     surface_area.csv      -> "C series" SA = 213.35 cm^2
#     simplified_volume.csv -> chamber "C" (S1), injection-dilution total volumes:
#         20230302: 1389.75 / 1270.10 / 1436.25 mL
#         20220907:  992.35 /  885.49 / 1007.67 mL   (mean total = 1163.6 mL)
#       pure chamber = mean total - LGR tubing(29) - LGR cell(70) = 1064.6 cm^3
#     additional_vol.csv    -> tubing = 29 cm^3 (Blueflux closed-loop kit)
# For THIS LI-7810 deployment, swap the analyzer cell to the LI-7810 (16 cm^3):
#   Vtot = pure_chamber(1064.6) + LI-7810 cell(16) + tubing(29) = 1109.6 cm^3.
# CAVEAT: the two Blueflux campaigns differ (~1.37 L vs ~0.96 L total); 1064.6
# is the mean. If the BCI chamber matched one campaign, set chamber_vol_cm3
# directly (1290 cm^3 for 2023-03, or 863 cm^3 for 2022-09, pure head-space).
# =============================================================================
CHAMBER_GEOMETRY_IS_PLACEHOLDER <- FALSE

chamber_area_cm2 <- 213.35     # class C "SA cm2"     (surface_area.csv)
chamber_vol_cm3  <- 1064.6     # class C pure head-space (Blueflux mean)

# Dead volume = LI-7810 optical cell + Blueflux closed-loop tubing
li7810_cell_cm3 <- 16          # LI-7810 optical cell
tubing_vol_cm3  <- 29          # Blueflux tubing dead volume (additional_vol.csv)
dead_vol_cm3    <- li7810_cell_cm3 + tubing_vol_cm3

Vtot_L <- (chamber_vol_cm3 + dead_vol_cm3) / 1000   # = 1.1096 L

# Chamber air temperature (C) and pressure (kPa) -> Tcham / Pcham
Tcham_default_C   <- 28.0    # <-- TODO BCI chamber air temperature (C)
Pcham_default_kPa <- 101.3   # <-- TODO BCI barometric pressure (kPa)

# --- Instrument precision (for goFlux / MDF) ---------------------------------
# c(CO2dry_ppm, CH4dry_ppb, H2O_ppm). LI-7810 spec; replace if you have measured.
li7810_prec <- c(3.5, 0.6, 45)

# --- best.flux model-selection criteria --------------------------------------
flux_criteria <- c("MAE", "AICc", "g.factor", "MDF")

# --- Campaign dates (used to subset the big .data file in 01_import.R) --------
measurement_dates <- c("2026-05-12")

message("Setup complete (BCI). Project: ", proj_dir)
if (isTRUE(CHAMBER_GEOMETRY_IS_PLACEHOLDER))
  message("  *** CHAMBER GEOMETRY IS A PLACEHOLDER -- set the class B dimensions ",
          "in 00_setup.R before 04 ***")
message(sprintf("  Chamber INPUTS: Vtot=%.3f L, Area=%.1f cm2, Tcham=%.1f C, Pcham=%.1f kPa",
                Vtot_L, chamber_area_cm2, Tcham_default_C, Pcham_default_kPa))
