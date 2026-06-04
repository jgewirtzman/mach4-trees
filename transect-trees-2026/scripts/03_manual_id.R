# =============================================================================
# 03_manual_id.R
# Interactive manual identification of measurement start/end times, showing CO2
# (top) and CH4 (bottom) stacked in one window via click.peak2.stacked()
# (helpers/click_peak2_stacked.R), a two-panel variant of goFlux::click.peak2().
#
# *** RUN THIS IN AN INTERACTIVE R SESSION (RStudio / R GUI) -- it needs a
#     clickable plot device. It does NOT run under Rscript. ***
#
# For each measurement a stacked plot pops up:
#   - TOP    = CO2dry_ppm   (you click on this panel)
#   - BOTTOM = CH4dry_ppb   (shown for reference; same window applied to it)
#   Blue lines = initial start/end from the auxfile.
#   1. Click ONCE on the true start of the rise (top/CO2 panel).
#   2. Click ONCE on the true end.
#   A stacked validation plot flashes (red = selected window); times are stored.
#
# Note: CO2dry_ppm / CH4dry_ppb are the LI-7810's native DRY mole fractions
# (the instrument applies the water-vapour correction internally); goFlux()
# additionally uses the H2O_ppm column for the head-space dilution term in 04.
#
# 40 windows -> do them in two batches of 20.
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/transect-trees-2026/scripts/00_setup.R"))
source(file.path(scripts_dir, "helpers", "click_peak2_stacked.R"))

load(file.path(rdata_dir, "licor_imported.RData"))   # licor
load(file.path(rdata_dir, "auxfile.RData"))          # aux, field_meta

dir.create(file.path(plots_dir, "click_peak"), recursive = TRUE,
           showWarnings = FALSE)

# --- Observation windows ------------------------------------------------------
# obs.length / shoulder are processing parameters from 00_setup.R (NOT the
# auxfile). They set the displayed trace: start - shoulder .. start + obs.length
# + shoulder  (here: 3 min before, 12 min after the recorded start time).
ow <- obs.win(
  inputfile  = licor,
  auxfile    = aux,
  obs.length = obs_length,
  shoulder   = shoulder_secs
)
message(length(ow), " observation windows created")

# --- Interactive identification (stacked CO2 + CH4), two batches of <= 20 -----
# Adjust plot.lim / plot.lim2 if the y-axis clips your traces.
# width/height = size of the pop-up plot window, in inches (shrink if too big).
win_w <- 7; win_h <- 6

manID.1 <- click.peak2.stacked(
  ow.list   = ow,
  seq       = seq(1, 20),
  gastype   = "CO2dry_ppm",     # clicked panel (top)
  gastype2  = "CH4dry_ppb",     # reference panel (bottom)
  ref.secs  = typical_obs_secs, # dashed ~3 min length guide
  plot.lim  = c(300, 5000),
  plot.lim2 = c(1500, 100000),
  width = win_w, height = win_h,
  save.plots= file.path(plots_dir, "click_peak", "batch1.pdf")
)

manID.2 <- click.peak2.stacked(
  ow.list   = ow,
  seq       = seq(21, length(ow)),
  gastype   = "CO2dry_ppm",
  gastype2  = "CH4dry_ppb",
  ref.secs  = typical_obs_secs,
  plot.lim  = c(300, 5000),
  plot.lim2 = c(1500, 100000),
  width = win_w, height = win_h,
  save.plots= file.path(plots_dir, "click_peak", "batch2.pdf")
)

# --- Combine and save ---------------------------------------------------------
manID <- dplyr::bind_rows(manID.1, manID.2)
save(manID, file = file.path(rdata_dir, "manID.RData"))

message("Manual ID complete: ", length(unique(manID$UniqueID)),
        " measurements identified")
message("Saved: ", file.path(rdata_dir, "manID.RData"))
message("Proceed to 04_flux_calculation.R")
