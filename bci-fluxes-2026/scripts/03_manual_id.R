# =============================================================================
# 03_manual_id.R  (BCI)
# Interactive manual identification of measurement start/end times, showing CO2
# (top) and CH4 (bottom) stacked via click.peak2.stacked()
# (helpers/click_peak2_stacked.R). Same procedure as transect-trees-2026.
#
# *** RUN THIS IN AN INTERACTIVE R SESSION (RStudio / R GUI) -- it needs a
#     clickable plot device. It does NOT run under Rscript. ***
#
# For each measurement a stacked plot pops up:
#   - TOP    = CO2dry_ppm   (you click on this panel)
#   - BOTTOM = CH4dry_ppb   (shown for reference; same window applied to it)
#   Blue line  = recorded start (from auxfile); dashed grey = +typical length.
#   1. Click the true START of the rise (top/CO2 panel).
#   2. Click the true END.
#
# 12 measurements (6 trees x {2m, Base}) -> one batch.
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/bci-fluxes-2026/scripts/00_setup.R"))
source(file.path(scripts_dir, "helpers", "click_peak2_stacked.R"))

load(file.path(rdata_dir, "licor_imported.RData"))   # licor
load(file.path(rdata_dir, "auxfile.RData"))          # aux, field_meta

dir.create(file.path(plots_dir, "click_peak"), recursive = TRUE,
           showWarnings = FALSE)

# --- Observation windows ------------------------------------------------------
ow <- obs.win(
  inputfile  = licor,
  auxfile    = aux,
  obs.length = obs_length,
  shoulder   = shoulder_secs
)
message(length(ow), " observation windows created")

# --- Interactive identification (stacked CO2 + CH4) --------------------------
# Adjust plot.lim / plot.lim2 if the y-axis clips your traces.
# width/height = size of the pop-up plot window, in inches (shrink if too big).
win_w <- 7; win_h <- 6

manID <- click.peak2.stacked(
  ow.list   = ow,
  seq       = seq(1, length(ow)),
  gastype   = "CO2dry_ppm",     # clicked panel (top)
  gastype2  = "CH4dry_ppb",     # reference panel (bottom)
  ref.secs  = typical_obs_secs, # dashed ~3 min length guide
  plot.lim  = c(300, 5000),
  plot.lim2 = c(1500, 100000),
  width = win_w, height = win_h,
  save.plots= file.path(plots_dir, "click_peak", "batch1.pdf")
)

save(manID, file = file.path(rdata_dir, "manID.RData"))
message("Manual ID complete: ", length(unique(manID$UniqueID)),
        " measurements identified")
message("Saved: ", file.path(rdata_dir, "manID.RData"))
message("Proceed to 04_flux_calculation.R")
