# =============================================================================
# 03_manual_id.R  (BASIC trees)
# Interactive manual identification of measurement start/end times, CO2 (top) +
# CH4 (bottom) stacked via click.peak2.stacked(). Same procedure as
# transect-trees-2026; ~128 measurements -> done in batches of 20.
#
# *** RUN IN AN INTERACTIVE R SESSION (RStudio / R GUI). Does NOT run under
#     Rscript -- it needs a clickable plot device. ***
#
# For each measurement: TOP = CO2dry_ppm (click START then END), BOTTOM =
# CH4dry_ppb (reference). Blue line = recorded start; dashed grey = +typical.
#
# You can re-run a single batch by editing `batches` below; results accumulate
# in manID_batches (saved after each batch so you don't lose progress).
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/transect-basic-2026/scripts/00_setup.R"))
source(file.path(scripts_dir, "helpers", "click_peak2_stacked.R"))

load(file.path(rdata_dir, "licor_imported.RData"))   # licor
load(file.path(rdata_dir, "auxfile.RData"))          # aux, field_meta

dir.create(file.path(plots_dir, "click_peak"), recursive = TRUE, showWarnings = FALSE)

# --- Observation windows ------------------------------------------------------
ow <- obs.win(inputfile = licor, auxfile = aux,
              obs.length = obs_length, shoulder = shoulder_secs)
message(length(ow), " observation windows created")

# --- Interactive identification in batches of 20 -----------------------------
win_w <- 7; win_h <- 6
batch_size <- 20
n <- length(ow)
batches <- split(seq_len(n), ceiling(seq_len(n) / batch_size))

# Resume support: keep a per-batch list on disk.
mb_path <- file.path(rdata_dir, "manID_batches.RData")
if (file.exists(mb_path)) load(mb_path) else manID_batches <- vector("list", length(batches))

for (b in seq_along(batches)) {
  if (!is.null(manID_batches[[b]])) {            # already done -> skip
    message("Batch ", b, " already done; skipping. (delete it in manID_batches to redo)")
    next
  }
  message("=== Batch ", b, " of ", length(batches), " (measurements ",
          min(batches[[b]]), "-", max(batches[[b]]), ") ===")
  manID_batches[[b]] <- click.peak2.stacked(
    ow.list   = ow,
    seq       = batches[[b]],
    gastype   = "CO2dry_ppm",     # clicked panel (top)
    gastype2  = "CH4dry_ppb",     # reference panel (bottom)
    ref.secs  = typical_obs_secs,
    plot.lim  = c(300, 5000),
    plot.lim2 = c(1500, 100000),
    width = win_w, height = win_h,
    save.plots = file.path(plots_dir, "click_peak", sprintf("batch%02d.pdf", b))
  )
  save(manID_batches, file = mb_path)            # checkpoint after each batch
}

# --- Combine and save ---------------------------------------------------------
manID <- as.data.frame(dplyr::bind_rows(manID_batches))
save(manID, file = file.path(rdata_dir, "manID.RData"))
message("Manual ID complete: ", length(unique(manID$UniqueID)), " measurements")
message("Saved manID.RData. Proceed to 04_flux_calculation.R")
