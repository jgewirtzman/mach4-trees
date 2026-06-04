# =============================================================================
# 03_manual_id.R  (BASIC trees) -- RESUMABLE by UniqueID
# Interactive manual identification, CO2 (top) + CH4 (bottom) stacked.
#
# *** RUN IN AN INTERACTIVE R SESSION (RStudio / R GUI). ***
#
# Resume model (robust to re-dating / re-ordering of the auxfile):
#   - Already-clicked measurements are kept in RData/manID_done.RData, keyed by
#     UniqueID. This script SKIPS any UniqueID already in manID_done and only
#     prompts for the ones still missing.
#   - New clicks accumulate in RData/manID_new_batches.RData (checkpointed after
#     every batch), then are merged with manID_done into manID.RData.
#
# NOTE: save.plots is OFF here. A failed plot-write previously aborted a batch
# and lost its clicks; skipping the PDF removes that failure point. The live
# validation flash still shows. (Re-make validation PDFs later from manID.)
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/transect-basic-2026/scripts/00_setup.R"))
source(file.path(scripts_dir, "helpers", "click_peak2_stacked.R"))

load(file.path(rdata_dir, "licor_imported.RData"))   # licor
load(file.path(rdata_dir, "auxfile.RData"))          # aux, field_meta

# --- Load already-done clicks (keyed by UniqueID) ----------------------------
done_path <- file.path(rdata_dir, "manID_done.RData")
if (file.exists(done_path)) { load(done_path) } else { manID_done <- NULL }
done_ids <- unique(manID_done$UniqueID)
message(length(done_ids), " measurements already clicked (from manID_done.RData)")

# --- Observation windows on the CURRENT auxfile ------------------------------
ow <- obs.win(inputfile = licor, auxfile = aux,
              obs.length = obs_length, shoulder = shoulder_secs)
ow_ids <- vapply(ow, function(x) unique(as.character(x$UniqueID))[1], character(1))

# --- Which measurements still need clicking ----------------------------------
todo <- which(!(ow_ids %in% done_ids))
message("To click now: ", length(todo), " measurements")
print(ow_ids[todo])

# --- Click the remaining ones, in batches of 20 ------------------------------
win_w <- 7; win_h <- 6; batch_size <- 20
todo_batches <- split(todo, ceiling(seq_along(todo) / batch_size))

nb_path <- file.path(rdata_dir, "manID_new_batches.RData")
if (file.exists(nb_path)) load(nb_path) else manID_new_batches <- vector("list", length(todo_batches))
if (length(manID_new_batches) != length(todo_batches))            # todo changed -> reset
  manID_new_batches <- vector("list", length(todo_batches))

for (b in seq_along(todo_batches)) {
  if (!is.null(manID_new_batches[[b]])) {
    message("New batch ", b, " already done; skipping."); next
  }
  message("=== NEW batch ", b, " of ", length(todo_batches),
          " (", length(todo_batches[[b]]), " measurements) ===")
  manID_new_batches[[b]] <- click.peak2.stacked(
    ow.list   = ow,
    seq       = todo_batches[[b]],
    gastype   = "CO2dry_ppm", gastype2 = "CH4dry_ppb",
    ref.secs  = typical_obs_secs,
    plot.lim  = c(300, 5000), plot.lim2 = c(1500, 100000),
    width = win_w, height = win_h,
    save.plots = NULL)                       # OFF: avoids the abort-and-lose-clicks bug
  save(manID_new_batches, file = nb_path)    # checkpoint after every batch
}

# --- Merge done + new -> manID -----------------------------------------------
manID <- as.data.frame(dplyr::bind_rows(manID_done, dplyr::bind_rows(manID_new_batches)))
save(manID, file = file.path(rdata_dir, "manID.RData"))
message("Manual ID complete: ", length(unique(manID$UniqueID)), " / ", length(ow),
        " measurements. Saved manID.RData. Proceed to 04_flux_calculation.R")
