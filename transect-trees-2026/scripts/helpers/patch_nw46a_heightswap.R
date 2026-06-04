# =============================================================================
# helpers/patch_nw46a_heightswap.R
# One-off correction: T-NW46a's 40 cm and 160 cm height labels were reversed in
# the original 5/24 field notes (confirmed by the 5/26 remeasurement: the steep
# 19:29 window is the 40 cm base, the flat 19:45 window is 160 cm).
#
# This swaps the labels in the PRODUCED artifacts only (manID + auxfile RData) --
# NOT by re-deriving from the source CSV (the local Downloads copy still has the
# old order). After running this, re-run 04 -> 07.
#
# NOTE: do NOT re-run 02 before 04, or it will overwrite this patch from the
# stale CSV. For a durable fix, also correct the heights in the source field
# notes ("40;80;160;500;1000" -> "160;80;40;500;1000").
# =============================================================================

proj <- "/Users/jongewirtzman/My Drive/Research/mach4-trees/transect-trees-2026"
load(file.path(proj, "RData", "auxfile.RData"))   # aux, field_meta
load(file.path(proj, "RData", "manID.RData"))      # manID

A <- "T-NW46a_40cm"; B <- "T-NW46a_160cm"
swap <- function(v) ifelse(v == A, B, ifelse(v == B, A, v))

# manID: relabel the clicked windows (gas data/flag stay; 19:29 steep -> 40 cm)
manID$UniqueID <- swap(manID$UniqueID)

# aux: relabel (start.time travels with each row, so "40 cm" -> 19:29)
aux$UniqueID <- swap(aux$UniqueID)

# field_meta: swap UniqueID AND height_cm together, keeping time_str/start.time
# on each physical row (the 19:29 row becomes 40 cm; the 19:45 row becomes 160 cm)
was40  <- field_meta$UniqueID == A
was160 <- field_meta$UniqueID == B
field_meta$UniqueID[was40]   <- B; field_meta$height_cm[was40]   <- 160
field_meta$UniqueID[was160]  <- A; field_meta$height_cm[was160]  <- 40

save(aux, field_meta, file = file.path(proj, "RData", "auxfile.RData"))
save(manID,           file = file.path(proj, "RData", "manID.RData"))

cat("Patched T-NW46a 40cm <-> 160cm in manID + auxfile RData.\n")
cat("Now re-run 04_flux_calculation.R -> 07_plots.R (NOT 02).\n")
