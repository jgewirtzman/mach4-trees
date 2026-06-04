# =============================================================================
# 13_sample_inventory.R
# Consolidated isotope SAMPLE INVENTORY for the tech running the vials.
# One row per physical sample (vial): t0 (ambient), t1, t2 for every flux that
# has isotopes, across BOTH datasets (climbing + basic). For each sample it
# gives the clock time, the estimated CO2 & CH4 concentration in the vial, and
# whether/how much to dilute, plus the measurement's flux rate for context.
#
# Dilution: a vial is flagged if its concentration exceeds the analyzer ceiling
# below. Suggested factor = ceil(conc / ceiling). EDIT these two numbers to
# match your instrument's working range, then re-run.
CO2_CEIL_PPM <- 2000      # flag CO2 vials above this (ppm)
CH4_CEIL_PPM <- 12        # flag CH4 vials above this (ppm)
#
# Writes: results/sample_inventory_isotopes.csv          (long: one row/vial)
#         results/sample_inventory_isotopes_wide.csv     (one row/measurement)
# Run with Rscript.
# =============================================================================

suppressMessages({ library(dplyr); library(tidyr) })
base <- "/Users/jongewirtzman/My Drive/Research/mach4-trees"
unescape_u <- function(x) vapply(x, function(s) {
  if (is.na(s) || !grepl("<U\\+[0-9A-Fa-f]{4,6}>", s)) return(s)
  for (cp in regmatches(s, gregexpr("<U\\+[0-9A-Fa-f]{4,6}>", s))[[1]]) {
    hex <- sub("<U\\+([0-9A-Fa-f]+)>", "\\1", cp); s <- sub(cp, intToUtf8(strtoi(hex, 16L)), s, fixed = TRUE) }
  enc2utf8(s) }, character(1), USE.NAMES = FALSE)

read_iso <- function(path, dataset) {
  f <- read.csv(path, stringsAsFactors = FALSE, encoding = "UTF-8")
  f <- f[f$isotope %in% c(TRUE, "TRUE"), ]
  f$dataset <- dataset
  f$start_posix <- as.POSIXct(f$start.time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  f
}
fb <- read_iso(file.path(base, "transect-basic-2026/results/transect_basic_2026_fluxes_isotopes.csv"), "basic")
fc <- read_iso(file.path(base, "transect-trees-2026/results/transect_trees_2026_fluxes_isotopes.csv"), "climbing")

# --- explode to one row per existing vial (t0 / t1 / t2) ---------------------
to_long <- function(f) {
  rows <- lapply(seq_len(nrow(f)), function(i) {
    r <- f[i, ]
    stages <- c("t0","t1","t2")
    emin   <- c(if (is.finite(r$iso_t0_min)) r$iso_t0_min else -0.83, r$iso_t1_min, r$iso_t2_min)
    co2    <- c(r$iso_t0_CO2_ppm, r$iso_t1_CO2_ppm, r$iso_t2_CO2_ppm)
    ch4    <- c(r$iso_t0_CH4_ppb, r$iso_t1_CH4_ppb, r$iso_t2_CH4_ppb)
    keep   <- is.finite(co2) | is.finite(ch4)          # a vial exists only if we have a value
    data.frame(
      dataset = r$dataset, TreeID = sub("^T-", "", r$TreeID),
      species = unescape_u(r$species), height_cm = r$height_cm,
      timepoint = stages[keep],
      date = format(r$start_posix, "%Y-%m-%d"),
      elapsed_min = round(emin[keep], 2),
      sample_clock = format(r$start_posix + emin[keep]*60, "%Y-%m-%d %H:%M:%S"),
      CO2_ppm = round(co2[keep]), CH4_ppm = round(ch4[keep]/1000, 2),
      CO2_flux_umol_m2_s = round(r$CO2_best.flux, 3),
      CH4_flux_nmol_m2_s = round(r$CH4_best.flux, 2),
      CH4_flux_ug_m2_h   = round(r$CH4_best.flux * 57.744, 1),
      review = ifelse(is.na(r$iso_review), "", r$iso_review),
      note = ifelse(is.na(r$iso_flag), "", r$iso_flag),
      stringsAsFactors = FALSE)
  })
  bind_rows(rows)
}
inv <- bind_rows(to_long(fb), to_long(fc))

# --- sample label + dilution guidance ----------------------------------------
inv <- inv %>%
  mutate(
    # campaign tag keeps IDs unique for trees sampled in BOTH campaigns
    Sample_ID = paste0(TreeID, "_", height_cm, "cm_", timepoint, "_",
                       ifelse(dataset == "climbing", "C", "B")),
    dilute_CO2 = ifelse(is.finite(CO2_ppm) & CO2_ppm > CO2_CEIL_PPM, "YES", "no"),
    dilute_CH4 = ifelse(is.finite(CH4_ppm) & CH4_ppm > CH4_CEIL_PPM, "YES", "no"),
    dil_factor_CO2 = ifelse(is.finite(CO2_ppm), pmax(1, ceiling(CO2_ppm / CO2_CEIL_PPM)), NA),
    dil_factor_CH4 = ifelse(is.finite(CH4_ppm), pmax(1, ceiling(CH4_ppm / CH4_CEIL_PPM)), NA),
    # if both gases are run from one vial, dilution must satisfy the larger factor
    dil_factor_max = pmax(dil_factor_CO2, dil_factor_CH4, na.rm = TRUE),
    DILUTE = ifelse(dilute_CO2 == "YES" | dilute_CH4 == "YES", "YES", "no")
  ) %>%
  arrange(dataset, TreeID, height_cm, match(timepoint, c("t0","t1","t2"))) %>%
  select(Sample_ID, dataset, date, TreeID, species, height_cm, timepoint,
         sample_clock, elapsed_min, CO2_ppm, CH4_ppm,
         DILUTE, dil_factor_max, dilute_CO2, dil_factor_CO2, dilute_CH4, dil_factor_CH4,
         CO2_flux_umol_m2_s, CH4_flux_nmol_m2_s, CH4_flux_ug_m2_h, review, note)

out_long <- file.path(base, "results", "sample_inventory_isotopes.csv")
dir.create(dirname(out_long), showWarnings = FALSE)
write.csv(inv, out_long, row.names = FALSE, fileEncoding = "UTF-8")

# --- wide: one row per measurement (the per-tree view) -----------------------
wide <- inv %>%
  select(dataset, TreeID, species, height_cm, timepoint,
         sample_clock, CO2_ppm, CH4_ppm, CH4_flux_ug_m2_h) %>%
  pivot_wider(names_from = timepoint,
              values_from = c(sample_clock, CO2_ppm, CH4_ppm),
              names_glue = "{timepoint}_{.value}") %>%
  arrange(dataset, TreeID, height_cm)
out_wide <- file.path(base, "results", "sample_inventory_isotopes_wide.csv")
write.csv(wide, out_wide, row.names = FALSE, fileEncoding = "UTF-8")

# --- console summary ---------------------------------------------------------
cat(sprintf("Sample inventory: %d vials across %d measurements (%d trees).\n",
            nrow(inv), nrow(fb)+nrow(fc), length(unique(inv$TreeID))))
cat(sprintf("  by timepoint:  t0=%d  t1=%d  t2=%d\n",
            sum(inv$timepoint=="t0"), sum(inv$timepoint=="t1"), sum(inv$timepoint=="t2")))
cat(sprintf("  need dilution (ceil CO2>%d ppm / CH4>%d ppm):  CO2 vials=%d  CH4 vials=%d\n",
            CO2_CEIL_PPM, CH4_CEIL_PPM, sum(inv$dilute_CO2=="YES"), sum(inv$dilute_CH4=="YES")))
cat("  max dilution factor needed:  CO2 x", max(inv$dil_factor_CO2, na.rm=TRUE),
    "  CH4 x", max(inv$dil_factor_CH4, na.rm=TRUE), "\n")
cat("\nWrote:\n  ", out_long, "\n  ", out_wide, "\n")
