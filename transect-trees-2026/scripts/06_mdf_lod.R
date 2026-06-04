# =============================================================================
# 06_mdf_lod.R
# Minimum Detectable Flux (MDF) / Limit of Detection (LOD) analysis, reproducing
# the methodology in:
#   whole_tree_flux/.../09_mdf_lod_comparison.R
#   ch4-data-filtering/scripts/03_mdf_computation.R
#
# Three MDF definitions, each at 99 / 95 / 90% confidence, in flux units
# (CO2 umol m-2 s-1, CH4 nmol m-2 s-1). For each: MDF = (noise) / t * flux.term
#   1. Manufacturer / goFlux : noise = manufacturer 1-sigma precision (no mult).
#                              This equals goFlux's native best.flux $MDF column.
#   2. Wassmann et al. 2018  : noise = z * global empirical SD
#                              (z = qnorm(.995/.975/.95); global SD = median of
#                               per-measurement Allan deviations -> instrument-specific).
#   3. Christiansen 2015      : noise = per-measurement Allan SD * 3 * t_crit
#                              (t_crit = qt(.995/.975/.95, df = n_meas_pts - 2)).
# A flux is "below MDF" when abs(best.flux) < MDF. Also reports noise floor & SNR.
#
# Can be run with Rscript (non-interactive).
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/transect-trees-2026/scripts/00_setup.R"))

load(file.path(rdata_dir, "flux_results.RData"))   # CO2_best, CH4_best
load(file.path(rdata_dir, "manID.RData"))          # manID (flagged windows)

# --- Manufacturer 1-sigma precision (LI-7810; from ch4-data-filtering) --------
prec_co2 <- li7810_prec[1]   # 3.5 ppm
prec_ch4 <- li7810_prec[2]   # 0.6 ppb

# --- Allan deviation: white-noise floor from a trending signal ----------------
# allan_sd(x) = sd(diff(x)) / sqrt(2)  (first difference removes the flux trend)
allan_sd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3) return(NA_real_)
  stats::sd(diff(x)) / sqrt(2)
}

# Per-measurement Allan deviation on the identified (flag == 1) points
allan <- manID %>%
  filter(flag == 1) %>%
  group_by(UniqueID) %>%
  summarise(
    allan_sd_CO2 = allan_sd(CO2dry_ppm),
    allan_sd_CH4 = allan_sd(CH4dry_ppb),
    n_meas_pts   = sum(!is.na(CH4dry_ppb)),
    .groups = "drop"
  ) %>%
  mutate(UniqueID = trimws(UniqueID))

# --- Pull flux.term / nb.obs / best.flux / native MDF from goFlux results ------
g_cols <- function(df, gas) {
  df %>%
    mutate(UniqueID = trimws(UniqueID)) %>%
    transmute(
      UniqueID,
      best.flux = best.flux,
      flux.term = flux.term,
      nb.obs    = nb.obs,
      MDF_native= MDF
    ) %>%
    rename_with(~ paste0(gas, "_", .x), -UniqueID)
}

mdf <- allan %>%
  full_join(g_cols(CO2_best, "CO2"), by = "UniqueID") %>%
  full_join(g_cols(CH4_best, "CH4"), by = "UniqueID")

# t = measurement length in seconds (1 Hz => nb.obs == seconds)
mdf <- mdf %>% mutate(t_sec = CO2_nb.obs)

# --- Global empirical SD (Wassmann): median per-measurement Allan dev ---------
gsd_co2 <- median(mdf$allan_sd_CO2, na.rm = TRUE)
gsd_ch4 <- median(mdf$allan_sd_CH4, na.rm = TRUE)
message(sprintf("Global Allan SD (Wassmann): CO2 = %.4g ppm, CH4 = %.4g ppb",
                gsd_co2, gsd_ch4))

# --- Confidence multipliers ---------------------------------------------------
z99 <- qnorm(0.995); z95 <- qnorm(0.975); z90 <- qnorm(0.95)

# --- Compute the three approaches x three confidence levels -------------------
mdf <- mdf %>%
  mutate(
    df_meas = pmax(n_meas_pts - 2, 1),
    t99 = qt(0.995, df_meas), t95 = qt(0.975, df_meas), t90 = qt(0.95, df_meas),

    # 1. Manufacturer / goFlux (no confidence multiplier; == native MDF)
    CO2_MDF_manuf = prec_co2 / t_sec * CO2_flux.term,
    CH4_MDF_manuf = prec_ch4 / t_sec * CH4_flux.term,

    # 2. Wassmann 2018: z * global empirical SD
    CO2_MDF_wass99 = z99 * gsd_co2 / t_sec * CO2_flux.term,
    CO2_MDF_wass95 = z95 * gsd_co2 / t_sec * CO2_flux.term,
    CO2_MDF_wass90 = z90 * gsd_co2 / t_sec * CO2_flux.term,
    CH4_MDF_wass99 = z99 * gsd_ch4 / t_sec * CH4_flux.term,
    CH4_MDF_wass95 = z95 * gsd_ch4 / t_sec * CH4_flux.term,
    CH4_MDF_wass90 = z90 * gsd_ch4 / t_sec * CH4_flux.term,

    # 3. Christiansen 2015: per-measurement Allan SD * 3 * t_crit
    CO2_MDF_chr99 = allan_sd_CO2 * 3 * t99 / t_sec * CO2_flux.term,
    CO2_MDF_chr95 = allan_sd_CO2 * 3 * t95 / t_sec * CO2_flux.term,
    CO2_MDF_chr90 = allan_sd_CO2 * 3 * t90 / t_sec * CO2_flux.term,
    CH4_MDF_chr99 = allan_sd_CH4 * 3 * t99 / t_sec * CH4_flux.term,
    CH4_MDF_chr95 = allan_sd_CH4 * 3 * t95 / t_sec * CH4_flux.term,
    CH4_MDF_chr90 = allan_sd_CH4 * 3 * t90 / t_sec * CH4_flux.term,

    # Empirical noise floor (1-sigma in flux units) and signal-to-noise ratio
    CO2_noise_floor = allan_sd_CO2 / t_sec * CO2_flux.term,
    CH4_noise_floor = allan_sd_CH4 / t_sec * CH4_flux.term,
    CO2_SNR = abs(CO2_best.flux) / CO2_noise_floor,
    CH4_SNR = abs(CH4_best.flux) / CH4_noise_floor,

    # "Below MDF" flags (abs(flux) < MDF) at each method's 99% level
    CO2_below_MDF_manuf = abs(CO2_best.flux) < CO2_MDF_manuf,
    CO2_below_MDF_wass99= abs(CO2_best.flux) < CO2_MDF_wass99,
    CO2_below_MDF_chr99 = abs(CO2_best.flux) < CO2_MDF_chr99,
    CH4_below_MDF_manuf = abs(CH4_best.flux) < CH4_MDF_manuf,
    CH4_below_MDF_wass99= abs(CH4_best.flux) < CH4_MDF_wass99,
    CH4_below_MDF_chr99 = abs(CH4_best.flux) < CH4_MDF_chr99
  )

# --- Report: % below MDF per method (99% level) ------------------------------
pct_below <- function(x) {
  ok <- !is.na(x); sprintf("%d / %d (%.0f%%)", sum(x[ok]), sum(ok), 100*mean(x[ok]))
}
message("\n=== Fraction below detection (|flux| < MDF, 99%) ===")
message("CO2 manuf/goFlux : ", pct_below(mdf$CO2_below_MDF_manuf))
message("CO2 Wassmann     : ", pct_below(mdf$CO2_below_MDF_wass99))
message("CO2 Christiansen : ", pct_below(mdf$CO2_below_MDF_chr99))
message("CH4 manuf/goFlux : ", pct_below(mdf$CH4_below_MDF_manuf))
message("CH4 Wassmann     : ", pct_below(mdf$CH4_below_MDF_wass99))
message("CH4 Christiansen : ", pct_below(mdf$CH4_below_MDF_chr99))

# --- Merge onto the compiled flux table and save -----------------------------
compiled_path <- file.path(results_dir, "transect_trees_2026_fluxes.csv")
if (file.exists(compiled_path)) {
  compiled <- readr::read_csv(compiled_path, show_col_types = FALSE) %>%
    mutate(UniqueID = trimws(UniqueID))
  out <- left_join(compiled, mdf %>% select(-CO2_best.flux, -CH4_best.flux,
                                            -CO2_flux.term, -CH4_flux.term,
                                            -CO2_nb.obs, -CH4_nb.obs),
                   by = "UniqueID")
  write.csv(out, file.path(results_dir, "transect_trees_2026_fluxes_with_mdf.csv"),
            row.names = FALSE)
  write.xlsx(out, file.path(results_dir, "transect_trees_2026_fluxes_with_mdf.xlsx"))
  message("\nSaved: ", file.path(results_dir, "transect_trees_2026_fluxes_with_mdf.csv"))
} else {
  message("\n(Run 05_compile_results.R first to produce the compiled table to merge onto.)")
}

save(mdf, gsd_co2, gsd_ch4, file = file.path(rdata_dir, "mdf_results.RData"))
write.csv(mdf, file.path(results_dir, "mdf_lod_table.csv"), row.names = FALSE)
message("Saved MDF table: ", file.path(results_dir, "mdf_lod_table.csv"))
