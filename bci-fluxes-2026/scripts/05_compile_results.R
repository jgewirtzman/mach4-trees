# =============================================================================
# 05_compile_results.R  (BCI)
# Merge CO2 & CH4 goFlux results with field metadata into one tidy table.
# Identical to transect-trees-2026/scripts/05_compile_results.R except output
# file names (bci_2026_fluxes.*).
# Can be run with Rscript (non-interactive).
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/bci-fluxes-2026/scripts/00_setup.R"))

load(file.path(rdata_dir, "flux_results.RData"))   # CO2_best, CH4_best
load(file.path(rdata_dir, "auxfile.RData"))        # field_meta

pick <- function(df, gas) {
  df %>%
    mutate(UniqueID = trimws(UniqueID)) %>%
    select(
      UniqueID,
      best.flux, model, quality.check,
      LM.flux, LM.r2, LM.p.val, LM.RMSE,
      HM.flux, HM.r2, HM.RMSE,
      LM.diagnose, HM.diagnose,
      nb.obs, flux.term, MDF, g.fact
    ) %>%
    rename_with(~ paste0(gas, "_", .x), -UniqueID)
}

flux <- full_join(pick(CO2_best, "CO2"), pick(CH4_best, "CH4"), by = "UniqueID")

final <- field_meta %>%
  mutate(UniqueID = trimws(UniqueID)) %>%
  left_join(flux, by = "UniqueID") %>%
  arrange(TreeID, desc(height_cm)) %>%
  mutate(Area_cm2 = chamber_area_cm2, Vtot_L = Vtot_L,
         Tcham_C = Tcham_default_C, Pcham_kPa = Pcham_default_kPa)

write.csv(final, file.path(results_dir, "bci_2026_fluxes.csv"), row.names = FALSE)
write.xlsx(final, file.path(results_dir, "bci_2026_fluxes.xlsx"))

message("=== Compiled: ", nrow(final), " measurements ===")
message("With CO2 flux: ", sum(!is.na(final$CO2_best.flux)), " / ", nrow(final))
message("With CH4 flux: ", sum(!is.na(final$CH4_best.flux)), " / ", nrow(final))
message("Saved: ", file.path(results_dir, "bci_2026_fluxes.csv"))
message("Proceed to 06_mdf_lod.R for the MDF/LOD analysis.")
