# =============================================================================
# 04_flux_calculation.R  (BCI)
# Calculate CO2 & CH4 fluxes with goFlux() and pick the best model with
# best.flux(). Identical to transect-trees-2026/scripts/04_flux_calculation.R.
#
# Units:  CO2 flux = umol m-2 s-1   |   CH4 flux = nmol m-2 s-1
# NOTE: requires the real class B chamber geometry in 00_setup.R (placeholder
# geometry gives placeholder flux magnitudes).
# Can be run with Rscript (non-interactive).
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/bci-fluxes-2026/scripts/00_setup.R"))

if (isTRUE(CHAMBER_GEOMETRY_IS_PLACEHOLDER))
  warning("CHAMBER GEOMETRY IS A PLACEHOLDER in 00_setup.R -- flux magnitudes ",
          "are NOT final until you set the real class B dimensions.")

load(file.path(rdata_dir, "manID.RData"))   # manID (from interactive 03)

message("=== Calculating CO2 fluxes ===")
CO2_flux <- goFlux(manID, "CO2dry_ppm", H2O_col = "H2O_ppm")

message("=== Calculating CH4 fluxes ===")
CH4_flux <- goFlux(manID, "CH4dry_ppb", H2O_col = "H2O_ppm")

message("Selecting best model (criteria: ", paste(flux_criteria, collapse = ", "), ")")
CO2_best <- best.flux(CO2_flux, flux_criteria)
CH4_best <- best.flux(CH4_flux, flux_criteria)

message("CO2: ", sum(CO2_best$model == "LM"), " LM / ",
        sum(CO2_best$model == "HM"), " HM")
message("CH4: ", sum(CH4_best$model == "LM"), " LM / ",
        sum(CH4_best$model == "HM"), " HM")

save(CO2_best, CH4_best, CO2_flux, CH4_flux,
     file = file.path(rdata_dir, "flux_results.RData"))
write.xlsx(CO2_best, file.path(results_dir, "CO2_best.xlsx"))
write.xlsx(CH4_best, file.path(results_dir, "CH4_best.xlsx"))

message("Saved flux results. Proceed to 05_compile_results.R")
