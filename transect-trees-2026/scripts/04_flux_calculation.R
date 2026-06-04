# =============================================================================
# 04_flux_calculation.R
# Calculate CO2 & CH4 fluxes with goFlux() and pick the best model with
# best.flux().  Mirrors whole_tree_flux ymf_05_flux_calculation.R.
#
# Units:  CO2 flux = umol m-2 s-1   |   CH4 flux = nmol m-2 s-1
# goFlux() returns MDF and g.fact per measurement; best.flux() uses flux_criteria.
# Can be run with Rscript (non-interactive).
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/transect-trees-2026/scripts/00_setup.R"))

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
