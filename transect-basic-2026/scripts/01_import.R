# =============================================================================
# 01_import.R  (BASIC trees)
# Subset the LI-7810 .data file to the campaign dates, then import with
# goFlux::import.LI7810(). Identical logic to transect-trees-2026/01_import.R.
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/transect-basic-2026/scripts/00_setup.R"))

if (!file.exists(licor_subset)) {
  message("Subsetting LI-7810 file to: ", paste(measurement_dates, collapse = ", "))
  date_alt <- paste(measurement_dates, collapse = "|")
  awk_prog <- sprintf(
    'NR<=7 {print; next} $1=="DATA" && $7 ~ /^(%s)$/ {print}', date_alt)
  cmd <- sprintf("awk -F'\\t' %s %s > %s",
                 shQuote(awk_prog), shQuote(licor_raw), shQuote(licor_subset))
  if (system(cmd) != 0) stop("awk subsetting failed")
  message("Wrote subset: ", licor_subset)
} else message("Using existing subset: ", licor_subset)

message("Importing with import.LI7810() ...")
licor <- import.LI7810(inputfile = licor_subset, date.format = "ymd",
                       timezone = tz_data, prec = li7810_prec)

message("Imported ", nrow(licor), " rows; time range: ",
        format(min(licor$POSIX.time)), " to ", format(max(licor$POSIX.time)),
        " (", tz_data, ")")
for (g in c("CO2dry_ppm", "CH4dry_ppb", "H2O_ppm")) {
  v <- licor[[g]]
  message(sprintf("  %-11s n=%d  range %.1f .. %.1f", g, sum(!is.na(v)),
                  min(v, na.rm = TRUE), max(v, na.rm = TRUE)))
}
save(licor, file = file.path(rdata_dir, "licor_imported.RData"))
message("Saved licor_imported.RData. Proceed to 02_build_auxfile.R")
